package daxing.vmapII_1000.variantsAnnotation;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import daxing.common.IOTool;
import daxing.common.PGF;
import daxing.common.RowTableTool;
import pgl.infra.range.Range;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.IntStream;

/**
 * using longest transcript
 */
public class GeneSiteAnnoDB {

    public static void extractGeneSiteInfoForAnnoDB(String vcfInputDir, String pgfFile, String nonOverlapGeneFile,
                                                    String ancestralDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(vcfInputDir, ".vcf.gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".siteAnno.txt.gz")).toArray(String[]::new);
        List<File> ancestralFiles= IOTool.getFileListInDirEndsWith(ancestralDir, ".gz");
        RowTableTool<String> nonOverlapGeneTable = new RowTableTool<>(nonOverlapGeneFile);
        Predicate<List<String>> uniqueGene = l -> l.get(3).equals("1");
        nonOverlapGeneTable.removeIf(uniqueGene.negate());
        Set<String> uniqueGeneSet=new HashSet<>(nonOverlapGeneTable.getColumn(0));
        PGF pgf = new PGF(pgfFile);
        Predicate<PGF.Gene> uniqueGeneP= gene -> uniqueGeneSet.contains(gene.getGeneName());
        pgf.removeIf(uniqueGeneP.negate());
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedReader brAnc = IOTool.getReader(ancestralFiles.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                List<String> temp;
                int chrID, pos;
                Range range;
                brAnc.readLine();
                Table<Integer,Integer, String> chrPosAncTable = HashBasedTable.create();
                while ((line=brAnc.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    chrID = Integer.parseInt(temp.get(0));
                    pos = Integer.parseInt(temp.get(1));
                    int geneIndex = pgf.getGeneIndex(chrID, pos);
                    if (geneIndex < 0) continue;
                    int longestTranscriptIndex = pgf.getLongestTranscriptIndex(geneIndex);
                    range = pgf.getTranscriptRange(geneIndex, longestTranscriptIndex);
                    if (!range.isContain(chrID, pos)) continue;
                    chrPosAncTable.put(chrID, pos, temp.get(2));
                }
                double refCount, altCount, heteroCount;
                double maf, altFre;
                String daf;
                String refAllele, altAllele, majorAllele, minorAllele, longestTranscript, ancestralAllele;
                StringBuilder sb = new StringBuilder();
                bw.write("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tTranscript\tAncestral\tDAF");
                bw.newLine();
                int total=0;
                while ((line=br.readLine()).startsWith("##")){}
                NumberFormat numberFormat = NumberFormat.getInstance();
                numberFormat.setMaximumFractionDigits(5);
                numberFormat.setGroupingUsed(false);
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    chrID = Integer.parseInt(temp.get(0));
                    pos = Integer.parseInt(temp.get(1));
                    int geneIndex = pgf.getGeneIndex(chrID, pos);
                    if (geneIndex < 0) continue;
                    int longestTranscriptIndex = pgf.getGene(geneIndex).getLongestTranscriptIndex();
                    range = pgf.getTranscriptRange(geneIndex, longestTranscriptIndex);
                    if (!range.isContain(chrID, pos)) continue;
                    longestTranscript = pgf.getTranscriptName(geneIndex, longestTranscriptIndex);
                    total++;
                    refCount =0 ;
                    altCount =0;
                    heteroCount = 0;
                    for (int i = 9; i < temp.size(); i++) {
                        if (temp.get(i).startsWith("./.")) continue;
                        if (temp.get(i).startsWith("0/0")){
                            refCount++;
                        }else if (temp.get(i).startsWith("0/1")){
                            heteroCount++;
                        }else if (temp.get(i).startsWith("1/1")){
                            altCount++;
                        }
                    }
                    refAllele = temp.get(3);
                    altAllele = temp.get(4);
                    altFre = (altCount*2+heteroCount)/((refCount+altCount+heteroCount)*2);
                    maf = altFre > 0.5 ? 1- altFre : altFre;
                    majorAllele = altFre > 0.5 ? altAllele : refAllele;
                    minorAllele = altFre <= 0.5 ? altAllele : refAllele;
                    ancestralAllele = chrPosAncTable.get(chrID, pos)==null ? "NA" : chrPosAncTable.get(chrID, pos);
                    daf = majorAllele.equals(ancestralAllele) ? numberFormat.format(maf) :
                            minorAllele.equals(ancestralAllele) ?
                            numberFormat.format(1-maf) : "NA";
                    sb.setLength(0);
                    sb.append(temp.get(2)).append("\t").append(chrID).append("\t").append(pos).append("\t");
                    sb.append(temp.get(3)).append("\t").append(temp.get(4)).append("\t");
                    sb.append(majorAllele).append("\t").append(minorAllele).append("\t");
                    sb.append(numberFormat.format(maf)).append("\t").append(longestTranscript).append("\t");
                    sb.append(ancestralAllele).append("\t").append(daf);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                System.out.println(outNames[e]+" have "+ total+ " snp in gene region (only longest transcript)");
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }
}
