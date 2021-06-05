package daxing.hybrid.dection;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

public class HybridDection {

    public static void getDepthHeterozygosity(String inputDir, String outDir_depth, String outDir_heterozygosity){
        List<File> fileList= IOTool.getVisibleFileListInDir(inputDir);
        String[] outFileName= fileList.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", ".txt.gz")).toArray(String[]::new);
//        IntStream.range(0, fileList.size()).parallel().forEach(e->getDepth(fileList.get(e), new File(outDir_depth, outFileName[e])));
        IntStream.range(0, fileList.size()).parallel().forEach(e->getHeterozygosity(fileList.get(e), new File(outDir_heterozygosity,
                outFileName[e])));

    }

    private static void getDepth(File inputFile, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            List<String> temp, tem, te;
            List<String> taxaList;
            String line;
            int depth;
            while ((line=br.readLine()).startsWith("##")) continue;
            temp=PStringUtils.fastSplit(line);
            taxaList=temp.subList(9,temp.size());
            StringBuilder sb=new StringBuilder();
            sb.append("Chr").append("\t").append("Pos").append("\t");
            sb.append(String.join("\t", taxaList));
            bw.write(sb.toString());
            bw.newLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t");
                for (int i = 9; i < temp.size(); i++) {
                    if (temp.get(i).startsWith("./.")){
                        sb.append("NA").append("\t");
                        continue;
                    }
                    tem=PStringUtils.fastSplit(temp.get(i), ":");
                    te=PStringUtils.fastSplit(tem.get(1), ",");
                    depth=Integer.parseInt(te.get(0))+Integer.parseInt(te.get(1));
                    sb.append(depth).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void merge(String inputDir, String outFile){
        List<File> files= IOTool.getVisibleFileListInDir(inputDir);
        RowTableTool<String> table0=new RowTableTool<>(files.get(0).getAbsolutePath());
        RowTableTool<String> table;
        for (int i = 1; i < files.size(); i++) {
            table= new RowTableTool<>(files.get(i).getAbsolutePath());
            table0.add(table);
        }
        table0.write(outFile, IOFileFormat.TextGzip);
    }

    public static void getHeterozygosity(File inputFile, File outFile){
        GenotypeGrid genotypeGrid=new GenotypeGrid(inputFile.getAbsolutePath(), GenoIOFormat.VCF_GZ);
        int taxonNum=genotypeGrid.getTaxaNumber();
        short chr= genotypeGrid.getChromosome(0);
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            sb.append("Taxon").append("\t").append("Chr").append("\t").append("HeterozygoteNumber").append("\t");
            sb.append("HomozygoteNumber");
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i <taxonNum; i++) {
                sb.setLength(0);
                sb.append(genotypeGrid.getTaxonName(i)).append("\t").append(chr).append("\t");
                sb.append(genotypeGrid.getHeterozygoteNumberByTaxon(i)).append("\t");
                sb.append(genotypeGrid.getHomozygoteNumberByTaxon(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
