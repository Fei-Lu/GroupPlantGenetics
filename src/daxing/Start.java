package daxing;

import daxing.common.ChrRange;
import daxing.common.ChrRanges;
import daxing.common.IOTool;
import gnu.trove.set.hash.TIntHashSet;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) {
        //        String rangeFile="/Users/xudaxing/Desktop/temp/002_DE_FTT_Range/DE_FTT_IntrogressionRange.txt.gz";
//        String snpAnnpDBDir="/Users/xudaxing/Desktop/temp/003_outDir";
//        String outDir="/Users/xudaxing/Desktop/temp/004_outDir";
//        ChrSNPAnnoDB.filterRange(rangeFile, snpAnnpDBDir, outDir);

//        File file =new File("/Users/xudaxing/Desktop/chr7A_vmap2.1_LR_0186_DE_0035.csv");
//        System.out.println(file.length());

        List<ChrRange> chrRanges=new ArrayList<>();
        chrRanges.add(new ChrRange("2B", 1, 12));
        chrRanges.add(new ChrRange("2B", 3, 12));
        chrRanges.add(new ChrRange("2B", 5, 7));
        chrRanges.add(new ChrRange("2B", 8, 9));
        ChrRanges chrRanges1= new ChrRanges(chrRanges);
        chrRanges1.sortBy(ChrRanges.SortType.POSITION);
        int[] index = chrRanges1.getIndicesContainsPosition("2B", 3);
        System.out.println();
    }

    public static void getMissing(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, ".vcf.gz");
        String[] outNames =files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            GenotypeGrid grid = new GenotypeGrid(files.get(e).getAbsolutePath(), GenoIOFormat.VCF_GZ);
            int chrID= Integer.parseInt(files.get(e).getName().substring(3,6));
            try (BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                StringBuilder sb = new StringBuilder();
                String[] taxaNames= grid.getTaxaNames();
                NumberFormat numberFormat = NumberFormat.getInstance();
                numberFormat.setGroupingUsed(false);
                numberFormat.setMaximumFractionDigits(5);
                bw.write("Taxa\tChrID\tMissingCount\tTotalCount");
                bw.newLine();
                double missingRate, missingNum;
                int totalCount;
                for (int i = 0; i <grid.getTaxaNumber(); i++) {
                    missingNum=grid.getMissingNumberByTaxon(i);
//                    missingRate=missingNum/grid.getSiteNumber();
                    totalCount=grid.getSiteNumber();
                    sb.setLength(0);
                    sb.append(taxaNames[i]).append("\t").append(chrID).append("\t").append(missingNum).append("\t");
                    sb.append(totalCount);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }

        });
    }
}