package daxing;

import daxing.common.IOTool;
import daxing.common.MD5;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) {
//        String exonSNPAnnoDir=args[0];
//        String exonVCFDir=args[1];
//        String taxa_InfoDB=args[2];
//        String triadFile=args[3];
//        String pgfFile=args[4];
//        String outDir =args[5];
//        int blockGeneNum=Integer.parseInt(args[6]);

        String file1="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/009_SIFT_GERP_NEW/C1.txt";
        String file2="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/008_SIFT_GERP/002_countMerge/C1.txt";
        MD5.checkTwoFileMD5(file1, file2);

//        String exonSNPAnnoDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_exon/002_exonVCFAnno/002_exonAnnotationByDerivedSift/001_exonSNPAnnotation";
//        String exonVCFDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_exon/001_exonVCF";
//        String taxa_InfoDB="/Users/xudaxing/Documents/deleteriousMutation/002_vmapII_taxaGroup/taxa_InfoDB.txt";
//        String triadFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/triadGenes1.1.txt";
//        String pgfFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/wheat_v1.1_Lulab.pgf";
//        String outDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/009_SIFT_GERP_NEW";
//        int blockGeneNum=20;
//        File subDir ;
////        for (int i = 0; i < SNPAnnotation.MethodCallDeleterious.values().length; i++) {
////            subDir= new File(outDir, PStringUtils.getNDigitNumber(3, i+1)+SNPAnnotation.MethodCallDeleterious.values()[i]);
////            subDir.mkdir();
////            ComplementaryGo.go(exonSNPAnnoDir, exonVCFDir, taxa_InfoDB, triadFile, blockGeneNum, pgfFile,
////                    subDir.getAbsolutePath(),
////                    SNPAnnotation.MethodCallDeleterious.values()[i]);
////        }
//        ComplementaryGo.go(exonSNPAnnoDir, exonVCFDir,taxa_InfoDB, triadFile, 20, pgfFile, outDir,
//                SNPAnnotation.MethodCallDeleterious.SIFT_GERP);

    }

    public static void biAllele(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, "vcf.gz");
        String[] outNames=
                files.stream().map(File::getName).map(s -> s.replaceAll("2.0.vcf.gz","2.0.vcf")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line, subLine;
                List<String> temp;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write("#CHROM\t"+line.substring(5));
                bw.newLine();
                int cnt=0;
                Set<String> altSet=new HashSet<>();
                while ((line=br.readLine())!=null){
                    subLine=line.substring(0,40);
                    temp= PStringUtils.fastSplit(subLine);
                    if (temp.get(4).length() > 1){
                        altSet.add(temp.get(4));
                    }
                    if (temp.get(4).equals("A") || temp.get(4).equals("T") ||
                            temp.get(4).equals("C") || temp.get(4).equals("G") ||
                            temp.get(4).equals("I") || temp.get(4).equals("D")){
                        bw.write(line);
                        bw.newLine();
                        cnt++;
                    }else {
                        altSet.add(temp.get(4));
                    }
                }
                bw.flush();
//                altSet.stream().forEach(System.out::println);
                System.out.println(new File(outDir, outNames[e]).getName()+  " have "+cnt+" biVariants");
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

}