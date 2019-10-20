package daxing.applets;

import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Predicate;

public class DeleteriousMutation {

    String vcfDir;
    public DeleteriousMutation(String vcfDir){
        this.vcfDir=vcfDir;
        this.assessFalesPositiveRate();
//        this.getMaf(vcfDir);
    }

    public void assessFalesPositiveRate(){
        double wheatSize=14066280851d;
        AtomicInteger total= new AtomicInteger(0);
        File[] files=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> p=File::isHidden;
        File[] f1= Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, 5);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(index->{
                int count=DeleteriousMutation.assessFalesPositiveRate(f1[index]);
                total.addAndGet(count);
            });
        }
        System.out.println("Fales Positive Rate: "+ total.get()/wheatSize+"("+total.get()+"/"+14066280851d+")");
    }

    private static int assessFalesPositiveRate(File vcfFile){
        int count=0;
        try (BufferedReader br = IOUtils.getTextReader(vcfFile.getAbsolutePath())) {
            String line;
            List<String> lineList;
            while ((line=br.readLine()).startsWith("##")){}
            lineList= PStringUtils.fastSplit(line);
            int index=lineList.indexOf("CS");
            List<String> infor;
            String genotype;
            while ((line=br.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                genotype=lineList.get(index).substring(0,3);
                if(genotype.equals("./.")) continue;
                if (!genotype.equals("0/0")){
                    count++;
                }
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return count;
    }

    public void getMaf(String vcfDir, String outFileDir){
        File[] files=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> p=File::isHidden;
        File[] f1= Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        String[] outNames=Arrays.stream(f1).map(File::getName).map(str->str.replaceAll("vcf$", "maf.txt"))
                .toArray(String[]::new);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, 5);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(index->{
                try (BufferedReader br =IOUtils.getTextReader(f1[index].getAbsolutePath());
                     BufferedWriter bw=IOUtils.getTextWriter(new File(outFileDir, outNames[index]).getAbsolutePath())) {
                    String line;
                    List<String> lineList;
                    List<String> infor;
                    double maf=-1d;
                    double maf_ABD=-1d;
                    double maf_ABorD=-1d;
                    StringBuilder sb;
                    bw.write("MAF\tMAF_ABD\tMAF_ABorD");
                    bw.newLine();
                    while ((line=br.readLine()).startsWith("##")){}
                    while ((line=br.readLine())!=null){
                        lineList=PStringUtils.fastSplit(line);
                        infor=PStringUtils.fastSplit(lineList.get(7));
                        maf=Double.parseDouble(PStringUtils.fastSplit(infor.get(6), "=").get(1));
                        maf_ABD=Double.parseDouble(PStringUtils.fastSplit(infor.get(7), "=").get(1));
                        maf_ABorD=Double.parseDouble(PStringUtils.fastSplit(infor.get(8), "=").get(1));
                        sb=new StringBuilder();
                        sb.append(maf).append("\t").append(maf_ABD).append("\t").append(maf_ABorD);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    bw.flush();
                }catch (Exception e){
                    e.printStackTrace();
                }
            });
        }
    }



//    public static void main(String[] args) {
//        new DeleteriousMutation(args[0]);
//    }

}
