package daxing.applets;

import daxing.common.VCF;
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
import java.util.stream.Collectors;

public class Vmap2QC {

    public Vmap2QC(){
        this.assessFalesPositiveRate("/data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII");
//        this.getMaf("vcfDir", "outFileDir");
//        this.mergeAddRefCSFordxy("chrAB.vcf",
//                "chrD.vcf",
//                "merged.vcf");
    }

    public void assessFalesPositiveRate(String vcfDir){
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
                int count= Vmap2QC.assessFalesPositiveRate(f1[index]);
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
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, 3);
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
                        infor=PStringUtils.fastSplit(lineList.get(7), ";");
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

    /**
     * 以merged后两个vcf为输入文件（chrAB.vcf, chrD.vcf），merge为一个vcf文件(taxon依次为AABBDD, AABB, DD), 并在taxons最后添加refCS
     * @param abVcfInputFile chrAB.vcf
     * @param dVcfInput chrD.vcf
     * @param mergedVCF 包含refCS的基因型, 0/0
     */
    public void mergeAddRefCSFordxy(String abVcfInputFile, String dVcfInput, String mergedVCF){
        try(BufferedReader br1=IOUtils.getTextReader(abVcfInputFile);
            BufferedReader br2=IOUtils.getTextReader(dVcfInput);
            BufferedWriter bw=IOUtils.getTextWriter(mergedVCF)){
            String line;
            List<String> lineList;
            String bread;
            String emmerTaxons;
            String tauschiiTaxons;
            while ((line=br1.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            lineList= PStringUtils.fastSplit(line);
            emmerTaxons=lineList.stream().skip(9+419).collect(Collectors.joining("\t"));
            while ((line=br2.readLine()).startsWith("##")){}
            lineList= PStringUtils.fastSplit(line);
            tauschiiTaxons=lineList.stream().skip(9+419).collect(Collectors.joining("\t"));
            bread=lineList.stream().limit(9+419).collect(Collectors.joining("\t"));
            bw.write(bread+"\t"+emmerTaxons+"\t"+tauschiiTaxons+"\tRefCS");
            bw.newLine();
            String emmerGenotype;
            String tauschiiGenotype;
            String refCSGenotype="0/0:50,50:0,7,44";
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < 187; i++) {
                sb.append("./.").append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            emmerGenotype=sb.toString();
            sb=new StringBuilder();
            for (int i = 0; i < 36; i++) {
                sb.append("./.").append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            tauschiiGenotype=sb.toString();
            while ((line=br1.readLine())!=null){
                sb=new StringBuilder();
                sb.append(line).append("\t").append(tauschiiGenotype).append("\t").append(refCSGenotype);
                bw.write(sb.toString());
                bw.newLine();
            }
            while ((line=br2.readLine())!=null){
                sb=new StringBuilder();
                lineList=PStringUtils.fastSplit(line);
                sb.append(lineList.stream().limit(9+419).collect(Collectors.joining("\t"))).append("\t");
                sb.append(emmerGenotype).append("\t");
                sb.append(lineList.stream().skip(9+419).collect(Collectors.joining("\t"))).append("\t");
                sb.append(refCSGenotype);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
        this.sort(mergedVCF);
    }

    private void sort(String mergedVCF){
        VCF vcf=new VCF(mergedVCF);
        vcf.sort();
        vcf.write(mergedVCF);
    }

//    public static void main(String[] args) {
//        new Vmap2QC();
//    }

}
