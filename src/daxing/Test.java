package daxing;

import daxing.common.ArrayTool;
import daxing.common.CollectionTool;
import daxing.common.RowTableTool;
import daxing.common.VCF;
import utils.IOUtils;
import utils.PStringUtils;
import utils.Tuple;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Test {

    public static Tuple<int[], int[]> readCgLog(String callGenotypeLogFile){
        Tuple<int[], int[]> t=null;
        int[] chrNum;
        int[] snpNum;
        try(BufferedReader bw= IOUtils.getTextReader(callGenotypeLogFile)){
            chrNum=bw.lines().filter(e->e.startsWith("Genotyping")).skip(1).limit(6)
                     .map(e->PStringUtils.fastSplit(e, " ").get(3)).mapToInt(Integer::valueOf).toArray();
            snpNum= Files.newBufferedReader(Paths.get(callGenotypeLogFile)).lines()
                         .filter(e->e.startsWith("Genotyping")).skip(1).limit(6)
                         .map(e->PStringUtils.fastSplit(e, " ").get(5)).mapToInt(Integer::valueOf).toArray();
            t=new Tuple(chrNum, snpNum);
            return t;
        }catch (Exception e){
            e.printStackTrace();
        }
        return t;
    }

    public static void calcalateGenotypedTaxonNum(String vcfInputDir, String outputDir, int[] chr){
        List<String> allVcf= VCF.getAllVcfInputPath(vcfInputDir,chr);
        Map<String, BufferedWriter> map=new HashMap<>();
        String key;
        BufferedWriter value;
        for (int i = 0; i <chr.length ; i++) {
            key=allVcf.get(i);
            value=IOUtils.getTextWriter(outputDir+"/chr"+PStringUtils.getNDigitNumber(3, chr[i])+".txt");
            map.put(key, value);
        }
        List<Map<Integer, Long>> l=allVcf.stream().map(VCF::new).map(VCF::getGenotypedTaxonNum)
                                         .map(ArrayTool::calcalateRepetitionNum).collect(Collectors.toList());
        IntStream.range(0,chr.length).parallel().forEach(e->{
            BufferedWriter bw=map.get(allVcf.get(e));
            Map<Integer, Long> m=l.get(e);
            TreeMap<Integer, Long> t=new TreeMap<>(m);
            try {
                bw.write("key"+"\t"+"value"+"\n");
                StringBuilder sb=new StringBuilder();
                for (Map.Entry<Integer, Long> ele: t.entrySet()){
                    sb.append(ele.getKey()).append("\t").append(ele.getValue());
                    bw.write(sb.toString());
                    bw.newLine();
                    sb=new StringBuilder();
                }
                bw.flush();
                bw.close();
            }catch (Exception ex){
                ex.printStackTrace();
            }
        });
    }

    public static void getRes(String inputDir){
        File[] files1=IOUtils.listRecursiveFiles(new File(inputDir));
        File[] files=IOUtils.listFilesEndsWith(files1, "txt");
        List<Map<Integer, Integer>> l= Arrays.stream(files).map(File::getAbsolutePath).map(RowTableTool<String>::new)
                .map(e->e.getMap(1)).map(CollectionTool::changeToIntMap).collect(Collectors.toList());
        Map<Integer, Integer> res=new HashMap<>();
        res.putAll(l.get(0));
        Integer key, value, tempValue;
        for (int i = 1; i < l.size(); i++) {
            for (Map.Entry<Integer, Integer> entry:l.get(i).entrySet()){
                key=entry.getKey();
                value=entry.getValue();
                if (res.containsKey(key)){
                    tempValue=res.get(key);
                    res.remove(key);
                    res.put(key,value+tempValue);
                }else {
                    res.put(key, value);
                }
            }
        }
        TreeMap<Integer, Integer> treeMap=new TreeMap<>(res);
//        TDoubleArrayList valuesList=new TDoubleArrayList();
////        List<String> namesList=new ArrayList<>();
////        for (Map.Entry<Integer, Integer> entry:treeMap.entrySet()){
////            Integer k=entry.getKey();
////            Integer v=entry.getValue();
////            valuesList.add(v);
////            namesList.add(String.valueOf(k));
////        }
////        double[] values=valuesList.toArray();
////        String[] names=namesList.toArray(new String[namesList.size()]);
////        BarPlot barPlot=new BarPlot(values, names);
////        barPlot.showGraph();
        try(BufferedWriter bw=IOUtils.getTextWriter(inputDir+"/all.txt")){
            bw.write("key"+"\t"+"value");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            for (Map.Entry<Integer, Integer> entry: treeMap.entrySet()){
                sb.append(entry.getKey()).append("\t").append(entry.getValue());
                bw.write(sb.toString());
                bw.newLine();
                sb=new StringBuilder();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }


//    public static void main(String[] args) {
//        File[] f=IOUtils.listRecursiveFiles(new File(args[0]));
//        File[] files=IOUtils.listFilesEndsWith(f, ".vcf");
//        VCF vcf1=new VCF(files[0]);
//        vcf1.removeVarianceWithHighGenotypeMissingRate(0.9);
//        System.out.println(files[0].getName()+"\t"+ vcf1.getData().size());
//        VCF vcf;
//        for (int i = 1; i <files.length ; i++) {
//            vcf=new VCF(files[i]);
//            vcf.removeVarianceWithHighGenotypeMissingRate(0.9);
//            System.out.println(files[i].getName()+"\t"+vcf.getData().size());
//            vcf1.addVCF(vcf);
//        }
//        vcf1.write(args[1]);
//    }
}
