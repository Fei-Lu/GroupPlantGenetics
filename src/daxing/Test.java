package daxing;

import daxing.common.ArrayTool;
import daxing.common.VCF;
import utils.IOUtils;
import utils.PStringUtils;
import utils.Tuple;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
            try {
                StringBuilder sb=new StringBuilder();
                for (Map.Entry<Integer, Long> ele: m.entrySet()){
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

//    public static void main(String[] args) {
//        int[] a=IntStream.range(0,45).toArray();
//        Test.calcalateGenotypedTaxonNum(args[0],args[1],a);
//    }
}
