package daxing.load.ancestralSite.complementary;

import daxing.common.IOTool;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class WindowLoadABD {


    public enum LoadType{
        Syn, Non, del
    }

    public static void mergeIndividual(String inputDir, String outDir){
        String[] taxa={"C1","C2","C11","C12"};
        List<File> fileList= IOTool.getVisibleFileRecursiveDir(inputDir);
        Predicate<File> p=file -> file.getName().contains("triadposA");
        BufferedReader brSyn, brNon, brDel;
        BufferedWriter bw;
        StringBuilder sb=new StringBuilder();
        try {
            for (int i = 0; i < taxa.length; i++) {
                int index=i;
                Predicate<File> pTaxa=file -> PStringUtils.fastSplit(file.getName(),".").get(0).equals(taxa[index]);
                List<File> files=fileList.stream().filter(p.and(pTaxa)).collect(Collectors.toList());
                brSyn=IOTool.getReader(files.get(2));
                brNon=IOTool.getReader(fileList.get(1));
                brDel=IOTool.getReader(fileList.get(0));
                bw=IOTool.getTextWriter(new File(outDir, taxa[i]+".triadA.10M_window_1M_step.txt.gz"));
                bw.write("Chr\tWindowStart\tWindowEnd\tTriadsNum\tMean_NormalizedTriadLoad_Syn" +
                        "\tMean_NormalizedTriadLoad_Non\tMean_NormalizedTriadLoad_Del");
                bw.newLine();
                String line;
                List<String> tempSyn, tempNon, tempDel;
                brSyn.readLine();
                brNon.readLine();
                brDel.readLine();
                while ((line=brSyn.readLine())!=null){
                    tempSyn= PStringUtils.fastSplit(line);
                    tempNon=PStringUtils.fastSplit(brNon.readLine());
                    tempDel=PStringUtils.fastSplit(brDel.readLine());
                    sb.setLength(0);
                    sb.append(tempSyn.get(0)).append("\t").append(tempSyn.get(1)).append("\t");
                    sb.append(tempSyn.get(2)).append("\t").append(tempSyn.get(3));
                    sb.append(tempSyn.get(4)).append("\t");
                    if (!(tempSyn.get(0).equals(tempNon.get(0)) && tempSyn.get(1).equals(tempNon.get(1)))){
                        System.out.println("check, exist error");
                        System.out.println(String.join("\t",tempSyn));
                        System.out.println(String.join("\t",tempNon));
                        System.out.println(String.join("\t",tempDel));
                        System.exit(1);
                    }
                    sb.append(tempNon.get(4)).append("\t");
                    if (!(tempSyn.get(0).equals(tempDel.get(0)) && tempSyn.get(1).equals(tempDel.get(1)))){
                        System.out.println("check, exist error");
                        System.exit(1);
                    }
                    sb.append(tempDel.get(4)).append("\t");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                brSyn.close();
                brNon.close();
                brDel.close();
                bw.flush();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }



}
