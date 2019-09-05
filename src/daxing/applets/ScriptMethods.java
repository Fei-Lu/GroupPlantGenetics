package daxing.applets;

import daxing.filterSNP.Cells;
import daxing.filterSNP.DepthInfo;
import daxing.filterSNP.Dot;
import utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class ScriptMethods {

    public static void getLastz(String wheatInputDir, String outgroupInputDir, String outMAFDir, String outSH){
        File[] files1= IOUtils.listRecursiveFiles(new File(wheatInputDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(outgroupInputDir));
        Predicate<File> p= e->e.getName().endsWith("fa");
        File[] f1= Arrays.stream(files1).filter(p).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p).toArray(File[]::new);
        try(BufferedWriter bw=IOUtils.getNIOTextWriter(outSH)){
            StringBuilder sb;
            sb=new StringBuilder();
            for (int i = 0; i < f1.length; i++) {
                for (int j = 0; j < f2.length; j++) {
                    sb.append("lastz ").append(f1[i].getAbsolutePath()).append(" ").append(f2[j].getAbsolutePath())
                            .append(" --notransition --ambiguous=iupac --step=20 --nogapped --format=maf > ").append(outMAFDir)
                            .append("/Ta_").append(f1[i].getName().substring(0, 6)).append("_vs_At_")
                            .append(f2[j].getName().substring(0,5)).append(".maf"+" &");
                    bw.write(sb.toString());
                    bw.newLine();
                    sb=new StringBuilder();
                }
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void getTopRows(String inputFile, int n, String outputFile){
        try{
            List<String> lineList= Files.newBufferedReader(Paths.get(inputFile)).lines().limit(n).collect(Collectors.toList());
            BufferedWriter bw=Files.newBufferedWriter(Paths.get(outputFile));
            for (int i = 0; i < lineList.size(); i++) {
                bw.write(lineList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    // chr pos depth sd输入文件，输出目录
    public static void filterSNP1(String inputFile, String outputDir){
        DepthInfo depthInfo=new DepthInfo(inputFile);
        List<Dot> dotList=depthInfo.getDotList();
        Cells cells=new Cells(30F, 20F, 100);
        int[] indexOfDepthSD;
        for (int i = 0; i < dotList.size(); i++) {
            indexOfDepthSD= cells.binarySearch(dotList.get(i));
            cells.getCell(indexOfDepthSD).add(dotList.get(i));
        }
        cells.write(outputDir);
    }

    // 根据输出的postition目录，筛选累计密度达到指定百分比的格子所对应的ChrPos信息
    public static void filterSNP2(String inputDirOfPosition, String outputFile, int num){
        File[] files=IOUtils.listRecursiveFiles(new File(inputDirOfPosition));
        List<String> chrPosList=new ArrayList<>();
        BufferedReader[] brs=new BufferedReader[num+1];
        try{
            for (int i = 0; i < brs.length; i++) {
                brs[i]=IOUtils.getTextReader(files[i].getAbsolutePath());
                String line;
                brs[i].readLine();
                List<String> lines=new ArrayList<>();
                while ((line=brs[i].readLine())!=null){
                    lines.add(line);
                }
                System.out.println(i+"\t"+lines.size());
                chrPosList.addAll(lines);
                brs[i].close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }

        BufferedWriter bw=IOUtils.getTextWriter(outputFile);
//        int[] randoms= ArrayTool.getRandomNonrepetitionArray(5000, 0, chrPosList.size());
//        Arrays.sort(randoms);
        try {
            bw.write("CHR"+"\t"+"POS"+"\t"+"AverageDepth"+"\t"+"SD"+"\n");
//            int index;
            for (int i = 0; i < chrPosList.size(); i++) {
//                index=Arrays.binarySearch(randoms, i);
//                if (index<0) continue;
                bw.write(chrPosList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
