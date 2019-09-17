package daxing.applets;

import daxing.common.ArrayTool;
import daxing.common.VCF;
import daxing.filterSNP.Cells;
import daxing.filterSNP.DepthInfo;
import daxing.filterSNP.Dot;
import format.position.ChrPos;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import utils.Benchmark;
import utils.IOUtils;
import utils.PStringUtils;

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
import java.util.stream.IntStream;

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

    public static void getTopRows(File inputFile, int n, File outputFile){
        try{
            List<String> lineList= Files.newBufferedReader(Paths.get(inputFile.getAbsolutePath())).lines().limit(n).collect(Collectors.toList());
            BufferedWriter bw=Files.newBufferedWriter(Paths.get(outputFile.getAbsolutePath()));
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

    public static void getTopRowsFromDir(String inputDir, int n, String outputDir){
        File[] files=IOUtils.listRecursiveFiles(new File(inputDir));
        String[] names= Arrays.stream(files).map(File::getName).toArray(String[]::new);
        IntStream.range(0, files.length).parallel().forEach(e->{
            ScriptMethods.getTopRows(files[e], n, new File(outputDir, names[e]));
        });
    }

    // chr pos depth sd输入文件，输出目录
    public static void getCellDensity(String inputFile, String outputDir){
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
    public static void getHighCumulativePos(String inputDirOfPosition, String outputFile, int num){
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

    public static void getSubsetForGraph(String inputFile, String outputFile, int num){
        try(BufferedReader br=IOUtils.getTextReader(inputFile);
            BufferedWriter bw=IOUtils.getTextWriter(outputFile)){
            String line;
            List<String> lineList=new ArrayList<>();
            String header=br.readLine();
            bw.write(header);
            bw.newLine();
            while ((line=br.readLine())!=null){
                lineList.add(line);
            }
            int[] randoms= ArrayTool.getRandomNonrepetitionArray(num, 0, lineList.size());
            Arrays.sort(randoms);
            int index;
            for (int i = 0; i < lineList.size(); i++) {
                index=Arrays.binarySearch(randoms, i);
                if (index<0) continue;
                bw.write(lineList.get(i));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    public static void addDepthGroup(String databaseFile, String queryFile, int queryNum, String outFile){
        try(BufferedReader brQuery=IOUtils.getTextReader(queryFile);
            BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            ChrPos[] database;
            database=Files.newBufferedReader(Paths.get(databaseFile)).lines().skip(1).map(PStringUtils::fastSplit).map(l->new ChrPos(Short.parseShort(l.get(0)), Integer.parseInt(l.get(1)))).sorted().toArray(ChrPos[]::new);
            String line;
            short chr;
            int pos;
            List<String> lineList;
            List<String> lines=new ArrayList<>();
            String header=brQuery.readLine();
            bw.write(header+"\t"+"Group");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            int[] chrs= Arrays.stream(database).map(ChrPos::getChromosome)
                    .mapToInt(Integer::valueOf).toArray();
            int[] poss=Arrays.stream(database).map(ChrPos::getPosition)
                    .mapToInt(Integer::intValue).toArray();
            int indexChr;
            int indexPOS;
            while ((line=brQuery.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                chr=Short.parseShort(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                indexChr= Arrays.binarySearch(chrs, chr);
                indexPOS=Arrays.binarySearch(poss, pos);
                if (indexChr>0 && indexPOS>0){
                    sb.append(line).append("\t").append(queryNum);
                    lines.add(sb.toString());
                }else {
                    sb.append(line).append("\t").append(0);
                    lines.add(sb.toString());
                }
                sb=new StringBuilder();
            }
            brQuery.close();
            for (int i = 0; i < lines.size(); i++) {
                bw.write(lines.get(i));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void getHapPos(File posAlleleDir, File hapPosDir){
        File[] files=IOUtils.listRecursiveFiles(posAlleleDir);
        BufferedReader[] brs=new BufferedReader[files.length];
        BufferedWriter[] bws=new BufferedWriter[files.length];
        for (int i = 0; i < files.length; i++) {
            brs[i]=IOUtils.getTextGzipReader(files[i].getAbsolutePath());
            bws[i]=IOUtils.getTextGzipWriter(files[i].getAbsolutePath());
        }
        List<ChrPos> chrPosList;
        StringBuilder sb=new StringBuilder();
        try{
            for (int i = 0; i < brs.length; i++) {
                chrPosList=brs[i].lines().skip(1).map(PStringUtils::fastSplit)
                        .map(e->new ChrPos(Short.parseShort(e.get(0)),Integer.parseInt(e.get(1)))).collect(Collectors.toList());
                sb=new StringBuilder();
                for (int j = 0; j < chrPosList.size(); j++) {
                    sb.append(chrPosList.get(j).getChromosome()).append("\t").append(chrPosList.get(j).getPosition());
                    bws[i].write(sb.toString());
                    bws[i].newLine();
                }
                bws[i].flush();
                bws[i].close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void caculateChr(String inputFile, String outFile){
        try(BufferedReader br=IOUtils.getTextReader(inputFile);
            BufferedWriter bw=IOUtils.getTextWriter(outFile)){
            String header, line;
            List<String> lineList;
            header=br.readLine();
            bw.write(header);
            bw.newLine();
            int[] chrID;
            int[] snpNum;
            int[] biallelicNum;
            int[] indelNum;
            int[] insertionNum;
            int[] delectionNum;
            while ((line=br.readLine())!=null){
                chrID=new int[2];
                snpNum=new int[2];
                biallelicNum=new int[2];
                indelNum=new int[2];
                insertionNum=new int[2];
                delectionNum=new int[2];
                StringBuilder sb=new StringBuilder();
                lineList=PStringUtils.fastSplit(line);
                chrID[0]=Integer.parseInt(lineList.get(0));
                snpNum[0]=Integer.parseInt(lineList.get(1));
                biallelicNum[0]=Integer.parseInt(lineList.get(2));
                indelNum[0]=Integer.parseInt(lineList.get(3));
                insertionNum[0]=Integer.parseInt(lineList.get(4));
                delectionNum[0]=Integer.parseInt(lineList.get(5));
                line=br.readLine();
                lineList=PStringUtils.fastSplit(line);
                chrID[1]=Integer.parseInt(lineList.get(0));
                snpNum[1]=Integer.parseInt(lineList.get(1));
                biallelicNum[1]=Integer.parseInt(lineList.get(2));
                indelNum[1]=Integer.parseInt(lineList.get(3));
                insertionNum[1]=Integer.parseInt(lineList.get(4));
                delectionNum[1]=Integer.parseInt(lineList.get(5));
                sb.append(VCF.chrToChrMap.get(chrID[0])).append("\t").append(snpNum[0]+snpNum[1]).append("\t").append(biallelicNum[0]+biallelicNum[1]).append("\t").append(indelNum[0]+indelNum[1]).append("\t").append(insertionNum[0]+insertionNum[1]).append("\t").append(delectionNum[0]+delectionNum[1]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void checkOutgroupAllele(String inputFile, String outDir){
        long start=System.nanoTime();
        String line;
        List<String> lineList;
        BufferedReader br;
        String header;
        TIntHashSet chrSet=new TIntHashSet();
        TIntArrayList chrList=new TIntArrayList();
        TIntArrayList posList=new TIntArrayList();
        List<String> refList=new ArrayList<>();
        List<String> outgroupList=new ArrayList<>();
        try{
            br=IOUtils.getTextReader(inputFile);
            header=br.readLine();
            while ((line=br.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                chrSet.add(Integer.parseInt(lineList.get(0)));
                chrList.add(Integer.parseInt(lineList.get(0)));
                posList.add(Integer.parseInt(lineList.get(1)));
                refList.add(lineList.get(2));
                outgroupList.add(lineList.get(3));
            }
            br.close();
            System.out.println(inputFile+" reading is completed");
            TIntArrayList chrs=new TIntArrayList(chrSet);
            int[] chrNum=chrs.toArray();
            Arrays.sort(chrNum);
            BufferedWriter[] bws=new BufferedWriter[chrNum.length];
            for (int i = 0; i < chrNum.length; i++) {
                bws[i]=IOUtils.getTextWriter(new File(outDir, chrNum[i]+".txt").getAbsolutePath());
            }
            for (int i = 0; i < chrNum.length; i++) {
                bws[i].write(header);
                bws[i].newLine();
            }
            int index=Integer.MIN_VALUE;
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < chrList.size(); i++) {
                index=Arrays.binarySearch(chrNum, chrList.get(i));
                sb.append(chrList.get(i)).append("\t").append(posList.get(i)).append("\t")
                        .append(refList.get(i)).append("\t").append(outgroupList.get(i));
                bws[index].write(sb.toString());
                bws[index].newLine();
                sb=new StringBuilder();
            }
            for (int i = 0; i < bws.length; i++) {
                bws[i].flush();
                bws[i].close();
            }
            System.out.println(inputFile+" completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

}
