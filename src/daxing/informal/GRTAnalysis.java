package daxing.informal;

import daxing.common.StringTool;
import utils.Benchmark;
import utils.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;

public class GRTAnalysis {

    GRTAnalysis(){
//        GRTAnalysis.getCutSiteIndex("/Users/xudaxing/Desktop/chr001.fa",
//                "/Users/xudaxing/Desktop/cutter1.txt","/Users/xudaxing/Desktop/cutter2.txt");
//        GRTAnalysis.mergeCutSite("/Users/xudaxing/Desktop/cutter1.txt",
//                "/Users/xudaxing/Desktop/cutter2.txt","/Users/xudaxing/Desktop/cutSite.txt");
//        GRTAnalysis.calculateGRTCutSite("/Users/xudaxing/Desktop/cutter1.txt",
//                "/Users/xudaxing/Desktop/cutter2.txt",
//                "/Users/xudaxing/Desktop/cutSite.txt",
//                "/Users/xudaxing/Desktop/range200_500.txt",
//                "/Users/xudaxing/Desktop/range500_1000.txt");
    }

    /**
     * 返回一个染色体上含有的GGATCC和CCGG切点index
     * @param chrFastaFile 染色体fasta文件
     * @param cutter1OutputFile BamH1(GGATCC)切点输出文件
     * @param cutter2OutputFile Msp1(CCGG)切点输出文件
     */
    public static void getCutSiteIndex(String chrFastaFile, String cutter1OutputFile, String cutter2OutputFile){
        String str=null;
        StringBuilder sb=new StringBuilder();
        String temp=null;
        try(BufferedReader br= IOUtils.getTextReader(chrFastaFile);
            BufferedWriter bw1=IOUtils.getTextWriter(cutter1OutputFile);
            BufferedWriter bw2=IOUtils.getTextWriter(cutter2OutputFile)){
            br.readLine();
            while ((temp=br.readLine())!=null){
                sb.append(temp);
            }
            str=sb.toString();
            int[] frequency1= StringTool.getFrequencyOfSubStr(str,"GGATCC");
            int[] frequency2=StringTool.getFrequencyOfSubStr(str,"CCGG");
            sb=new StringBuilder();
            for(int i=0;i<frequency1.length;i++){
                sb.append(frequency1[i]);
                sb.append("\n");
            }
            bw1.write(sb.toString());
            sb=new StringBuilder();
            for(int i=0;i<frequency2.length;i++){
                sb.append(frequency2[i]);
                sb.append("\n");
            }
            bw2.write(sb.toString());
            bw1.flush();
            bw2.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * 将BamH1和Msp1的切点进行融合后，只保留200-1000bp的片段（即去除非200-1000bp的切点）
     * @param cutter1InputFile BamH1(GGATCC)切点输入文件
     * @param cutter2InputFile Msp1(CCGG)切点输入文件
     * @param cutSiteOutputFile merge后的切点输出文件
     */
    public static void mergeCutSite(String cutter1InputFile, String cutter2InputFile, String cutSiteOutputFile){
        try(BufferedReader br1=IOUtils.getTextReader(cutter1InputFile);
            BufferedReader br2=IOUtils.getTextReader(cutter2InputFile);
            BufferedWriter bw=IOUtils.getTextWriter(cutSiteOutputFile)){
            String temp;
            List<Integer> GATCC=new ArrayList<>();
            while ((temp=br1.readLine())!=null){
                GATCC.add(Integer.valueOf(temp));
            }
            while ((temp=br2.readLine())!=null){
                GATCC.add(Integer.valueOf(temp));
            }
            Collections.sort(GATCC);
            System.out.println(GATCC.size());
            int previousRange, nextRange;
            for(int i=1;i<GATCC.size()-1;i++){
                previousRange=GATCC.get(i)-GATCC.get(i-1);
                nextRange=GATCC.get(i+1)-GATCC.get(i);
                if((previousRange<200||previousRange>1000) && (nextRange<200||nextRange>1000)){
                    GATCC.remove(i);
                    i--;
                }
            }
            System.out.println(GATCC.size());
            for(Integer e:GATCC){
                bw.write(String.valueOf(e));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * GRT切点组合（即5'端是BamH1切点，3'端是Msp1切点）产生的200-500bp和500-1000bp的片段，注意一个BamH1切点可以对应连续不同的Msp1切点
     * @param cutter1InputFile BamH1(GGATCC)切点输入文件
     * @param cutter2InputFile Msp1(CCGG)切点输入文件
     * @param cutSiteInputFile merge后的切点输入文件
     * @param range200_500OutFile 产生200-500bp的GRT切点组合输出文件
     * @param range500_1000OutFile 产生500-1000bp的GRT切点组合输出文件
     */
    public static void calculateGRTCutSite(String cutter1InputFile, String cutter2InputFile, String cutSiteInputFile, String range200_500OutFile, String range500_1000OutFile){
        try(BufferedReader br1=IOUtils.getTextReader(cutter1InputFile);
            BufferedReader br2=IOUtils.getTextReader(cutter2InputFile);
            BufferedReader br=IOUtils.getTextReader(cutSiteInputFile);
            BufferedWriter bw1=IOUtils.getTextWriter(range200_500OutFile);
            BufferedWriter bw2=IOUtils.getTextWriter(range500_1000OutFile)){
            Set<Integer> l1=new HashSet<>();
            Set<Integer> l2=new HashSet<>();
            Set<Integer> bamH1CutSite=new HashSet<>();
            Set<Integer> msp1CutSite=new HashSet<>();
            String temp;
            while ((temp=br1.readLine())!=null){
                l1.add(Integer.valueOf(temp));
            }
            while ((temp=br2.readLine())!=null){
                l2.add(Integer.valueOf(temp));
            }
            while ((temp=br.readLine())!=null){
                bamH1CutSite.add(Integer.valueOf(temp));
            }
            msp1CutSite.addAll(bamH1CutSite);
            bamH1CutSite.retainAll(l1);
            msp1CutSite.retainAll(l2);
            List<Integer> list1=new ArrayList<>(bamH1CutSite);
            List<Integer> list2=new ArrayList<>(msp1CutSite);
            Collections.sort(list1);
            Collections.sort(list2);
            List<int[]> range200_500=new ArrayList<>();
            List<int[]> range500_1000=new ArrayList<>();
            int[] a;
            long start=System.nanoTime();
            int cnt1=0, cnt2=0;
            for(int i=0;i<list1.size();i++){
                for(int j=0;j<list2.size();j++){
                    int range=list2.get(j)-list1.get(i);
                    if(range<200) continue;
                    else if (range>1000) break;
                    else if(range<500){
                        a=new int[2];
                        a[0]=list1.get(i);
                        a[1]=list2.get(j);
                        range200_500.add(a);
                        cnt1++;
                        if(cnt1%1000 == 0){
                            System.out.println(cnt1+" cutsites covering 200-500bp have been found");
                        }
                    }
                    else {
                        a=new int[2];
                        a[0]=list1.get(i);
                        a[1]=list2.get(j);
                        range500_1000.add(a);
                        cnt2++;
                        if(cnt2%1000 == 0){
                            System.out.println(cnt2+" cutsites covering 500-1000bp have been found");
                        }
                    }
                }
            }
            for(int i=0;i<range200_500.size();i++){
                bw1.write(String.valueOf(range200_500.get(i)[0]));
                bw1.write("\t");
                bw1.write(String.valueOf(range200_500.get(i)[1]));
                bw1.newLine();
            }
            for(int i=0;i<range500_1000.size();i++){
                bw2.write(String.valueOf(range500_1000.get(i)[0]));
                bw2.write("\t");
                bw2.write(String.valueOf(range500_1000.get(i)[1]));
                bw2.newLine();
            }
            bw1.flush();bw2.flush();
            bw1.close();bw2.close();
            System.out.println("completed in " + String.format("%.4f", Benchmark.getTimeSpanMinutes(start)) + " minutes");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * GRT切点组合（即5'端是BamH1切点，3'端是Msp1切点）产生的200-500bp和500-1000bp的片段，注意一个BamH1切点仅对应一个不同的Msp1切点
     * @param cutter1InputFile BamH1(GGATCC)切点输入文件
     * @param cutter2InputFile Msp1(CCGG)切点输入文件
     * @param cutSiteInputFile merge后的切点输入文件
     * @param range200_500OutFile 产生200-500bp的GRT切点组合输出文件
     * @param range500_1000OutFile 产生500-1000bp的GRT切点组合输出文件
     */
    public static void calculateGRTCutSite2(String cutter1InputFile, String cutter2InputFile, String cutSiteInputFile, String range200_500OutFile, String range500_1000OutFile){
        try(BufferedReader br1=IOUtils.getTextReader(cutter1InputFile);
            BufferedReader br2=IOUtils.getTextReader(cutter2InputFile);
            BufferedReader br=IOUtils.getTextReader(cutSiteInputFile);
            BufferedWriter bw1=IOUtils.getTextWriter(range200_500OutFile);
            BufferedWriter bw2=IOUtils.getTextWriter(range500_1000OutFile)){
            Set<Integer> l1=new HashSet<>();
            Set<Integer> l2=new HashSet<>();
            Set<Integer> bamH1CutSite=new HashSet<>();
            Set<Integer> msp1CutSite=new HashSet<>();
            String temp;
            while ((temp=br1.readLine())!=null){
                l1.add(Integer.valueOf(temp));
            }
            while ((temp=br2.readLine())!=null){
                l2.add(Integer.valueOf(temp));
            }
            while ((temp=br.readLine())!=null){
                bamH1CutSite.add(Integer.valueOf(temp));
            }
            msp1CutSite.addAll(bamH1CutSite);
            bamH1CutSite.retainAll(l1);
            msp1CutSite.retainAll(l2);
            List<Integer> list1=new ArrayList<>(bamH1CutSite);
            List<Integer> list2=new ArrayList<>(msp1CutSite);
            Collections.sort(list1);
            Collections.sort(list2);
            List<int[]> range200_500=new ArrayList<>();
            List<int[]> range500_1000=new ArrayList<>();
            int[] a;
            long start=System.nanoTime();
            int cnt1=0, cnt2=0;
            for(int i=0;i<list1.size();i++){
                for(int j=0;j<list2.size();j++){
                    int range=list2.get(j)-list1.get(i);
                    if(range<200) continue;
                    else if (range>1000) break;
                    else if(range<500){
                        a=new int[2];
                        a[0]=list1.get(i);
                        a[1]=list2.get(j);
                        range200_500.add(a);
                        cnt1++;
                        if(cnt1%1000 == 0){
                            System.out.println(cnt1+" cutsites covering 200-500bp have been found");
                        }
                        break;
                    }
                    else {
                        a=new int[2];
                        a[0]=list1.get(i);
                        a[1]=list2.get(j);
                        range500_1000.add(a);
                        cnt2++;
                        if(cnt2%1000 == 0){
                            System.out.println(cnt2+" cutsites covering 500-1000bp have been found");
                        }
                        break;
                    }
                }
            }
            for(int i=0;i<range200_500.size();i++){
                bw1.write(String.valueOf(range200_500.get(i)[0]));
                bw1.write("\t");
                bw1.write(String.valueOf(range200_500.get(i)[1]));
                bw1.newLine();
            }
            for(int i=0;i<range500_1000.size();i++){
                bw2.write(String.valueOf(range500_1000.get(i)[0]));
                bw2.write("\t");
                bw2.write(String.valueOf(range500_1000.get(i)[1]));
                bw2.newLine();
            }
            bw1.flush();bw2.flush();
            bw1.close();bw2.close();
            System.out.println("completed in " + String.format("%.4f", Benchmark.getTimeSpanMinutes(start)) + " minutes");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

//    public static void main(String[] args){
//        new GRTAnalysis();
//    }
}
