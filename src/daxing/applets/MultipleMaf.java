package daxing.applets;

import daxing.common.ChrConvertionRule;
import daxing.common.ChrRange;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import org.apache.commons.lang.StringUtils;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.function.Predicate;

public class MultipleMaf {

    public static List<ChrRange> getRange(String multileMafFile, String abd){
        List<ChrRange> rangeInterfaceList=new ArrayList<>();
        try (BufferedReader br1 = IOUtils.getTextReader(multileMafFile)) {
            String line;
            String[] temp;
            List<String> lines;
            ChrRange range;
            int refStart, len;
            String chr;
            StringBuilder sb;
            while ((line=br1.readLine()).startsWith("#")){}
            do {
                lines=new ArrayList<>();
                while ((line=br1.readLine()).startsWith("s")){
                    lines.add(line);
                }
                temp= StringUtils.split(lines.get(0), " .");
                sb=new StringBuilder();
                chr=sb.append(temp[2].substring(3)).append(abd).toString();
                refStart=Integer.parseInt(temp[3])+1;
                len=Integer.parseInt(temp[4]);
                range=new ChrRange(chr, refStart, refStart+len+1);
                rangeInterfaceList.add(range);
            }while ((line=br1.readLine())!=null && !(line.equals("##eof maf")));

        }catch (Exception e){
            e.printStackTrace();
        }
        return rangeInterfaceList;
    }

    /**
     * 输入文件是 chrA.subgenome.jgerpFile chrB.subgenome.jgerpFile chrD.subgenome.jgerpFile
     * @param jgerpFile
     * @param recordRange
     * @param outFile 没有排序
     */
    public static void getGerpFile(String jgerpFile, List<ChrRange> recordRange, String outFile, ChrConvertionRule c){
        try (BufferedReader br = IOUtils.getTextReader(jgerpFile);
             BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            bw.write("Chr\tPos\tGerpNeutralRate\tGerpScore");
            bw.newLine();
            String line;
            String[] temp;
            List<String> lines=new ArrayList<>();
            StringBuilder sb;
            ChrRange range;
            String chr;
            int pos;
            int vcfChr, vcfPos;
            for (int i = 0; i < recordRange.size(); i++) {
                range=recordRange.get(i);
                for (int j = 0; j < range.getLength(); j++) {
                    chr=range.getChr();
                    pos=range.getStart()+j;
                    while ((line=br.readLine())!=null){
                        temp=StringUtils.split(line);
                        sb=new StringBuilder();
                        vcfChr=c.getVCFChrFromRefChrPos(chr, pos);
                        vcfPos=c.getVCFPosFromRefChrPos(chr, pos);
                        sb.append(vcfChr).append("\t").append(vcfPos).append("\t").append(temp[0]).append("\t").append(temp[1]);
                        bw.write(sb.toString());
                        bw.newLine();
                        break;
                    }
                }
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     *
     * @param gerpFileDir
     * @param outDir 没有排序
     */
    public static void splitGerpFile(String gerpFileDir, String outDir){
        File[] files=new File(gerpFileDir).listFiles();
        Predicate<File> p=File::isHidden;
        File[] fs= Arrays.stream(files).filter(p.negate()).filter(str->str.getName().endsWith("gerp++"))
                .sorted().toArray(File[]::new);
        for (int e = 0; e < fs.length; e++) {
            int[] chrs=MultipleMaf.getChrs(fs[e].getAbsolutePath());
            Integer key;
            BufferedWriter value;
            Map<Integer, BufferedWriter> chrBuffwriterMap=new HashMap<>();
            for (int i = 0; i < chrs.length; i++) {
                key=chrs[i];
                value=IOUtils.getTextWriter(outDir+"/chr"+ PStringUtils.getNDigitNumber(3, key)+".gerpScore");
                chrBuffwriterMap.put(key, value);
            }
            try {
                String line;
                String[] temp;
                try (BufferedReader br = IOUtils.getTextReader(fs[e].getAbsolutePath())) {
                    String header=br.readLine();
                    for(Map.Entry<Integer, BufferedWriter> entry: chrBuffwriterMap.entrySet()){
                        entry.getValue().write(header);
                        entry.getValue().newLine();
                    }
                    while ((line=br.readLine())!=null){
                        temp=StringUtils.split(line);
                        key=Integer.parseInt(temp[0]);
                        value=chrBuffwriterMap.get(key);
                        value.write(line);
                        value.newLine();
                    }
                }catch (Exception exce){
                    exce.printStackTrace();
                }
                for (Map.Entry<Integer, BufferedWriter> entry: chrBuffwriterMap.entrySet()){
                    entry.getValue().flush();
                    entry.getValue().close();
                }
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    private static int[] getChrs(String gerpFile){
        int[] chrs=null;
        try (BufferedReader br = IOUtils.getTextReader(gerpFile)) {
            String line;
            br.readLine();
            String[] temp;
            TIntHashSet chrSet=new TIntHashSet();
            while ((line=br.readLine())!=null){
                temp=StringUtils.split(line);
                chrSet.add(Integer.parseInt(temp[0]));
            }
            TIntArrayList chrList=new TIntArrayList(chrSet);
            chrList.sort();
            chrs=chrList.toArray();
        }catch (Exception e){
            e.printStackTrace();
        }
        return chrs;
    }

//    public static void main(String[] args) {
//        ChrConvertionRule c=new ChrConvertionRule(Paths.get(""));
//        List<RangeInterface> rangeInterfaceList=MultipleMaf.getRange("", "", c);
//        MultipleMaf.getGerpFile("", rangeInterfaceList, "");
//        MultipleMaf.splitGerpFile("", "");
//    }
}
