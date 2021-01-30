package daxing.applets;

import daxing.common.ChrConvertionRule;
import daxing.common.DateTime;
import daxing.common.NumberTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang.text.StrBuilder;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;

public class Est_sfs {

    public static void getEstsfsInputFile(String chrOutDir,
                                          String mergedDir){
//        List<File> dbFiles=IOUtils.getVisibleFileListInDir(annotationDBDir);
//        String[] outNames=dbFiles.stream().map(File::getName).map(str->str.substring(0,6)+"_data.file.txt").toArray(String[]::new);
//        IntStream.range(0, dbFiles.size()).forEach(e->getEstsfsInputFile(dbFiles.get(e), new File(chrOutDir,
//                outNames[e]), chrConvertionRule));
        fastMergeVCFtoLineage(chrOutDir, mergedDir);
    }

    private static void getEstsfsInputFile(File dbFile, File outFile, String chrConvertionRule){
        String[] bases={"A","C", "G","T"};
        int[] counts;
        ChrConvertionRule c=new ChrConvertionRule(Paths.get(chrConvertionRule));
        try (BufferedReader br = IOUtils.getTextGzipReader(dbFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line, chr, minor, major;
            double maf=-1;
            List<String> temp;
//            bw.write("chr\tpos\tACGT\n");
            br.readLine();
            StrBuilder sb;
            String[] countsStr;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                int chrID=Integer.parseInt(temp.get(0));
                chr=c.getRefChrFromVCFChr(chrID);
                int pos=c.getRefPosFromVCFChrPos(chrID, Integer.parseInt(temp.get(1)));
                major=temp.get(4);
                minor=temp.get(5);
                int indexOfMinor=Arrays.binarySearch(bases, minor);
                int indexOfMajor= Arrays.binarySearch(bases, major);
                maf=Double.parseDouble(temp.get(6));
                int countOfMinor=(int) (NumberTool.format(maf, 2)*100);
                int countOfMajor=100-countOfMinor;
                counts=new int[4];
                counts[indexOfMinor]=countOfMinor;
                counts[indexOfMajor]=countOfMajor;
                sb=new StrBuilder(32);
                countsStr= Arrays.stream(counts).boxed().map(s->String.valueOf(s)).toArray(String[]::new);
                sb.append(chr).append("\t").append(pos).append("\t").append(String.join(",", countsStr)).toString();
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void fastMergeVCFtoLineage(String inputDir, String outDir){
        System.out.println(DateTime.getDateTimeOfNow()+ " start");
        File[] files=new File(inputDir).listFiles();
        Predicate<File> hidden=File::isHidden;
        Predicate<File> p= hidden.negate().and(f->f.getName().toLowerCase().startsWith("chr"));
        File[] f=Arrays.stream(files).filter(p).sorted().toArray(File[]::new);
        TIntArrayList[] abd=new TIntArrayList[3];
        abd[0]=new TIntArrayList(WheatLineage.valueOf("A").getChrID());
        abd[1]=new TIntArrayList(WheatLineage.valueOf("B").getChrID());
        abd[2]=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        Predicate<File> ap=fa->abd[0].contains(Integer.parseInt(fa.getName().substring(3,6)));
        Predicate<File> bp=fa->abd[1].contains(Integer.parseInt(fa.getName().substring(3,6)));
        Predicate<File> dp=fa->abd[2].contains(Integer.parseInt(fa.getName().substring(3,6)));
        File[][] abd_lineageFile=new File[3][];
        abd_lineageFile[0]= Arrays.stream(f).filter(ap).toArray(File[]::new);
        abd_lineageFile[1]= Arrays.stream(f).filter(bp).toArray(File[]::new);
        abd_lineageFile[2]= Arrays.stream(f).filter(dp).toArray(File[]::new);
        String[] outNames={"A.data.txt", "B.data.txt", "D.data.txt"};
        BufferedWriter[] bws=new BufferedWriter[3];
        for (int i = 0; i < bws.length; i++) {
            bws[i]=IOUtils.getTextWriter(new File(outDir, outNames[i]).getAbsolutePath());
        }
        try {
            BufferedReader br;
            StringBuilder sb;
            String line;
            for (int i = 0; i < abd_lineageFile.length; i++) {
                for (int j = 0; j < abd_lineageFile[i].length; j++) {
                    sb=new StringBuilder(1000);
                    br=IOUtils.getTextReader(abd_lineageFile[i][j].getAbsolutePath());
                    while ((line=br.readLine())!=null){
                        bws[i].write(line);
                        bws[i].newLine();
                    }
                    br.close();
                }
                bws[i].flush();
                bws[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(DateTime.getDateTimeOfNow()+" end");
    }

}
