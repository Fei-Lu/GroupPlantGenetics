package daxing.applets;

import daxing.common.StringTool;
import daxing.common.WheatLineage;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class Shell {

    /**
     * lastz /data1/home/daxing/reference/triticum_aestivum/bychr/chr005.fa /data1/home/daxing/reference/triticum_urartu/bychr/chr1A.fa --notransition --step=20  --nogapped --ambiguous=iupac --format=maf > /data1/home/daxing/ancestralAllele/MAF/ta_tu/wheatTauschii/Ta_chr005_vs_Au_chr1A.maf 2>/data1/home/daxing/ancestralAllele/MAF/ta_tu/wheatTauschiiLog/Ta_chr005_vs_Au_chr1ALog.txt &
     *
     * @param wheatInputDir
     * @param outgroupInputDir
     * @param outMAFDir
     * @param outSH
     */
    public static void getLastz(String wheatInputDir, String outgroupInputDir, String outMAFDir, String logFileDir, String outSH){
        File[] files1= IOUtils.listRecursiveFiles(new File(wheatInputDir));
        int[] d= WheatLineage.valueOf("D").getChrID();
        List<Integer> l= Arrays.stream(d).boxed().collect(Collectors.toList());
        Predicate<File> dPredicate= e-> l.contains(StringTool.getNumFromString(e.getName()));
        File[] files2=IOUtils.listRecursiveFiles(new File(outgroupInputDir));
        Predicate<File> p= e->e.getName().endsWith("fa");
        File[] f1= Arrays.stream(files1).filter(p).filter(dPredicate).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p).toArray(File[]::new);
        try(BufferedWriter bw=IOUtils.getTextWriter(outSH)){
            StringBuilder sb;
            sb=new StringBuilder();
            for (int i = 0; i < f1.length; i++) {
                for (int j = 0; j < f2.length; j++) {
                    sb.append("lastz ").append(f1[i].getAbsolutePath()).append(" ").append(f2[j].getAbsolutePath())
                            .append(" --notransition --step=20  --nogapped --ambiguous=iupac --format=maf > ")
                            .append(outMAFDir).append("/Ta_").append(f1[i].getName(), 0, 6).append("_vs_Au_")
                            .append(f2[j].getName(), 0, 5).append(".maf 2>").append(logFileDir+"/Ta_")
                            .append(f1[i].getName(), 0, 6).append("_vs_Au_").append(f2[j].getName(), 0, 5)
                            .append("Log.txt &");
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

    /**
     * plink --vcf chrTauschiiRecoded.vcf --allow-extra-chr --chr-set 21 --make-bed --r2  triangle  --threads 1 --out f
     * @param vcfInputDir
     * @param outDir
     * @param outShellFile
     */
    public static void getPlink(String vcfInputDir, String outDir, String outShellFile){
        File[] files=IOUtils.listRecursiveFiles(new File(vcfInputDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        String[] outNames= Arrays.stream(f1).map(File::getName).map(str->str.replaceAll("/.vcf$", ""))
                .toArray(String[]::new);
        try (BufferedWriter bw = IOUtils.getTextWriter(outShellFile)) {
            StringBuilder sb;
            for (int i = 0; i < f1.length; i++) {
                sb=new StringBuilder();
                sb.append("plink --vcf --allow-extra-chr ").append(f1[i].getAbsolutePath()).append(" --chr-set 21 ")
                        .append("--make-bed ").append("--r2 --matrix --out ")
                        .append(new File(outDir, outNames[i]).getAbsolutePath());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * parallel for Shell script
     * @param inputSh
     * @param threadNum
     * @param outputSh
     */
    public static void getShellParallelScript(String inputSh, int threadNum, String outputSh){
        try (BufferedReader br = IOUtils.getTextReader(inputSh);
             BufferedWriter bw=IOUtils.getTextWriter(outputSh)) {
            int count=0;
            String temp;
            StringBuilder sb;
            bw.write("echo $(date '+%Y-%m-%d %H:%M:%S') start");
            bw.newLine();
            while ((temp=br.readLine())!=null){
                count++;
                if (count % threadNum==0){
                    sb=new StringBuilder();
                    sb.append(temp).append("\n");
                    sb.append("wait").append("\n");
                    bw.write(sb.toString());
                    continue;
                }
                bw.write(temp);
                bw.newLine();
            }
            bw.write("echo $(date '+%Y-%m-%d %H:%M:%S') end");
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
