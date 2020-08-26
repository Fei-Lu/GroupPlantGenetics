package daxing.applets;

import com.google.common.io.Files;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.WheatLineage;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class PopsStaticsVmapII {

    private static String[] p2={"Cultivar","Landrace_Europe"};
    private static String[] p3_daxing={"WildEmmer","DomesticatedEmmer","FreeThreshTetraploid","Ae.tauschii"};
    private static String[] p3_aoyue={"Wild_emmer","Domesticated_emmer","Free_threshing_tetraploid","Ae.tauschii"};
    private static String header="CHROM\tBIN_MIDDLE\tTYPE_PARAMETER\tVALUE\tpopulation\tpop";

    public static void popsStaticsVmapII(String fdInputDir, String fstInputDir, String piInputDir, String outDir){
        List<File> fdFiles= IOTool.getVisibleFileRecursiveDir(fdInputDir).stream().filter(getFdPredict()).collect(Collectors.toList());
        List<File> fstFiles=IOTool.getVisibleFileRecursiveDir(fstInputDir).stream().filter(getFstPredict()).collect(Collectors.toList());
        List<File> piFiles=IOTool.getVisibleFileRecursiveDir(piInputDir).stream().filter(getPiPredict()).collect(Collectors.toList());
        String[] subdir={"001_fd","002_fst","003_pi", "004_all", "005_merge"};
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> dirs=IOTool.getVisibleDir(outDir);
        transformFdToPlotFormat(fdFiles, dirs.get(0).getAbsolutePath());
        transformFstToPlotFormat(fstFiles, dirs.get(1).getAbsolutePath());
        transformPiToPlotFormat(piFiles, dirs.get(2).getAbsolutePath());
        mergeFdFstPi(dirs.get(0).getAbsolutePath(),dirs.get(1).getAbsolutePath(), dirs.get(2).getAbsolutePath(),
                dirs.get(3).getAbsolutePath());
        mergeFdFstPi(dirs.get(3), dirs.get(4));
    }

    private static void mergeFdFstPi(File inputDir, File outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir.getAbsolutePath());
        Predicate<File> d=f->f.getName().substring(4,5).equals("D");
        List<File> abFiles=files.stream().filter(d.negate()).collect(Collectors.toList());
        List<File> dFiles=files.stream().filter(d).collect(Collectors.toList());
        String[] subdir={"001_ab","002_d"};
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> subDirs=IOTool.getVisibleDir(outDir.getAbsolutePath());
        String outName;
        List<File> temp;
        for (int i = 0; i < abFiles.size(); i=i+3) {
            outName=abFiles.get(i).getName().substring(0,5)+".txt";
            temp=new ArrayList<>();
            for (int j = 0; j < 3; j++) {
                temp.add(abFiles.get(i+j));
            }
            RowTableTool.mergeRowTables(temp, new File(subDirs.get(0), outName).getAbsolutePath());
        }
        try {
            for (int i = 0; i < dFiles.size(); i++) {
                outName=dFiles.get(i).getName().substring(0,5)+".txt";
                Files.copy(dFiles.get(i), new File(subDirs.get(1), outName));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void mergeFdFstPi(String fdDir, String fstDir, String piDir, String outDir){
        List<File> fdFiles= IOUtils.getVisibleFileListInDir(fdDir);
        List<File> fstFiles= IOUtils.getVisibleFileListInDir(fstDir);
        List<File> piFiles= IOUtils.getVisibleFileListInDir(piDir);
        String[] outNames=fdFiles.stream().map(File::getName).toArray(String[]::new);
        IntStream.range(0, fdFiles.size()).forEach(e->mergeFdFstPi(fdFiles.get(e), fstFiles.get(e), piFiles.get(e),
                new File(outDir, outNames[e]).getAbsolutePath()));
    }

    private static void mergeFdFstPi(File fdFile, File fstFile, File piFile, String outFile){
        try (BufferedReader br1 = IOTool.getReader(fdFile);
             BufferedReader br2 =IOTool.getReader(fstFile);
             BufferedReader br3 =IOTool.getReader(piFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            br1.readLine();
            br2.readLine();
            br3.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            List<String> temp;
            while ((line=br1.readLine())!=null){
                bw.write(line);
                bw.newLine();
            }
            while ((line=br2.readLine())!=null){
                bw.write(line);
                bw.newLine();
            }
            while ((line=br3.readLine())!=null){
                bw.write(line);
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void transformPiToPlotFormat(List<File> piFiles, String outDir){
        Predicate<File> p_diploid=f->PStringUtils.fastSplit(f.getName(), "_chr").get(0).equals(p3_aoyue[3]);
        Predicate<File>[] p_tetraploid=new Predicate[3];
        Predicate<File>[] p_hexaploid=new Predicate[2];
        for (int i = 0; i < p3_aoyue.length-1; i++) {
            int index=i;
            p_tetraploid[i]=f->PStringUtils.fastSplit(f.getName(), "_chr").get(0).equals(p3_aoyue[index]);
        }
        for (int i = 0; i < p2.length; i++) {
            int index=i;
            p_hexaploid[i]=f->PStringUtils.fastSplit(f.getName(), "_chr").get(0).equals(p2[index]);
        }
        List<File> diploid=piFiles.stream().filter(p_diploid).collect(Collectors.toList());
        List<File> wildEmmer=piFiles.stream().filter(p_tetraploid[0]).collect(Collectors.toList());
        List<File> domesticatedEmmer=piFiles.stream().filter(p_tetraploid[1]).collect(Collectors.toList());
        List<File> freeThresh=piFiles.stream().filter(p_tetraploid[2]).collect(Collectors.toList());
        List<File> cultivarFiles=piFiles.stream().filter(p_hexaploid[0]).collect(Collectors.toList());
        List<File> landraceFiles=piFiles.stream().filter(p_hexaploid[1]).collect(Collectors.toList());
        List<File> tetraploidFiles=new ArrayList<>();
        tetraploidFiles.addAll(wildEmmer);
        tetraploidFiles.addAll(domesticatedEmmer);
        tetraploidFiles.addAll(freeThresh);
        Collections.sort(tetraploidFiles, Comparator.comparing(File::getName));
        String[] abdChrs= cultivarFiles.stream().map(File::getName).map(s->PStringUtils.fastSplit(s,"_chr").get(1))
                .map(s->s.substring(0,2)).sorted().toArray(String[]::new);
        String chr, outFileName;
        int index;
        String fileName, pop;
        for (int i = 0; i < tetraploidFiles.size(); i++) {
            fileName=tetraploidFiles.get(i).getName();
            pop=PStringUtils.fastSplit(fileName, "_chr").get(0);
            chr=PStringUtils.fastSplit(fileName, "_chr").get(1).substring(0,2);
            index=Arrays.binarySearch(abdChrs, chr);
            outFileName="chr"+chr+"_"+pop+".txt";
            mergePiFiles(cultivarFiles.get(index), p2[0], landraceFiles.get(index), p2[1], tetraploidFiles.get(i), pop,
                    new File(outDir, outFileName));
        }
        for (int i = 0; i < diploid.size(); i++) {
            fileName=diploid.get(i).getName();
            pop=PStringUtils.fastSplit(fileName, "_chr").get(0);
            chr=PStringUtils.fastSplit(fileName, "_chr").get(1).substring(0,2);
            index=Arrays.binarySearch(abdChrs, chr);
            outFileName="chr"+chr+"_"+pop+".txt";
            mergePiFiles(cultivarFiles.get(index), p2[0], landraceFiles.get(index), p2[1], diploid.get(i), pop,
                    new File(outDir, outFileName));
        }

    }

    private static void mergePiFiles(File hexaploidFile1, String p2File1, File hexaploidFile2, String p2File2,
                                     File file3, String pop3, File outFile){
        try (BufferedReader brHexaploid1 = IOTool.getReader(hexaploidFile1);
             BufferedReader brHexaploid2 = IOTool.getReader(hexaploidFile2);
             BufferedReader br = IOTool.getReader(file3);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            brHexaploid1.readLine();
            brHexaploid2.readLine();
            br.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            double pi, piHexaploid1, piHexaploid2;
            double piRatio1, piRatio2;
            int binMiddle;
            List<String> temp;
            StringBuilder sb;
            while ((line=brHexaploid1.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                piHexaploid1=Double.parseDouble(temp.get(4));
                line=brHexaploid2.readLine();
                temp=PStringUtils.fastSplit(line);
                piHexaploid2=Double.parseDouble(temp.get(4));
                line=br.readLine();
                temp=PStringUtils.fastSplit(line);
                pi=Double.parseDouble(temp.get(4));
                piRatio1=pi/piHexaploid1;
                piRatio2=pi/piHexaploid2;
                binMiddle=(Integer.parseInt(temp.get(1))+Integer.parseInt(temp.get(2)))/2;
                sb=new StringBuilder();
                sb.append(temp.get(0)).append("\t").append(binMiddle).append("\t").append("piRatio").append("\t");
                sb.append(piRatio1).append("\t").append(p2File1).append("\t").append(pop3);
                bw.write(sb.toString());
                bw.newLine();
                sb=new StringBuilder();
                sb.append(temp.get(0)).append("\t").append(binMiddle).append("\t").append("piRatio").append("\t");
                sb.append(piRatio2).append("\t").append(p2File2).append("\t").append(pop3);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void transformFstToPlotFormat(List<File> fstFiles, String outDir){
        Predicate<File> d=f->PStringUtils.fastSplit(f.getName(), "_chr").get(1).substring(1,2).equals("D");
        List<File> abFiles=fstFiles.stream().filter(d.negate()).collect(Collectors.toList());
        List<File> dFiles=fstFiles.stream().filter(d).collect(Collectors.toList());
        List<String> abFileNames= abFiles.stream().map(File::getName)
                .map(s->PStringUtils.fastSplit(s, ".windowed").get(0)).collect(Collectors.toList());
        List<String> dFileNames= dFiles.stream().map(File::getName)
                .map(s->PStringUtils.fastSplit(s, ".windowed").get(0)).collect(Collectors.toList());
        List<String> abChr= WheatLineage.abLineage();
        List<String> dChr=WheatLineage.valueOf("D").getChr();
        int indexA, indexB, index1, index2;
        String str1, str2, outFileName;
        for (int i = 0; i < abChr.size(); i++) {
            for (int j = 0; j < p3_aoyue.length-1; j++) {
                str1=p2[0]+"_VS_"+p3_aoyue[j]+"_chr"+abChr.get(i);
                str2=p3_aoyue[j]+"_VS_"+p2[0]+"_chr"+abChr.get(i);
                index1=abFileNames.indexOf(str1);
                index2=abFileNames.indexOf(str2);
                indexA= index1 > -1 ? index1 : index2;
                str1=p2[1]+"_VS_"+p3_aoyue[j]+"_chr"+abChr.get(i);
                str2=p3_aoyue[j]+"_VS_"+p2[1]+"_chr"+abChr.get(i);
                index1=abFileNames.indexOf(str1);
                index2=abFileNames.indexOf(str2);
                indexB= index1 > -1 ? index1 : index2;
                outFileName="chr"+abChr.get(i)+"_"+p3_aoyue[j]+".txt";
                mergeFstFiles(abFiles.get(indexA), p2[0], abFiles.get(indexB), p2[1], new File(outDir, outFileName), p3_aoyue[j]);
            }
        }
        for (int i = 0; i < dChr.size(); i++) {
            str1=p2[0]+"_VS_Ae.tauschii_chr"+dChr.get(i);
            str2="Ae.tauschii_VS_"+ p2[0]+"_chr"+dChr.get(i);
            index1=dFileNames.indexOf(str1);
            index2=dFileNames.indexOf(str2);
            indexA=index1 > -1 ? index1 : index2;
            str1=p2[1]+"_VS_Ae.tauschii_chr"+dChr.get(i);
            str2="Ae.tauschii_VS_"+ p2[1]+"_chr"+dChr.get(i);
            index1=dFileNames.indexOf(str1);
            index2=dFileNames.indexOf(str2);
            indexB=index1 > -1 ? index1 : index2;
            outFileName="chr"+dChr.get(i)+"_Ae.tauschii+.txt";
            mergeFstFiles(dFiles.get(indexA), p2[0], dFiles.get(indexB), p2[1], new File(outDir, outFileName), "Ae.tauschii");
        }
    }

    private static void mergeFstFiles(File fstFile1, String p2File1, File fstFile2, String p2File2, File outFile,
                                      String pop){
        try (BufferedReader br1 = IOTool.getReader(fstFile1);
             BufferedReader br2 =IOTool.getReader(fstFile2);
             BufferedWriter bw=IOTool.getWriter(outFile)) {
            br1.readLine();
            br2.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            List<String> temp;
            int binMiddle=-1;
            StringBuilder sb;
            while ((line=br1.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                binMiddle=(Integer.parseInt(temp.get(1))+Integer.parseInt(temp.get(2)))/2;
                sb=new StringBuilder();
                sb.append(temp.get(0)).append("\t").append(binMiddle).append("\t").append("fst");
                sb.append("\t").append(temp.get(4)).append("\t").append(p2File1).append("\t").append(pop);
                bw.write(sb.toString());
                bw.newLine();
            }
            while ((line=br2.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                binMiddle=(Integer.parseInt(temp.get(1))+Integer.parseInt(temp.get(2)))/2;
                sb=new StringBuilder();
                sb.append(temp.get(0)).append("\t").append(binMiddle).append("\t").append("fst");
                sb.append("\t").append(temp.get(4)).append("\t").append(p2File2).append("\t").append(pop);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void transformFdToPlotFormat(List<File> fdFiles, String outDir){
        Predicate<File> d=f->f.getName().substring(4,5).equals("D");
        List<File> abFiles=fdFiles.stream().filter(d.negate()).collect(Collectors.toList());
        List<File> dFiles=fdFiles.stream().filter(d).collect(Collectors.toList());
        Comparator<File> c=Comparator.comparing(f->f.getName());
        Collections.sort(abFiles, c);
        Collections.sort(dFiles, c);
        String[] p2Names={"Cultivar","Landrace_Europe"};
        String[] p3Names={"Domesticated_emmer","Free_threshing_tetraploid","Wild_emmer"};
        String outName;
        for (int i = 0; i < abFiles.size(); i=i+6) {
            outName=abFiles.get(i).getName().substring(0,6);
            mergeFdFiles(abFiles.get(i), p2Names[0], abFiles.get(i+3), p2Names[1], new File(outDir,
                    outName+p3Names[0]+".txt"), p3Names[0]);
            mergeFdFiles(abFiles.get(i+1), p2Names[0], abFiles.get(i+4), p2Names[1], new File(outDir,
                    outName+p3Names[1]+".txt"), p3Names[1]);
            mergeFdFiles(abFiles.get(i+2),p2Names[0], abFiles.get(i+5), p2Names[1], new File(outDir,
                    outName+p3Names[2]+ ".txt"), p3Names[2]);
        }
        for (int i = 0; i < dFiles.size(); i=i+2) {
            outName=dFiles.get(i).getName().substring(0,6);
            mergeFdFiles(dFiles.get(i), p2Names[0], dFiles.get(i+1), p2Names[1], new File(outDir, outName+"Ae" +
                    ".tauschii.txt"), p3_aoyue[3]);
        }
    }

    private static void mergeFdFiles(File fdFile1, String p2File1, File fdFile2, String p2File2, File outFile,
                                     String pop){
        try (BufferedReader br1 = IOTool.getReader(fdFile1);
             BufferedReader br2 =IOTool.getReader(fdFile2);
             BufferedWriter bw=IOTool.getWriter(outFile)) {
            br1.readLine();
            br2.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            List<String> temp;
            double d, fd;
            int binMiddle=-1;
            StringBuilder sb;
            while ((line=br1.readLine())!=null){
                temp=PStringUtils.fastSplit(line, ",");
                if (temp.get(8).equals("nan") || temp.get(9).equals("nan")) continue;
                d=Double.parseDouble(temp.get(8));
                fd=Double.parseDouble(temp.get(9));
                if (d < 0 || fd <0 || fd >1 ) continue;
                binMiddle=Integer.parseInt(temp.get(3));
                sb=new StringBuilder();
                sb.append(temp.get(0)).append("\t").append(binMiddle).append("\t").append("fd");
                sb.append("\t").append(temp.get(9)).append("\t").append(p2File1).append("\t").append(pop);
                bw.write(sb.toString());
                bw.newLine();
            }
            while ((line=br2.readLine())!=null){
                temp=PStringUtils.fastSplit(line, ",");
                if (temp.get(8).equals("nan") || temp.get(9).equals("nan")) continue;
                d=Double.parseDouble(temp.get(8));
                fd=Double.parseDouble(temp.get(9));
                if (d < 0 || fd <0 || fd >1 ) continue;
                binMiddle=Integer.parseInt(temp.get(3));
                sb=new StringBuilder();
                sb.append(temp.get(0)).append("\t").append(binMiddle).append("\t").append("fd");
                sb.append("\t").append(temp.get(9)).append("\t").append(p2File2).append("\t").append(pop);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private static Predicate<File> getFdPredict(){
        List<String> p2p3=new ArrayList<>();
        for (int i = 0; i < p2.length; i++) {
            for (int j = 0; j < p3_daxing.length; j++) {
                p2p3.add(p2[i]+"_"+p3_daxing[j]+".csv");
            }
        }
        Predicate<File> p=f->p2p3.contains(PStringUtils.fastSplit(f.getName(),"geno_").get(1));
        return p;
    }

    private static Predicate<File> getPiPredict(){
        List<String> piFiles = new ArrayList<>();
        piFiles.addAll(Arrays.asList(p2));
        piFiles.addAll(Arrays.asList(p3_aoyue));
        Predicate<File> p=f->piFiles.contains(PStringUtils.fastSplit(f.getName(),"_chr").get(0));
        return p;
    }

    private static Predicate<File> getFstPredict(){
        List<String> p2p3=new ArrayList<>();
        for (int i = 0; i < p2.length; i++) {
            for (int j = 0; j < p3_aoyue.length; j++) {
                p2p3.add(p2[i]+"_VS_"+p3_aoyue[j]);
                p2p3.add(p3_aoyue[j]+"_VS_"+p2[i]);
            }
        }
        Predicate<File> p=f->p2p3.contains(PStringUtils.fastSplit(f.getName(),"_chr").get(0));
        return p;
    }
}
