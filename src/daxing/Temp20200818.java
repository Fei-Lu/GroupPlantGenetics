package daxing;

import daxing.common.*;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 * 个体load来源于毕傲月(包括杂合基因型load)
 */
public class Temp20200818 {

//    public static void main(String[] args){
//        String inputDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/003_retainTriad/005_hexaploid_test";
//        String pgfFile="/Users/xudaxing/Documents/deleteriousMutation/001_vmap2.1Before20200525/002_analysis/014_deleterious/wheat_v1.1_Lulab_geneHC.pgf";
//        String outDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/003_retainTriad/006_hexaploid_test_expected";
//        addTriadPos(inputDir, pgfFile, outDir);

//        String inputDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_retainTriad/006_hexaploid_test_expected";
//        String outDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_retainTriad/009_retainTriadHexaploid_expected_pos_slidingWindow/001_triadPos";
//        forSlidingWindow(inputDir, outDir);
//
//        String inputDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_retainTriad/009_retainTriadHexaploid_expected_pos_slidingWindow/001_triadPos";
//        String outDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_retainTriad/009_retainTriadHexaploid_expected_pos_slidingWindow/002_slidWindow";
//        String synNonDel="syn";
//        splitChrAndSlidWindow(inputDir,outDir,"syn");
//        splitChrAndSlidWindow(inputDir,outDir,"non");
//        splitChrAndSlidWindow(inputDir,outDir,"del");
//    }

    public static void splitChrAndSlidWindow(String inputDir, String outDir, String synNonDel){
        String[] synNonDelArray={"syn","non","del"};
        int columnIndex= Arrays.asList(synNonDelArray).indexOf(synNonDel)+2;
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        String[] outFileNames= files.stream().map(File::getName).
                map(s->s.replaceAll("txt.gz",synNonDel+".10M_slidWindow_1M_step.txt.gz")).toArray(String[]::new);
        String[] outDirNames=files.stream().map(File::getName).map(s->s.replaceAll(".txt.gz",".splitChr.txt.gz")).toArray(String[]::new);
        File[] subDir=new File[outDirNames.length];
        File[] subSlidWindowDir=new File[outDirNames.length];
        for (int i = 0; i < outDirNames.length; i++) {
            subDir[i]=new File(outDir, outDirNames[i]);
            subDir[i].mkdir();
            subSlidWindowDir[i]=new File(outDir, outDirNames[i]+"_slidWindow");
            subSlidWindowDir[i].mkdir();
        }
        VCF vcf;
        for (int i = 0; i < files.size(); i++) {
            vcf=new VCF(files.get(i));
            vcf.writeVcfToSplitedChr(subDir[i].getAbsolutePath());
        }
        for (int i = 0; i < subDir.length; i++) {
            PlotTools.slidingWindow(subDir[i].getAbsolutePath(),subSlidWindowDir[i].getAbsolutePath(),columnIndex);
        }
        for (int i = 0; i < subSlidWindowDir.length; i++) {
            RowTableTool.meregMultipleTable(subSlidWindowDir[i].getAbsolutePath(), "Chr",new File(outDir,
                    outFileNames[i]).getAbsolutePath());
        }
        for (int i = 0; i < subDir.length; i++) {
            for(File file: subDir[i].listFiles()){
                if (!file.isDirectory()){
                    file.delete();
                }
            }
            subDir[i].delete();
            for (File file: subSlidWindowDir[i].listFiles()){
                if (!file.isDirectory()){
                    file.delete();
                }
            }
            subSlidWindowDir[i].delete();
        }

    }

    public static void forSlidingWindow(String inputDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        String[] outNamesA= files.stream().map(File::getName).
                map(s -> s.replaceAll(".txt.gz","A.txt.gz")).toArray(String[]::new);
        String[] outNamesB= files.stream().map(File::getName).
                map(s -> s.replaceAll(".txt.gz","B.txt.gz")).toArray(String[]::new);
        String[] outNamesD= files.stream().map(File::getName).
                map(s -> s.replaceAll(".txt.gz","D.txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->forSlidingWindow(files.get(e), "A", new File(outDir,outNamesA[e])));
        IntStream.range(0, files.size()).forEach(e->forSlidingWindow(files.get(e), "B", new File(outDir,outNamesB[e])));
        IntStream.range(0, files.size()).forEach(e->forSlidingWindow(files.get(e), "D", new File(outDir,outNamesD[e])));
    }

    public static void forSlidingWindow(File inputFile, String triadPosSubgenome, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getTextWriter(outFile)) {
            String line;
            List<String> temp;
            br.readLine();
            String header="Chr\tPos\tNormalizedDerivedSyn\tNormalizedDerivedNon\tNormalizedDerivedDel\tTriadID";
            bw.write(header);
            bw.newLine();
            String chr, pos, triadID;
            double normalizedDerivedSyn, normalizedDerivedNon,normalizedDerivedDel;
//            1000*10*derivedSyn[i]/(cdsLen[i]*snpNum[i]);
            double cdsLen, ancestralNum, derivedSynNum, derivedNonsynNum, derivedDelNum;
            String[] subs={"A","B","D"};
            int indexSub= Arrays.binarySearch(subs, triadPosSubgenome);
            StringBuilder sb=new StringBuilder();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chr=temp.get(2).substring(7,9);
                pos=temp.get(4+indexSub);
                cdsLen=Integer.parseInt(temp.get(1));
                ancestralNum=Integer.parseInt(temp.get(14));
                if (ancestralNum==0){
                    continue;
                }
                derivedSynNum=Integer.parseInt(temp.get(8));
                derivedNonsynNum=Integer.parseInt(temp.get(10));
                derivedDelNum=Integer.parseInt(temp.get(12));
                normalizedDerivedSyn=derivedSynNum*10000/(ancestralNum*cdsLen);
                normalizedDerivedNon=derivedNonsynNum*10000/(ancestralNum*cdsLen);
                normalizedDerivedDel=derivedDelNum*10000/(ancestralNum*cdsLen);
                triadID=temp.get(0);
                sb.setLength(0);
                sb.append(chr).append("\t").append(pos).append("\t");
                sb.append(NumberTool.format(normalizedDerivedSyn, 5)).append("\t");
                sb.append(NumberTool.format(normalizedDerivedNon,5)).append("\t");
                sb.append(NumberTool.format(normalizedDerivedDel, 5)).append("\t");
                sb.append(triadID);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void addTriadPos(String inputDir, String pgfFile, String outDir){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        PGF pgf=new PGF(pgfFile);
        String[] outName= files.stream().map(File::getName).map(s->s.replaceAll(".txt", "pos.txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->addTriadPos(files.get(e), pgf, new File(outDir, outName[e])));
    }

    private static void addTriadPos(File inputFile, PGF pgf, File outFile){
        pgf.sortGeneByName();
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getTextWriter(outFile)) {
            String line;
            String[] geneNames;
            List<String> temp;
            List<String>[] tem;
            temp= PStringUtils.fastSplit(br.readLine());
            temp.add(3,"Pos\tTriadPosA\tTriadPosB\tTriadPosD");
            bw.write(String.join("\t",temp));
            bw.newLine();
            int[] geneIndex, posOnChrID, posOnChr, chrIDStart, chrIDEnd, chrID;
            int[] pos;
            StringBuilder sbA=new StringBuilder();
            StringBuilder sbB=new StringBuilder();
            StringBuilder sbD=new StringBuilder();
            while ((line=br.readLine())!=null){
                tem=new List[3];
                tem[0]=PStringUtils.fastSplit(line);
                tem[1]=PStringUtils.fastSplit(br.readLine());
                tem[2]=PStringUtils.fastSplit(br.readLine());
                geneNames=new String[3];
                for (int i = 0; i < 3; i++) {
                    geneNames[i]=tem[i].get(2);
                }
                geneIndex=new int[3];
                for (int i = 0; i < 3; i++) {
                    geneIndex[i]=pgf.getGeneIndex(geneNames[i]);
                }
                chrID = new int[3];
                chrIDStart=new int[3];
                chrIDEnd=new int[3];
                posOnChrID=new int[3];
                posOnChr=new int[3];
                pos=new int[3];
                for (int i = 0; i < 3; i++) {
                    chrID[i]=pgf.getGene(geneIndex[i]).getGeneRange().chr;
                    chrIDStart[i]=pgf.getGene(geneIndex[i]).getGeneRange().start;
                    chrIDEnd[i]=pgf.getGene(geneIndex[i]).getGeneRange().end;
                    posOnChrID[i]=(chrIDStart[i]+chrIDEnd[i])/2;
                    posOnChr[i]= RefV1Utils.getPosOnChromosome(chrID[i], posOnChrID[i]);
                    pos[i]=posOnChr[i];
                }
                sbA.setLength(0);
                sbB.setLength(0);
                sbD.setLength(0);
//                temp.add(3,"Pos\tTriadPosA\tTriadPosB\tTriadPosD");
                sbA.append(pos[0]).append("\t").append(pos[0]).append("\t").append(pos[1]).append("\t").append(pos[2]);
                sbB.append(pos[1]).append("\t").append(pos[0]).append("\t").append(pos[1]).append("\t").append(pos[2]);
                sbD.append(pos[2]).append("\t").append(pos[0]).append("\t").append(pos[1]).append("\t").append(pos[2]);
                tem[0].add(3, sbA.toString());
                tem[1].add(3, sbB.toString());
                tem[2].add(3, sbD.toString());
                for (int i = 0; i < 3; i++) {
                    bw.write(String.join("\t",tem[i]));
                    bw.newLine();
                }
                System.out.println();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
