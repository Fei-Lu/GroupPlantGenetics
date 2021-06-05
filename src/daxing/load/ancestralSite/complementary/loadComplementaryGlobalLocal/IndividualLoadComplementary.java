package daxing.load.ancestralSite.complementary.loadComplementaryGlobalLocal;

import daxing.common.*;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * 在六倍体retainTriad文件中，使用triadA pos进行滑窗，查看亚基因组load的互补
 * load包括杂合子
 */
public class IndividualLoadComplementary {

    public static void start(){
//        String triadInputDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/006_test_C1_C2_C11_C12/003_triad";
        String triadInputDir="/Users/xudaxing/Desktop/IndivComplementary2/001_triad";
        String pgfFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/001_vmap2.1Before20200525/002_analysis/014_deleterious/wheat_v1.1_Lulab_geneHC.pgf";
        String triadGeneFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/001_vmap2.1Before20200525/002_analysis/014_deleterious/triadGenes1.1_cdsLen_geneHC.txt";
//        String outDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/006_test_C1_C2_C11_C12/004_nonmalizedAndSlidingWindow";
        String outDir="/Users/xudaxing/Desktop/IndivComplementary2/002_normalizedAndSlidingWindow";
        String[] subDirs={"001_triad.pos","002_forSlidingWindow","003_slidingWindow","004_mergeIndividualDel"};
        String sub="A";
        File[] subDirFiles=new File[subDirs.length];
        for (int i = 0; i < subDirs.length; i++) {
            subDirFiles[i]=new File(outDir, subDirs[i]);
            subDirFiles[i].mkdir();
        }
        addTriadPos(triadInputDir, pgfFile, subDirFiles[0].getAbsolutePath());
        forSlidingWindow(subDirFiles[0].getAbsolutePath(), triadGeneFile, subDirFiles[1].getAbsolutePath());
        slidingWindow(subDirFiles[1].getAbsolutePath(),sub, 10_000_000, 1_000_000, subDirFiles[2].getAbsolutePath());
        mergeTaxonDel(subDirFiles[2].getAbsolutePath(),
                new File(subDirFiles[3],"triadPos"+sub+".IndividualMerged.txt.gz").getAbsolutePath());
    }

    public static void addTriadPos(String inputDir, String pgfFile, String outDir){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        PGF pgf=new PGF(pgfFile);
        String[] outName= files.stream().map(File::getName).map(s->s.replaceAll(".txt", "pos.txt")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->addTriadPos(files.get(e), pgf, new File(outDir, outName[e])));
    }

    /**
     * triadPos 为基因中点
     * @param inputFile
     * @param pgf
     * @param outFile
     */
    private static void addTriadPos(File inputFile, PGF pgf, File outFile){
        pgf.sortGeneByName();
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
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

    public static void forSlidingWindow(String inputDir, String triadGeneFile, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        String[] outNamesA= files.stream().map(File::getName).
                map(s -> s.replaceAll(".txt.gz","A.txt.gz")).toArray(String[]::new);
        String[] outNamesB= files.stream().map(File::getName).
                map(s -> s.replaceAll(".txt.gz","B.txt.gz")).toArray(String[]::new);
        String[] outNamesD= files.stream().map(File::getName).
                map(s -> s.replaceAll(".txt.gz","D.txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->forSlidingWindow(files.get(e), "A", triadGeneFile, new File(outDir
                ,outNamesA[e])));
        IntStream.range(0, files.size()).forEach(e->forSlidingWindow(files.get(e), "B", triadGeneFile, new File(outDir
                ,outNamesB[e])));
        IntStream.range(0, files.size()).forEach(e->forSlidingWindow(files.get(e), "D", triadGeneFile, new File(outDir
                ,outNamesD[e])));
    }

    public static void forSlidingWindow(File inputFile, String triadPosSubgenome, String triadGeneFile, File outFile){
        Triads triads =new Triads(triadGeneFile);
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            br.readLine();
            String header="TriadsPosChr\tPos\tNormalizedDerivedSyn\tNormalizedDerivedNon\tNormalizedDerivedDel" +
                    "\tTriadID\tChr";
            bw.write(header);
            bw.newLine();
            String chr, triadsPosChr, pos, triadID;
            double normalizedDerivedSyn, normalizedDerivedNon,normalizedDerivedDel;
            double cdsLen, ancestralNum, derivedSynNum, derivedNonsynNum, derivedDelNum;
            double heterSynNum, heterNonNum, heterDelNum;
            String[] subs={"A","B","D"};
            int indexSub= Arrays.binarySearch(subs, triadPosSubgenome);
            StringBuilder sb=new StringBuilder();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                triadID=temp.get(0);
                chr=temp.get(2).substring(7,9);
                triadsPosChr=getTriadPosChr(triads.getTraidGenes(triadID),triadPosSubgenome, chr);
                pos=temp.get(4+indexSub);
                cdsLen=Integer.parseInt(temp.get(1));
                ancestralNum=Integer.parseInt(temp.get(7))+Integer.parseInt(temp.get(10));
                sb.setLength(0);
                sb.append(triadsPosChr).append("\t").append(pos).append("\t");
                derivedSynNum=Integer.parseInt(temp.get(8));
                heterSynNum=Integer.parseInt(temp.get(9));
                derivedNonsynNum=Integer.parseInt(temp.get(11));
                heterNonNum=Integer.parseInt(temp.get(12));
                derivedDelNum=Integer.parseInt(temp.get(14));
                heterDelNum=Integer.parseInt(temp.get(15));

//                normalizedDerivedSyn=(derivedSynNum+heterSynNum*0.5)*10000/(ancestralNum*cdsLen);
//                normalizedDerivedNon=(derivedNonsynNum+heterNonNum*0.5)*10000/(ancestralNum*cdsLen);
//                normalizedDerivedDel=(derivedDelNum+heterDelNum*0.5)*10000/(ancestralNum*cdsLen);

                normalizedDerivedSyn=(derivedSynNum);
                normalizedDerivedNon=(derivedNonsynNum);
                normalizedDerivedDel=(derivedDelNum);

                if (ancestralNum==0){
                    sb.append("NaN").append("\t").append("NaN").append("\t").append("NaN").append("\t");
                }else {
                    sb.append(NumberTool.format(normalizedDerivedSyn, 5)).append("\t");
                    sb.append(NumberTool.format(normalizedDerivedNon,5)).append("\t");
                    sb.append(NumberTool.format(normalizedDerivedDel, 5)).append("\t");
                }
                sb.append(triadID).append("\t").append(chr);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static String getTriadPosChr(String[] triadsABD, String triadPosSubgenome, String chr){
        String[] abd={"A","B","D"};
        int index=Arrays.binarySearch(abd, triadPosSubgenome);
        String chrTriad=triadsABD[index];
        String subID=chrTriad.substring(7,8);
        String sub=chr.substring(1,2);
        StringBuilder sb=new StringBuilder();
        sb.append(subID).append(sub);
        return sb.toString();
    }

    /**
     *
     * @param inputDir
     * @param triadPosition A B D
     * @param windowSize
     * @param stepSize
     * @param outDir
     */
    public static void slidingWindow(String inputDir, String triadPosition, int windowSize, int stepSize,
                                     String outDir){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        Predicate<File> filesAP=file -> file.getName().contains("triadpos"+triadPosition);
        List<File> filesA=files.stream().filter(filesAP).collect(Collectors.toList());
        String[] outNames=filesA.stream().map(File::getName).map(s -> s.replaceAll("txt.gz",
                windowSize/1000000+"MbWindow_"+stepSize/1000000)+"MbStep.txt.gz").toArray(String[]::new);
        IntStream.range(0, filesA.size()).forEach(e->slidingWindow(filesA.get(e),triadPosition, windowSize, stepSize,
         new File(outDir, outNames[e])));
    }

    public static void slidingWindow(File inputFile, String triadPosition, int windowSize, int stepSize,
                                     File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            int chrNum= RowTableTool.getColumnSet(inputFile.getAbsolutePath(),0).size();
            String line, chr;
            List<String> temp;
            int pos, chrIndex;
            double[] synNonDelValue;
            String[] chrs=RefV1Utils.getChromosomes();
            int[] chrASize=new int[chrs.length];
            for (int i = 0; i <chrs.length; i=i+3) {
                chrASize[i]=RefV1Utils.getChromosomeLength(chrs[i]);
                chrASize[i+1]=RefV1Utils.getChromosomeLength(chrs[i]);
                chrASize[i+2]=RefV1Utils.getChromosomeLength(chrs[i]);
            }
            SlidingWindowForLoadComplement[] slidingWindows=new SlidingWindowForLoadComplement[chrNum];
            for (int i = 0; i < slidingWindows.length; i++) {
                slidingWindows[i]=new SlidingWindowForLoadComplement(chrASize[i],windowSize, stepSize);
            }
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chr=temp.get(0);
                chrIndex= Arrays.binarySearch(chrs, chr);
                pos=Integer.parseInt(temp.get(1));
                if (temp.get(2).equals("NaN")) continue;
                synNonDelValue=new double[3];
                synNonDelValue[0]=Double.parseDouble(temp.get(2));
                synNonDelValue[1]=Double.parseDouble(temp.get(3));
                synNonDelValue[2]=Double.parseDouble(temp.get(4));
                slidingWindows[chrIndex].addValue(pos, synNonDelValue);
            }
            bw.write("TriadsPosChr\tWindow_Start\tWindowEnd\tGeneNum\tMean_NormalizedSynLoad" +
                    "\tMean_NormalizedNonsynLoad" +
                    "\tMean_NormalizedDeleteriousLoad");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            int windowNum, windowStart, windowEnd, geneNum;
            double mean_NormalizedSynLoad, mean_NormalizedNonLoad, mean_NormalizedDelLoad;
            TDoubleArrayList synValueList, nonValueList, delValueList;
            for (int i = 0; i < slidingWindows.length; i++) {
                windowNum=slidingWindows[i].getWindowNum();
                chr=chrs[i];
                for (int j = 0; j < windowNum; j++) {
                    windowStart=slidingWindows[i].getWindowStarts()[j];
                    windowEnd=slidingWindows[i].getWindowEnds()[j];
                    geneNum=slidingWindows[i].getCountInWindow(j);
                    synValueList=slidingWindows[i].getWindow1Value(j);
                    nonValueList=slidingWindows[i].getWindow2Value(j);
                    delValueList=slidingWindows[i].getWindow3Value(j);
                    mean_NormalizedSynLoad=synValueList.sum()/geneNum;
                    mean_NormalizedNonLoad=nonValueList.sum()/geneNum;
                    mean_NormalizedDelLoad=delValueList.sum()/geneNum;
                    sb.setLength(0);
                    sb.append(chr).append("\t");
                    sb.append(windowStart).append("\t").append(windowEnd).append("\t");
                    sb.append(geneNum).append("\t");
                    if (Double.isNaN(mean_NormalizedSynLoad)){
                        sb.append("NaN").append("\t");
                    }else {
                        sb.append(NumberTool.format(mean_NormalizedSynLoad,6)).append("\t");
                    }
                    if (Double.isNaN(mean_NormalizedNonLoad)){
                        sb.append("NaN").append("\t");
                    }else {
                        sb.append(NumberTool.format(mean_NormalizedNonLoad, 6)).append("\t");
                    }
                    if (Double.isNaN(mean_NormalizedDelLoad)){
                        sb.append("NaN");
                    }else {
                        sb.append(NumberTool.format(mean_NormalizedDelLoad, 6));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void mergeTaxonDel(String inputDir, String outFile){
        List<File> files=IOTool.getVisibleDir(inputDir);
        RowTableTool<String> table0=new RowTableTool<>(files.get(0).getAbsolutePath());
        List<String> meanSynCount=table0.getColumn(4);
        List<String> meanNonCount=table0.getColumn(5);
        List<String> meanDelCount=table0.getColumn(6);
        table0.removeColumn("GeneNum");
        table0.removeColumn("Mean_NormalizedSynLoad");
        table0.removeColumn("Mean_NormalizedNonsynLoad");
        table0.removeColumn("Mean_NormalizedDeleteriousLoad");
        RowTableTool<String> table1;
        String taxonName=PStringUtils.fastSplit(files.get(0).getName(), ".").get(0);
//        String columnName=taxonName+".MeanNormalizedDeleteriousLoad";

        table0.addColumn(taxonName+".MeanSynDerivedCount", meanSynCount);
        table0.addColumn(taxonName+".MeanNonDerivedCount", meanNonCount);
        table0.addColumn(taxonName+".MeanDelDerivedCount", meanDelCount);

//        for (int i = 1; i < files.size(); i++) {
//            table1=new RowTableTool<>(files.get(i).getAbsolutePath());
//            taxonName=PStringUtils.fastSplit(files.get(i).getName(), ".").get(0);
//            table0.addColumn(taxonName+".MeanNormalizedDeleteriousLoad", table1.getColumn(6));
//        }

        for (int i = 1; i < files.size(); i++) {
            table1=new RowTableTool<>(files.get(i).getAbsolutePath());
            taxonName=PStringUtils.fastSplit(files.get(i).getName(), ".").get(0);
            table0.addColumn(taxonName+".MeanSynDerivedCount", table1.getColumn(4));
            table0.addColumn(taxonName+".MeanNonDerivedCount", table1.getColumn(5));
            table0.addColumn(taxonName+".MeanDelDerivedCount", table1.getColumn(6));
        }

        table0.write(outFile, IOFileFormat.TextGzip);
    }

}
