package daxing.load.complementary;

import daxing.common.*;
import daxing.load.complementary.selection.TriadsBlockScore;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class TriadsBlockUtils {

    /**
     *
     * @param triadGeneFile
     * @param pgfFile
     * @param blockGeneNum 上(下)游基因数目, total: blockGeneNum * 2 +1
     * @param outFile
     */
    public static void writeTriadsBlock(String triadGeneFile, String pgfFile, int blockGeneNum, String outFile){
        PGF pgf=new PGF(pgfFile);
        try (BufferedReader br = IOTool.getReader(triadGeneFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            bw.write("TriadID\tTriadsGeneName\tCDSLen\tBlockGeneNum\tBlockGeneListA\tBlockGeneListB\tBlockGeneListD");
            bw.newLine();
            br.readLine();
            String line, triadsID, refChr;
            String[] geneName;
            int[] cdsLen, chrIDArray;
            List<String> temp, chrGeneNameList;
            TriadsBlock triadsBlock=null;
            List<String>[] blockGeneName;
            int geneIndex, startIndex, endIndex, maxEndIndex, triadsBlockNum=0;
            Map<String,Integer> geneNameIndexMap;
            System.out.println("start writing triadsBlock...");
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                if (!StringTool.isNumeric(temp.get(6))) continue;
                if (!StringTool.isNumeric(temp.get(7))) continue;
                if (!StringTool.isNumeric(temp.get(8))) continue;
                triadsID=temp.get(0);
                geneName=new String[3];
                cdsLen=new int[3];
                for (int i = 0; i < geneName.length; i++) {
                    geneName[i]=temp.get(1+i);
                    cdsLen[i]=Integer.parseInt(temp.get(6+i));
                }
                blockGeneName=new List[3];
                for (int j = 0; j < blockGeneName.length; j++) {
                    blockGeneName[j]=new ArrayList<>();
                }
                geneNameIndexMap=new HashMap<>();
                for (int j = 0; j < geneName.length; j++) {
                    refChr=geneName[j].substring(7,9);
                    chrIDArray= RefV1Utils.getChrIDsOfSubgenome(refChr);
                    chrGeneNameList=pgf.getGeneNameOnChr(chrIDArray);
                    for (int i = 0; i < chrGeneNameList.size(); i++) {
                        geneNameIndexMap.put(chrGeneNameList.get(i), i);
                    }
                    geneIndex=geneNameIndexMap.get(geneName[j]);
                    maxEndIndex=chrGeneNameList.size()-1;
                    startIndex= geneIndex < blockGeneNum/2 ? 0 : geneIndex - blockGeneNum/2;
                    endIndex= (maxEndIndex - geneIndex) < blockGeneNum/2 ? maxEndIndex : geneIndex + blockGeneNum/2;
                    for (int k = startIndex; k < endIndex+1; k++) {
                        blockGeneName[j].add(chrGeneNameList.get(k));
                    }
                }
                triadsBlock=new TriadsBlock(triadsID, geneName, cdsLen, blockGeneName);
                bw.write(triadsBlock.toString());
                bw.newLine();
                triadsBlockNum++;
                if (triadsBlockNum%2000==0){
                    System.out.println("writing "+triadsBlockNum+" to "+new File(outFile).getName());
                }
            }
            bw.flush();
            System.out.println("total triadsBlock num is "+triadsBlockNum);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static TriadsBlock[] readTriadBlock(String triadsBlockFile){
        List<TriadsBlock> triadsBlockList=new ArrayList<>();
        TriadsBlock triadsBlock;
        try (BufferedReader br = IOTool.getReader(triadsBlockFile)) {
            br.readLine();
            String line;
            List<String> temp, tem;
            String triadsID;
            String[] geneName;
            int[] triadsCDSLen;
            List<String>[] blockGeneName;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                triadsID=temp.get(0);
                geneName=new String[3];
                triadsCDSLen=new int[3];
                blockGeneName=new List[3];
                for (int i = 0; i < blockGeneName.length; i++) {
                    blockGeneName[i]=new ArrayList<>();
                }
                tem=PStringUtils.fastSplit(temp.get(1), ",");
                for (int i = 0; i < tem.size(); i++) {
                    geneName[i]=tem.get(i);
                }
                tem=PStringUtils.fastSplit(temp.get(2),",");
                for (int i = 0; i < tem.size(); i++) {
                    triadsCDSLen[i]=Integer.parseInt(tem.get(i));
                }
                blockGeneName[0].addAll(PStringUtils.fastSplit(temp.get(4), ","));
                blockGeneName[1].addAll(PStringUtils.fastSplit(temp.get(5), ","));
                blockGeneName[2].addAll(PStringUtils.fastSplit(temp.get(6), ","));
                triadsBlock=new TriadsBlock(triadsID, geneName, triadsCDSLen, blockGeneName);
                triadsBlockList.add(triadsBlock);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        TriadsBlock[] triadsBlocks=new TriadsBlock[triadsBlockList.size()];
        for (int i = 0; i < triadsBlocks.length; i++) {
            triadsBlocks[i]=triadsBlockList.get(i);
        }
        return triadsBlocks;
    }

    /**
     *
     * @param triadsBlockInputFile no ChrRange
     * @param pgfFile
     * @param triadsBlockOutFile have ChrRange
     */
    public static void writeTriadsBlock(String triadsBlockInputFile, String pgfFile, String triadsBlockOutFile){
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByName();
        TriadsBlock[] triadsBlockArray= TriadsBlockUtils.readTriadBlock(triadsBlockInputFile);
        try (BufferedWriter bw = IOTool.getWriter(triadsBlockOutFile)) {
            StringBuilder sb=new StringBuilder();
            sb.setLength(0);
            sb.append("TriadID\tTriadsGeneName\tCDSLen\tBlockGeneNum\tBlockGeneListA\tBlockGeneListB\tBlockGeneListD");
            sb.append("\t").append("ChrRange");
            bw.write(sb.toString());
            bw.newLine();
            ChrRange[] chrRanges;
            for (int i = 0; i < triadsBlockArray.length; i++) {
                sb.setLength(0);
                chrRanges=TriadsBlock.getChrRange(triadsBlockArray[i],pgf);
                sb.append(triadsBlockArray[i].toString()).append("\t");
                for (int j = 0; j < chrRanges.length; j++) {
                    sb.append(chrRanges[j].toString()).append(";");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static TriadsBlock[] readFromTriadsBlockChrRange(String triadsBlockChrRangeFile){
        List<TriadsBlock> triadsBlockList=new ArrayList<>();
        TriadsBlock triadsBlock;
        ChrRange[] chrRanges;
        ChrRange chrRange;
        try (BufferedReader br = IOTool.getReader(triadsBlockChrRangeFile)) {
            br.readLine();
            String line, chr;
            int start, end;
            List<String> temp, tem, te, t;
            String triadsID;
            String[] geneName;
            int[] triadsCDSLen;
            List<String>[] blockGeneName;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                triadsID=temp.get(0);
                geneName=new String[3];
                triadsCDSLen=new int[3];
                blockGeneName=new List[3];
                for (int i = 0; i < blockGeneName.length; i++) {
                    blockGeneName[i]=new ArrayList<>();
                }
                tem=PStringUtils.fastSplit(temp.get(1), ",");
                for (int i = 0; i < tem.size(); i++) {
                    geneName[i]=tem.get(i);
                }
                tem=PStringUtils.fastSplit(temp.get(2),",");
                for (int i = 0; i < tem.size(); i++) {
                    triadsCDSLen[i]=Integer.parseInt(tem.get(i));
                }
                blockGeneName[0].addAll(PStringUtils.fastSplit(temp.get(4), ","));
                blockGeneName[1].addAll(PStringUtils.fastSplit(temp.get(5), ","));
                blockGeneName[2].addAll(PStringUtils.fastSplit(temp.get(6), ","));
                chrRanges=new ChrRange[WheatLineage.values().length];
                tem=PStringUtils.fastSplit(temp.get(7), ";");
                for (int i = 0; i <tem.size(); i++) {
                    te=PStringUtils.fastSplit(tem.get(i), ":");
                    chr=te.get(0);
                    t=PStringUtils.fastSplit(te.get(1), ",");
                    start=Integer.parseInt(t.get(0));
                    end=Integer.parseInt(t.get(1));
                    chrRange=new ChrRange(chr, start, end);
                    chrRanges[i]=chrRange;
                }
                triadsBlock=new TriadsBlock(triadsID, geneName, triadsCDSLen, blockGeneName, chrRanges);
                triadsBlockList.add(triadsBlock);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        TriadsBlock[] triadsBlocks=new TriadsBlock[triadsBlockList.size()];
        for (int i = 0; i < triadsBlocks.length; i++) {
            triadsBlocks[i]=triadsBlockList.get(i);
        }
        return triadsBlocks;
    }

    public static TriadsBlockScore[] buildFromTriadsBlockChrRange(String triadsBlockChrRangeFile){
        List<TriadsBlockScore> triadsBlockList=new ArrayList<>();
        TriadsBlockScore triadsBlock;
        ChrRange[] chrRanges;
        ChrRange chrRange;
        try (BufferedReader br = IOTool.getReader(triadsBlockChrRangeFile)) {
            br.readLine();
            String line, chr;
            int start, end;
            List<String> temp, tem, te, t;
            String triadsID;
            String[] geneName;
            int[] triadsCDSLen;
            List<String>[] blockGeneName;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                triadsID=temp.get(0);
                geneName=new String[3];
                triadsCDSLen=new int[3];
                blockGeneName=new List[3];
                for (int i = 0; i < blockGeneName.length; i++) {
                    blockGeneName[i]=new ArrayList<>();
                }
                tem=PStringUtils.fastSplit(temp.get(1), ",");
                for (int i = 0; i < tem.size(); i++) {
                    geneName[i]=tem.get(i);
                }
                tem=PStringUtils.fastSplit(temp.get(2),",");
                for (int i = 0; i < tem.size(); i++) {
                    triadsCDSLen[i]=Integer.parseInt(tem.get(i));
                }
                blockGeneName[0].addAll(PStringUtils.fastSplit(temp.get(4), ","));
                blockGeneName[1].addAll(PStringUtils.fastSplit(temp.get(5), ","));
                blockGeneName[2].addAll(PStringUtils.fastSplit(temp.get(6), ","));
                chrRanges=new ChrRange[WheatLineage.values().length];
                tem=PStringUtils.fastSplit(temp.get(7), ";");
                for (int i = 0; i <tem.size(); i++) {
                    te=PStringUtils.fastSplit(tem.get(i), ":");
                    chr=te.get(0);
                    t=PStringUtils.fastSplit(te.get(1), ",");
                    start=Integer.parseInt(t.get(0));
                    end=Integer.parseInt(t.get(1));
                    chrRange=new ChrRange(chr, start, end);
                    chrRanges[i]=chrRange;
                }
                triadsBlock=new TriadsBlockScore(triadsID, geneName, triadsCDSLen, blockGeneName, chrRanges);
                triadsBlockList.add(triadsBlock);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        TriadsBlockScore[] triadsBlockScoreArray=new TriadsBlockScore[triadsBlockList.size()];
        for (int i = 0; i < triadsBlockScoreArray.length; i++) {
            triadsBlockScoreArray[i]=triadsBlockList.get(i);
        }
        return triadsBlockScoreArray;
    }
}
