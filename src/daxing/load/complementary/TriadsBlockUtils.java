package daxing.load.complementary;

import daxing.common.IOTool;
import daxing.common.PGF;
import daxing.common.StringTool;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class TriadsBlockUtils {

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
                for (int j = 0; j < geneName.length; j++) {
                    refChr=geneName[j].substring(7,9);
                    chrIDArray= RefV1Utils.getChrIDsOfSubgenome(refChr);
                    chrGeneNameList=pgf.getGeneNameOnChr(chrIDArray);
                    geneIndex= Collections.binarySearch(chrGeneNameList, geneName[j]);
                    maxEndIndex=chrGeneNameList.size()-1;
                    startIndex= geneIndex < blockGeneNum ? 0 : geneIndex - blockGeneNum;
                    endIndex= (maxEndIndex - geneIndex) < blockGeneNum ? maxEndIndex : geneIndex + blockGeneNum;
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
}
