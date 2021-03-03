package daxing.load.complementary.selection;

import daxing.common.ChrRange;
import daxing.common.IOTool;
import daxing.load.complementary.TriadsBlock;
import daxing.load.complementary.TriadsBlockUtils;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

public class DominanceAdditiveSelection {

    TriadsBlockScore[] triadsBlockScoreArray;

    public DominanceAdditiveSelection(String triadsBlockChrRangeFile, String xpclrScoreFile){
        triadsBlockScoreArray= TriadsBlockUtils.buildFromTriadsBlockChrRange(triadsBlockChrRangeFile);
        initializeXPCLR(xpclrScoreFile);
    }

    public enum Sub {
        A(0), B(1), D(2);

        int index;

        Sub(int i) {
            this.index=i;
        }

        public int getIndex() {
            return index;
        }

        public static Sub newInstanceFromIndex(int index){
            switch (index){
                case 0:
                    return A;
                case 1:
                    return B;
                case 2:
                    return D;
            }
            return null;
        }
    }

    private void initializeXPCLR(String xpclrScoreFile){
        Map<ChrRange, Double>[] chrRangeScoreMapArray=groupABD(xpclrScoreFile);
        Sub sub;
        Comparator<TriadsBlockScore> subComparator;
        Map<ChrRange, Double> chrRangeDoubleMap;
        TIntArrayList indexList;
        for (int i = 0; i < chrRangeScoreMapArray.length; i++) {
            chrRangeDoubleMap=chrRangeScoreMapArray[i];
            sub=Sub.newInstanceFromIndex(i);
            subComparator=getSortComparatorBy(sub);
            Arrays.sort(triadsBlockScoreArray, subComparator);
            for (Map.Entry<ChrRange,Double> entry: chrRangeDoubleMap.entrySet()){
                indexList=this.getChrRangeIndex(entry.getKey(), subComparator);
                for (int j = 0; j < indexList.size(); j++) {
                    this.triadsBlockScoreArray[indexList.get(j)].addXPCLRScore(sub, entry.getValue());
                }
            }
        }
    }

    private Map<ChrRange, Double>[] groupABD(String xpclrScoreFile){
        Map<ChrRange, Double>[] chrRangeScoreMapArray=new Map[Sub.values().length];
        for (int i = 0; i < chrRangeScoreMapArray.length; i++) {
            chrRangeScoreMapArray[i]=new HashMap<>();
        }
        try (BufferedReader br = IOTool.getReader(xpclrScoreFile)) {
            br.readLine();
            String line, refChr;
            int start, end;
            List<String> temp;
            ChrRange chrRange;
            Sub sub;
            double xpclrScore;
            while ((line= br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                refChr=temp.get(5);
                start=Integer.parseInt(temp.get(6));
                end=start+100_000;
                chrRange=new ChrRange(refChr, start, end);
                sub=Sub.valueOf(refChr.substring(1,2));
                xpclrScore=Double.parseDouble(temp.get(7));
                chrRangeScoreMapArray[sub.getIndex()].put(chrRange, xpclrScore);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return chrRangeScoreMapArray;
    }

    /**
     * 根据triads亚基因组的position进行排序
     * @param sub
     */
    private static Comparator<TriadsBlockScore> getSortComparatorBy(Sub sub){
        return Comparator.comparing(triadsBlockScore -> triadsBlockScore.getChrRanges()[sub.getIndex()]);
    }

    private TIntArrayList getChrRangeIndex(ChrRange chrRange, Comparator<TriadsBlockScore> subComparator){
        Sub sub=Sub.valueOf(chrRange.getChr().substring(1,2));
        TriadsBlockScore triadsBlockScore=new TriadsBlockScore(chrRange);
        TIntArrayList indexList=new TIntArrayList();
        int index, indexLeft, indexRight;
        int hit=Arrays.binarySearch(triadsBlockScoreArray, triadsBlockScore, subComparator);
        index = hit < 0 ? -hit-1 : hit;
        if (index < triadsBlockScoreArray.length){
            indexList.add(index);
        }
        indexLeft=index;
        indexRight=index;
        while ((indexLeft > 0)){
            indexLeft--;
            if (!triadsBlockScoreArray[indexLeft].getChrRanges()[sub.getIndex()].isOverlapped(chrRange)) break;
            indexList.add(indexLeft);
        }
        while ((indexRight < triadsBlockScoreArray.length-1)){
            indexRight++;
            if (!triadsBlockScoreArray[indexRight].getChrRanges()[sub.getIndex()].isOverlapped(chrRange)) break;
            indexList.add(indexRight);
        }
        indexList.sort();
        return indexList;
    }

    public void writeTriadsBlockScore(String outFile){
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            sb.append("TriadsBlockID\tSub\tXPCLR_Mean");
            bw.write(sb.toString());
            bw.newLine();
            Comparator<TriadsBlockScore> comparator=Comparator.comparing(TriadsBlock::getTriadsID);
            Arrays.sort(triadsBlockScoreArray, comparator);
            String triadsBlockID;
            double[] xpclrScoreABD;
            for (TriadsBlockScore triadsBlockScore:triadsBlockScoreArray){
                triadsBlockID=triadsBlockScore.getTriadsID();
                xpclrScoreABD=triadsBlockScore.calculateXPCLRScore();
                for (Sub sub : Sub.values()){
                    sb.setLength(0);
                    sb.append(triadsBlockID).append("\t");
                    sb.append(sub).append("\t");
                    sb.append(xpclrScoreABD[sub.getIndex()]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
