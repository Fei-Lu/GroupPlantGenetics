package daxing.temp;

import com.google.common.collect.Table;
import daxing.common.ChrRange;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * this class used to record population introgression region
 */
public class IntrogressionRegion {

    ChrRange[] chrRanges;
    double[][][] fdModifyBins;//dim1:region, dim2: P2, dim3: P3

    public enum P2{
        CL(0),LR_AF(1),LR_AM(2),LR_CSA(3),LR_EA(4),LR_EU(5),LR_WA(6);

        int index;
        P2(int index){
            this.index=index;
        }

        public int getIndex() {
            return index;
        }
    }

    public enum P3{
        WE(0),DE(1),FTT(2),AT(3);

        int index;
        P3(int index){
            this.index=index;
        }

        static P3 newInstanceFrom(int index){
            switch (index){
                case 0:
                    return WE;
                case 1:
                    return DE;
                case 2:
                    return FTT;
                case 3:
                    return AT;
            }
            return null;
        }

        public int getIndex() {
            return index;
        }
    }

    public IntrogressionRegion(String populationIntrogressionFile){
        ChrRange[] chrRanges=initialize(populationIntrogressionFile);
        Arrays.sort(chrRanges);
        this.chrRanges=chrRanges;
        fdModifyBins = new double[chrRanges.length][][];
        for (int i = 0; i < fdModifyBins.length; i++) {
            fdModifyBins[i] = new double[P2.values().length][];
            for (int j = 0; j < fdModifyBins[i].length; j++) {
                fdModifyBins[i][j]= new double[P3.values().length];
                Arrays.fill(fdModifyBins[i][j], -1);
            }
        }
        try (BufferedReader br = IOTool.getReader(populationIntrogressionFile)) {
            String line;
            List<String> temp;
            br.readLine();
            String chr;
            int start, end, index;
            ChrRange chrRange;
            P2 p2;
            P3 p3;
            double bin;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chr = temp.get(0);
                start = Integer.parseInt(temp.get(1));
                end = Integer.parseInt(temp.get(2));
                chrRange = new ChrRange(chr, start, end);
                index = Arrays.binarySearch(this.getChrRanges(), chrRange);
                assert index >=0 : "check index";
                p2=P2.valueOf(temp.get(4));
                p3=P3.valueOf(temp.get(5));
                bin = Double.parseDouble(temp.get(6));
                fdModifyBins[index][p2.index][p3.index]=bin;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private ChrRange[] initialize(String populationIntrogressionFile){
        Table<String, String, String> table=RowTableTool.getTable(populationIntrogressionFile, 0,1,2);
        List<ChrRange> chrRangeList= new ArrayList<>(table.size());
        String chr;
        int start, end;
        ChrRange chrRange;
        for (Table.Cell<String,String,String> cell: table.cellSet()){
            chr=cell.getRowKey();
            start=Integer.parseInt(cell.getColumnKey());
            end=Integer.parseInt(cell.getValue());
            chrRange = new ChrRange(chr, start, end);
            chrRangeList.add(chrRange);
        }
        return chrRangeList.toArray(new ChrRange[0]);
    }

    public ChrRange[] getChrRanges() {
        return chrRanges;
    }

    public double[][][] getFdModifyBins() {
        return fdModifyBins;
    }

    public int getChrRangeIndexFrom(String refChr, int refPos){
        ChrRange chrRange = new ChrRange(refChr, refPos, refPos+1);
        int hit = Arrays.binarySearch(this.getChrRanges(), chrRange);
        int index = hit;
        if (index < -1) {
            index = -index-2;
            if (this.isWithinThisChrRange(index, refChr, refPos)) return index;
        }
        return hit;
    }

    public double[] getP3fdModifyBins(int chrID, int pos, P2 p2){
        String refChr= RefV1Utils.getChromosome(chrID, pos);
        int refPos=RefV1Utils.getPosOnChromosome(chrID, pos);
        int index = getChrRangeIndexFrom(refChr, refPos);
        assert index >=0:"check chrID pos";
        return this.fdModifyBins[index][p2.index];
    }

    private boolean isWithinThisChrRange(int index, String refChr, int refPos){
        if (refChr!=this.chrRanges[index].getChr()) return false;
        if (refPos < this.chrRanges[index].getStart()) return false;
        if (refPos >= this.chrRanges[index].getEnd()) return false;
        return  true;
    }
}
