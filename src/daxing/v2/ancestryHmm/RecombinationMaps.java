package daxing.v2.ancestryHmm;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import daxing.common.factors.WheatLineage;
import daxing.common.utiles.IOTool;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

public class RecombinationMaps{

    String[] refChrs;
    int[] chrRefCodes;
    int[] positions;
    double[] geneticPositions;

    public static BiMap<String, Integer> getRefChr2RecodeBiMap(){
        BiMap<String, Integer> chr2RecodeMap = HashBiMap.create();
        List<String> refChrList = WheatLineage.abdLineage();
        for (int i = 0; i < refChrList.size(); i++) {
            chr2RecodeMap.put(refChrList.get(i), i+1);
        }
        return chr2RecodeMap;
    }

    public RecombinationMaps(String geneticsMapFile){
        List<String> refChrList= new ArrayList<>();
        TIntArrayList chrRefCodes = new TIntArrayList();
        TIntArrayList positionList = new TIntArrayList();
        TDoubleArrayList geneticPositionList = new TDoubleArrayList();
        BiMap<String, Integer> chr2RecodeMap=RecombinationMaps.getRefChr2RecodeBiMap();
        try (BufferedReader br = IOTool.getReader(geneticsMapFile)) {
            String line, chr;
            List<String> temp;
            int position;
            double geneticPosition;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                chr = temp.get(1).substring(3);
                position = Integer.parseInt(temp.get(2));
                geneticPosition = Double.parseDouble(temp.get(3));
                refChrList.add(chr);
                chrRefCodes.add(chr2RecodeMap.get(chr));
                positionList.add(position);
                geneticPositionList.add(geneticPosition);
            }
            this.refChrs = refChrList.toArray(new String[refChrList.size()]);
            this.chrRefCodes = chrRefCodes.toArray();
            this.positions = positionList.toArray();
            this.geneticPositions = geneticPositionList.toArray();
        } catch (IOException e) {
            e.printStackTrace();
        }
        this.sort();
    }

    public int getMarkerNum(){
        return this.refChrs.length;
    }

    public int getRefChrCode(String chr){
        return RecombinationMaps.getRefChr2RecodeBiMap().get(chr);
    }

    Swapper swapper = (i, i1) -> {
        String c1;
        int p1;
        double gp1;

        c1= refChrs[i];
        refChrs[i]= refChrs[i1];
        refChrs[i1]=c1;

        p1 = positions[i];
        positions[i] = positions[i1];
        positions[i1] = p1;

        gp1 = geneticPositions[i];
        geneticPositions[i] = geneticPositions[i1];
        geneticPositions[i1] = gp1;
    };

    IntComparator intComparator = (a, b) -> {
        if (refChrs[a].equals(refChrs[b])){
            if (positions[a]==positions[b]){
                return 0;
            }else if (positions[a] < positions[b]){
                return -1;
            }
            return 1;
        }else {
            return refChrs[a].compareTo(refChrs[b]);
        }
    };

    private void sort(){
        GenericSorting.quickSort(0, refChrs.length, intComparator, swapper);
    }

    private int findFirstIndex(String refChr){
        int targetCode = this.getRefChrCode(refChr);
        int low = 0;
        int high = this.getMarkerNum()-1;
        while (low <= high){
            int mid = low + (high - low)/2;
            if ((mid==0 || this.chrRefCodes[mid-1] < targetCode) && this.chrRefCodes[mid]==targetCode){
                return mid;
            }else if (targetCode > this.chrRefCodes[mid]){
                low=mid + 1;
            }else {
                high = mid -1;
            }
        }
        return -1;
    }

    private int findLastIndex(String refChr){
        int targetCode = this.getRefChrCode(refChr);
        int low = 0;
        int high = this.getMarkerNum()-1;
        while (low <= high){
            int mid = low + (high - low)/2;
            if ((mid== this.getMarkerNum()-1 || this.chrRefCodes[mid+1] > targetCode) && this.chrRefCodes[mid]==targetCode){
                return mid;
            }else if (targetCode < this.chrRefCodes[mid]){
                high = mid -1;
            }else {
                low=mid + 1;
            }
        }
        return -1;
    }

    public int[] findFirstLastIndex(String refChr){
        int[] firstLast = new int[2];
       Arrays.fill(firstLast, -1);
       firstLast[0] = this.findFirstIndex(refChr);
       firstLast[1] = this.findLastIndex(refChr);
       return firstLast;
    }

    /**
     *
     * @param refChr
     * @param position
     * @return negative value if position not in database
     */
    public int findPositionIndex(String refChr, int position){
        int[] firstLastIndex = this.findFirstLastIndex(refChr);
        return Arrays.binarySearch(this.positions, firstLastIndex[0], firstLastIndex[1], position);
    }

    public double getGeneticPositions(String refChr, int position){
        int hit = this.findPositionIndex(refChr, position);
        int positionIndex1, positionIndex2;
        int positionIndex=-1;
        if (hit < 0){
            if (hit == -1){
                positionIndex = -hit-1;
            }else {
                positionIndex1 = -hit-2;
                positionIndex2 =-hit-1;
                if ((position-this.positions[positionIndex1]) < (this.positions[positionIndex2]-position)){
                    positionIndex = positionIndex1;
                }else {
                    positionIndex = positionIndex2;
                }
            }
        }
        return this.geneticPositions[positionIndex];
    }

    public double getGeneticPositions(int chrID, int position){
        String refChr = RefV1Utils.getChromosome(chrID, position);
        int refPos = RefV1Utils.getPosOnChromosome(chrID, position);
        int hit = this.findPositionIndex(refChr, refPos);
        int positionIndex1, positionIndex2;
        int positionIndex=-1;
        if (hit < 0){
            if (hit == -1){
                positionIndex = -hit-1;
            }else {
                positionIndex1 = -hit-2;
                positionIndex2 =-hit-1;
                if ((refPos-this.positions[positionIndex1]) < (this.positions[positionIndex2]-refPos)){
                    positionIndex = positionIndex1;
                }else {
                    positionIndex = positionIndex2;
                }
            }
        }
        return this.geneticPositions[positionIndex];
    }
}
