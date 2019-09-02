package daxing.common;

import format.position.ChrPos;
import gnu.trove.list.array.TIntArrayList;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * VCF coordinates are 1-based;
 * 0-based coordinates in ChrNameMap class
 */
public class ChrConvertionRule {
    private  int[] chrID;
    private  int[] endIndex;
    private  String[] oriChrName;
    private  int[] startIndexOnOriChr;
    private  int[] endIndexOnOriChr;

    public ChrConvertionRule(Path chrConvertionRuleFile){
        this.initialize(chrConvertionRuleFile);
    }

    private void initialize(Path inputFile){
        TIntArrayList chrIDlist=new TIntArrayList();
        TIntArrayList endIndexList=new TIntArrayList();
        List<String> oriChrNameList=new ArrayList<>();
        TIntArrayList startIndexOnOriChrList=new TIntArrayList();
        TIntArrayList endIndexOnOriChiList=new TIntArrayList();
        try(BufferedReader br= Files.newBufferedReader(inputFile)){
            br.readLine();
            String line;
            List<String> lineList;
            while ((line=br.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                chrIDlist.add(Byte.parseByte(lineList.get(0)));
                endIndexList.add(Integer.parseInt(lineList.get(2)));
                oriChrNameList.add(lineList.get(3).replace("chr",""));
                startIndexOnOriChrList.add(Integer.parseInt(lineList.get(4)));
                endIndexOnOriChiList.add(Integer.parseInt(lineList.get(5)));
            }
            this.chrID=chrIDlist.toArray();
            this.endIndex=endIndexList.toArray();
            this.oriChrName=oriChrNameList.stream().toArray(String[]::new);
            this.startIndexOnOriChr=startIndexOnOriChrList.toArray();
            this.endIndexOnOriChr=endIndexOnOriChiList.toArray();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public int[] getChrID() {
        return chrID;
    }

    public int[] getEndIndex() {
        return endIndex;
    }

    public int[] getEndIndexOnOriChr() {
        return endIndexOnOriChr;
    }

    public int[] getStartIndexOnOriChr() {
        return startIndexOnOriChr;
    }

    public String[] getOriChrName() {
        return oriChrName;
    }

    /**
     *
     * @param chrID 1,  3,  5,  0,  43,  44   et al.
     * @return      1A, 1B, 1D, Un, Mit, Chl  et al.
     */
    public String getOriChrNameFromChrID(int chrID){
        int index=Arrays.binarySearch(this.getChrID(), chrID);
        if(index < 0){
            System.out.println("chrID must between 0~44");
            System.exit(1);
        }
        return this.getOriChrName()[index];
    }

    /**
     *
     * @param chr 1A, 1B, 1D, Un, Mit, Chl  et al.
     * @return    [1,2]  [3,4]  [5,6]  [0,-1]  [43,-1]  [44,-1]   et al.
     */
    public int[] getChrIDFromOriChrName(String chr){
        int[] res={-1, -1};
        Set<String> s=new HashSet<>(CollectionTool.changeToList(this.oriChrName));
        String[] oriChrNameArray=s.toArray(new String[s.size()]);
        Arrays.sort(oriChrNameArray);
        int index=Arrays.binarySearch(oriChrNameArray, chr);
        if(index < 0){
            System.out.println(chr+" is not a chromosome of wheat");
            System.exit(1);
        }else if (index<=20){
            index=2*index+1;
            res[0]=index;
            res[1]=index+1;
            return res;
        }else if (index==21){
            res[0]=44;
            return res;
        }else if (index==22){
            res[0]=43;
            return res;
        }
        res[0]=0;
        return res;

    }

    /**
     * 0-based coordinates in ChrNameMap class
     * @param chr
     * @param positionInOriChrName coordinates are 1-based
     * @return ChrID
     */
    public int getChrIDFromOriChrName(String chr, int positionInOriChrName){
        int coordinateBased_0=positionInOriChrName-1;
        int[] indexArray=this.getChrIDFromOriChrName(chr);
        if (coordinateBased_0<this.getEndIndexOnOriChr()[indexArray[0]]){
            return indexArray[0];
        }if (indexArray[0]!=-1 && coordinateBased_0>=this.getEndIndexOnOriChr()[indexArray[1]]){
            System.out.println(positionInOriChrName +" is larger than "+chr+" size");
            System.exit(1);
        }
        return indexArray[1];
    }

    /**
     *
     * @param chr
     * @param positionInOriChrName coordinates are 1-based
     * @return position which is 0-based coordinates in ChrNameMap class
     */
    public int getPositionFromOriChrName(String chr, int positionInOriChrName){
        int coordinateBased_0=positionInOriChrName-1;
        int[] indexArray=this.getChrIDFromOriChrName(chr);
        if (coordinateBased_0<this.getEndIndexOnOriChr()[indexArray[0]]){
            return coordinateBased_0;
        }else if (indexArray[1]!=-1 && coordinateBased_0>=this.getEndIndexOnOriChr()[indexArray[1]]){
            System.out.println(positionInOriChrName +" is larger than "+chr+" size");
            System.exit(1);
        }
        return coordinateBased_0-this.getEndIndexOnOriChr()[indexArray[0]];
    }

    /**
     *
     * @param chrID    1,  3,  5,  0,  43,  44   et al.
     * @param position coordinates are 0-based in ChrNameMap class
     * @return position which is 1-based coordinates in OriChr
     */
    public int getOriChrPositin(int chrID, int position){
        if (position>=this.getEndIndex()[chrID]){
            System.out.println(position+" is larger than "+chrID+" chromosome which has a size "+((this.getEndIndex()[chrID])-1));
            System.exit(1);
        }
        int[] even2_42= IntStream.iterate(2, n->n+2).limit(21).toArray();
        int res=Arrays.binarySearch(even2_42, chrID);
        if (res<0){
            return position+1;
        }
        return this.getStartIndexOnOriChr()[chrID]+position+1;
    }

    /**
     * 将vcf的ChrPos(1-based)转换为ref坐标（1-based）
     * @param vcfChrPos
     * @return reference position of VCF ChrPos
     */
    public int getRefPositionFromVCF(ChrPos vcfChrPos){
        int chr=vcfChrPos.getChromosome();
        int pos_O_based=vcfChrPos.getPosition()-1;
        return this.getOriChrPositin(chr, pos_O_based);
    }

    /**
     *
     * @param chr 1A, Un, et al.
     * @param pos coordinates are 1-based
     * @return coordinates in vcf
     */
    public int getVCFPositionFromOriChrName(String chr, int pos){
        int chrID=this.getChrIDFromOriChrName(chr, pos);
        return pos-this.getStartIndexOnOriChr()[chrID];
    }


    public static Map<Integer, String> getChrID_OriChrMap(){
        Map<Integer,String> ChrID_OriChrMap=new HashMap<>();
        List<Integer> numOfChr= IntStream.range(1,43).boxed().collect(Collectors.toList());
        List<Integer> int1_7= IntStream.range(1,8).boxed().collect(Collectors.toList());
        List<Integer> chrList=new ArrayList<>();
        for(int i=0;i<6;i++){
            chrList.addAll(int1_7);
        }
        Collections.sort(chrList);
        String abd=String.join("", Collections.nCopies(7,"AABBDD"));
        ChrID_OriChrMap.put(0, "Un");
        for(int i=0;i<numOfChr.size();i++){
            ChrID_OriChrMap.put(numOfChr.get(i),String.valueOf(chrList.get(i))+abd.charAt(i));
        }
        ChrID_OriChrMap.put(43, "Mit");
        ChrID_OriChrMap.put(44, "Chl");
        return ChrID_OriChrMap;
    }
}
