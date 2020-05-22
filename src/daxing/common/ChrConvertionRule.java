package daxing.common;

import pgl.infra.pos.ChrPos;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.IntStream;

/**
 * VCF coordinates are 1-based;
 * 0-based coordinates in ChrNameMap class
 * @author Daxing Xu
 */
public class ChrConvertionRule {
    private  int[] chrID;
    private  int[] endIndex;
    private  String[] oriChrName;
    private  int[] startIndexOnOriChr;
    private  int[] endIndexOnOriChr;

    public ChrConvertionRule(Path chrConvertionRuleFile){
        this.initialize(chrConvertionRuleFile.toString());
    }

    public ChrConvertionRule(String chrConvertionRuleFile){
        this.initialize(chrConvertionRuleFile);
    }

    private void initialize(String inputFile){
        TIntArrayList chrIDlist=new TIntArrayList();
        TIntArrayList endIndexList=new TIntArrayList();
        List<String> oriChrNameList=new ArrayList<>();
        TIntArrayList startIndexOnOriChrList=new TIntArrayList();
        TIntArrayList endIndexOnOriChiList=new TIntArrayList();
        try(BufferedReader br= IOUtils.getTextReader(inputFile)){
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
    public String getRefChrFromVCFChr(int chrID){
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
    public int[] getVCFChrFromRefChr(String chr){
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
     *
     * @param chr 1A, 1B, 1D, Un, Mit, Chl  et al.
     * @param positionInRef_1_based coordinates are ref_1_based
     * @return chrID in vcf
     */
    public int getVCFChrFromRefChrPos(String chr, int positionInRef_1_based){
        int[] indexArray=this.getVCFChrFromRefChr(chr);
        int coordinateBased_0=positionInRef_1_based-1;
        int endIndexOnOriChr=this.getEndIndexOnOriChr()[indexArray[0]];
        if (coordinateBased_0 < endIndexOnOriChr){
            return indexArray[0];
        }else if (indexArray[0]==0 || indexArray[0]==43 || indexArray[0]==44){
            System.out.println(positionInRef_1_based +" is larger than "+chr+" size: "+endIndexOnOriChr);
            System.exit(1);
        }else if (coordinateBased_0>=this.getEndIndexOnOriChr()[indexArray[1]]){
            System.out.println(positionInRef_1_based +" is larger than "+chr+" size: "+this.getEndIndexOnOriChr()[indexArray[1]]);
            System.exit(1);
        }
        return indexArray[1];
    }

    /**
     *
     * @param chr 1A, 1B, 1D, Un, Mit, Chl  et al.
     * @param positionInOriChrName coordinates are ref_1_based
     * @return position in vcf
     */
    public int getVCFPosFromRefChrPos(String chr, int positionInOriChrName){
        int coordinateBased_0=positionInOriChrName-1;
        int[] indexArray=this.getVCFChrFromRefChr(chr);
        int endIndexOnOriChr=this.getEndIndexOnOriChr()[indexArray[0]];
        if (coordinateBased_0 < endIndexOnOriChr){
            return positionInOriChrName;
        }else if (indexArray[0]==0 || indexArray[0]==43 || indexArray[0]==44){
            System.out.println(positionInOriChrName +" is larger than "+chr+" size: "+endIndexOnOriChr);
            System.exit(1);
        }else if (coordinateBased_0>=this.getEndIndexOnOriChr()[indexArray[1]]){
            System.out.println(positionInOriChrName +" is larger than "+chr+" size: "+this.getEndIndexOnOriChr()[indexArray[1]]);
            System.exit(1);
        }
        return positionInOriChrName-endIndexOnOriChr;
    }

    /**
     *
     * @param chr 1A, 1B, 1D, Un, Mit, Chl  et al.
     * @param positionInOriChrName coordinates are ref_1_based
     * @return chrPos in vcf
     */
    public ChrPos getVCFChrPosFromRefChrPos(String chr, int positionInOriChrName){
        short chrID=(short) this.getVCFChrFromRefChrPos(chr, positionInOriChrName);
        int pos=this.getVCFPosFromRefChrPos(chr, positionInOriChrName);
        return new ChrPos(chrID, pos);
    }

    /**
     *
     * @param chrID   0,  1,  3,  5,  43,  44   et al.
     * @param position position in vcf
     * @return coordinates are ref_1_based
     */
    public int getRefPosFromVCFChrPos(int chrID, int position){
        if (position>this.getEndIndex()[chrID]){
            System.out.println(position+" is larger than "+chrID+" chromosome which has a size "+((this.getEndIndex()[chrID])));
            System.exit(1);
        }
        int[] even2_42= IntStream.iterate(2, n->n+2).limit(21).toArray();
        int res=Arrays.binarySearch(even2_42, chrID);
        if (res<0){
            return position;
        }
        return this.getStartIndexOnOriChr()[chrID]+position;
    }

    /**
     *
     * @param chrPos in vcf
     * @return
     */
    public int getRefPosFromVCFChrPos(ChrPos chrPos){
        return this.getRefPosFromVCFChrPos(chrPos.getChromosome(), chrPos.getPosition());
    }

    /**
     *
     * @param chr 1A, 1B, 1D, Un, Mit, Chl  et al.
     * @return size of chr
     */
    public int getChrSize(String chr){
        int chrSize=Integer.MIN_VALUE;
        int[] chrID=this.getVCFChrFromRefChr(chr);
        if (chrID[1]>chrID[0]){
            chrSize=this.getEndIndexOnOriChr()[chrID[1]];
        }else {
            chrSize=this.getEndIndexOnOriChr()[chrID[0]];
        }
        return chrSize;
    }
}
