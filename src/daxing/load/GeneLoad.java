package daxing.load;

import gnu.trove.list.array.TByteArrayList;

import java.util.HashMap;
import java.util.Map;

public class GeneLoad{

    //0: 0/0 1: 1/1
    String geneName;
    TByteArrayList synGenotype;
    TByteArrayList nonsynGenotype;
    TByteArrayList hgDeleteriousGenop;

    public static Map<String, Byte> genotypeToByteMap1 =initializeGenotypeToByteMap1();
    private static Map<String, Byte> initializeGenotypeToByteMap1(){
        String[] genotype={"0/0", "1/1"};
        Map<String,Byte> map=new HashMap<>();
        map.put(genotype[0], (byte) 0);
        map.put(genotype[1], (byte) 1);
        return map;
    }

    public static Map<String, Byte> genotypeToByteMap2 =initializeGenotypeToByteMap2();
    private static Map<String, Byte> initializeGenotypeToByteMap2(){
        String[] genotype={"0/0", "1/1"};
        Map<String,Byte> map=new HashMap<>();
        map.put(genotype[0], (byte) 1);
        map.put(genotype[1], (byte) 0);
        return map;
    }

    public GeneLoad(String geneName){
        this.geneName=geneName;
        this.synGenotype =new TByteArrayList();
        this.nonsynGenotype=new TByteArrayList();
        this.hgDeleteriousGenop=new TByteArrayList();
    }

    public String getGeneName() {
        return geneName;
    }

    public void addGenotype(byte[] indexGenotype){
        if (indexGenotype[0]==0){
            this.synGenotype.add(indexGenotype[1]);
        }else if(indexGenotype[0]==1){
            this.nonsynGenotype.add(indexGenotype[1]);
        }else if (indexGenotype[0]==2){
            this.hgDeleteriousGenop.add(indexGenotype[1]);
        }
    }

    public int getSynNum(){
        return this.synGenotype.size();
    }

    public int getNonsynNum(){
        return this.nonsynGenotype.size();
    }

    public int getHGDeleteriousNum(){
        return this.hgDeleteriousGenop.size();
    }

    public int getSynDerivedNum(){
        int sum=0;
        for (int i = 0; i < this.synGenotype.size(); i++) {
            sum+=this.synGenotype.get(i);
        }
        return sum;
    }

    public int getNonsynDerivedNum(){
        int sum=0;
        for (int i = 0; i < this.nonsynGenotype.size(); i++) {
            sum+=this.nonsynGenotype.get(i);
        }
        return sum;
    }

    public int getHGDeleteriousDerivedNum(){
        int sum=0;
        for (int i = 0; i < this.hgDeleteriousGenop.size(); i++) {
            sum+=this.hgDeleteriousGenop.get(i);
        }
        return sum;
    }

    public static byte caculateGenotype(String genotype, boolean isRefAlleleAncestral){
        if (isRefAlleleAncestral){
            return genotypeToByteMap1.get(genotype);
        }
        return genotypeToByteMap2.get(genotype);
    }
}
