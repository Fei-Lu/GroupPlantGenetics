package daxing.individualIntrogression;

import daxing.common.IOTool;
import daxing.common.LoadType;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Only record derived
 */
public class IndividualChrPosLoad{

    String taxonName;
    int chr;
    TIntArrayList posList;
    List<LoadType> synNonDelList;
    TByteArrayList ifHeter;
    TByteArrayList ifHomozygousDerived;

    public static Map<String, Byte> genotypeToByteMap1 =initializeGenotypeToByteMap1();
    private static Map<String, Byte> initializeGenotypeToByteMap1(){
        String[] genotype={"0/0", "1/1", "0/1"};
        Map<String,Byte> map=new HashMap<>();
        map.put(genotype[0], (byte) 0);
        map.put(genotype[1], (byte) 1);
        map.put(genotype[2], (byte) 2);
        return map;
    }

    public static Map<String, Byte> genotypeToByteMap2 =initializeGenotypeToByteMap2();
    private static Map<String, Byte> initializeGenotypeToByteMap2(){
        String[] genotype={"0/0", "1/1", "0/1"};
        Map<String,Byte> map=new HashMap<>();
        map.put(genotype[0], (byte) 1);
        map.put(genotype[1], (byte) 0);
        map.put(genotype[2], (byte) 2);
        return map;
    }

    public IndividualChrPosLoad(String taxonName, int chr){
        this.taxonName=taxonName;
        this.chr=chr;
        this.posList=new TIntArrayList();
        this.synNonDelList=new ArrayList<>();
        this.ifHeter=new TByteArrayList();
        this.ifHomozygousDerived=new TByteArrayList();
    }

    public int getChr() {
        return chr;
    }

    public List<LoadType> getSynNonDelList() {
        return synNonDelList;
    }

    public String getTaxonName() {
        return taxonName;
    }

    public TByteArrayList getIfHeter() {
        return ifHeter;
    }

    public TByteArrayList getIfHomozygousDerived() {
        return ifHomozygousDerived;
    }

    public TIntArrayList getPosList() {
        return posList;
    }

    public void addChrPos(int chr, int pos, LoadType synNonDel, boolean ifHeter, boolean ifHomozygousDerived){
        if (chr!=this.getChr()){
            System.out.println("please check your parameter chr: "+chr);
            System.exit(1);
        }
        this.posList.add(pos);
        this.synNonDelList.add(synNonDel);
        byte ifHeterByte= (byte) (ifHeter ? 1 : 0);
        byte ifHomozygousDerivedByte= (byte) (ifHomozygousDerived ? 1 : 0);
        this.ifHeter.add(ifHeterByte);
        this.ifHomozygousDerived.add(ifHomozygousDerivedByte);
    }

    public void write(String outDir){
        int chr=this.getChr();
        String taxonName=this.getTaxonName();
        TIntArrayList posList=this.getPosList();
        List<LoadType> loadTypes=this.getSynNonDelList();
        TByteArrayList ifHeter=this.getIfHeter();
        TByteArrayList ifHomozygousDerived=this.getIfHomozygousDerived();
        File outFile=new File(outDir, "chr"+ PStringUtils.getNDigitNumber(3, chr)+"."+taxonName+".txt.gz");
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Chr\tPos\tLoadType\tIfHeter\tIfHomozygousDerived");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < posList.size(); i++) {
                sb.setLength(0);
                sb.append(chr).append("\t").append(posList.get(i)).append("\t").append(loadTypes.get(i));
                sb.append("\t").append(ifHeter.get(i)).append("\t").append(ifHomozygousDerived.get(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 0: 0/0 1: 1/1 2: 0/1 0为ancestral, 1为derived, 2为杂合
     * @param genotype
     * @param isRefAlleleAncestral
     * @return
     */
    public static byte caculateGenotype(String genotype, boolean isRefAlleleAncestral){
        if (isRefAlleleAncestral){
            return genotypeToByteMap1.get(genotype);
        }
        return genotypeToByteMap2.get(genotype);
    }

}
