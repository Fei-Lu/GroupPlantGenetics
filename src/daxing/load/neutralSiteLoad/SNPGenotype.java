package daxing.load.neutralSiteLoad;

import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.utils.PStringUtils;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class SNPGenotype extends BiSNP {

    TIntArrayList genotypeList; //0: 0/0 1: 1/1 2: 0/1 3 ./.

    private static Map<String, Byte> genotypeToByteMap1 =initializeGenotypeToByteMap1();
    private static Map<String, Byte> initializeGenotypeToByteMap1(){
        String[] genotype={"0/0", "1/1", "0/1", "./."};
        Map<String,Byte> map=new HashMap<>();
        map.put(genotype[0], (byte) 0);
        map.put(genotype[1], (byte) 1);
        map.put(genotype[2], (byte) 2);
        map.put(genotype[3], (byte) 3);
        return map;
    }

    private static Map<String, Byte> genotypeToByteMap2 =initializeGenotypeToByteMap2();
    private static Map<String, Byte> initializeGenotypeToByteMap2(){
        String[] genotype={"0/0", "1/1", "0/1", "./."};
        Map<String,Byte> map=new HashMap<>();
        map.put(genotype[0], (byte) 1);
        map.put(genotype[1], (byte) 0);
        map.put(genotype[2], (byte) 2);
        map.put(genotype[3], (byte) 3);
        return map;
    }

    private SNPGenotype(short chr, int pos, char refBase, char altBase, String info, char ancestralBase,
                        TIntArrayList genotypeList){
        super(chr, pos, refBase, altBase, info);
        this.genotypeList=genotypeList;
        setAncestral(ancestralBase);
    }

    public static SNPGenotype getSNPGenotype(String vcfLine, char ancestralBase, String variant_type, List<String> genotypeList){
        List<String> temp= PStringUtils.fastSplit(vcfLine);
        short chr=Short.parseShort(temp.get(0));
        int pos=Integer.parseInt(temp.get(1));
        char refBase=temp.get(3).charAt(0);
        char altBase=temp.get(4).charAt(0);
        TIntArrayList res=new TIntArrayList();
        List<String> genotype=genotypeList.stream().map(s->s.substring(0,3)).collect(Collectors.toList());
        for (int i = 0; i < genotype.size(); i++) {
            if (ancestralBase==refBase){
                res.add(genotypeToByteMap1.get(genotype.get(i)));
            }else if (ancestralBase==altBase){
                res.add(genotypeToByteMap2.get(genotype.get(i)));
            }
        }
        return new SNPGenotype(chr, pos, refBase, altBase, variant_type, ancestralBase, res);

    }

    private void setAncestral(char ancestralBase){
        char refAlleleBase=this.reference.getAlleleBase();
        char altAlleleBase=this.alternative.getAlleleBase();
        if (refAlleleBase==ancestralBase){
            this.setReferenceAlleleType(AlleleType.Ancestral);
        } else if (altAlleleBase==ancestralBase){
            this.setAlternativeAlleleType(AlleleType.Ancestral);
        }
    }


}
