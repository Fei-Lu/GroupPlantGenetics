package daxing.load.neutralSiteLoad;

import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;

public class DynamicSNPGenotypeDB {

    List<SNPGenotype> snpList;

    public DynamicSNPGenotypeDB(){
        snpList =new ArrayList<>();
    }

    public void addSNPGenotype(SNPGenotype snpGenotype){
        this.snpList.add(snpGenotype);
    }

    public int size(){
        return snpList.size();
    }

    public int getTaxonNum(){
        return this.snpList.get(0).genotypeList.size();
    }

    public int getGenotype(int snpIndex, int taxonIndex){
        return snpList.get(snpIndex).genotypeList.get(taxonIndex);
    }

    public TIntArrayList[] getAllTaxonGenotype(){
        TIntArrayList[] taoxnGenotype=new TIntArrayList[this.getTaxonNum()];
        for (int i = 0; i < taoxnGenotype.length; i++) {
            taoxnGenotype[i]=new TIntArrayList();
        }
        for (int i = 0; i < snpList.size(); i++) {
            for (int j = 0; j < this.getTaxonNum(); j++) {
                taoxnGenotype[j].add(this.getGenotype(i, j));
            }
        }
        return taoxnGenotype;
    }

    public int[] countAllTaxonDerived(){
        TIntArrayList[] taoxnGenotype=this.getAllTaxonGenotype();
        int[] taxonDerivedCount=new int[taoxnGenotype.length];
        for (int i = 0; i < taoxnGenotype.length; i++) {
            for (int j = 0; j < taoxnGenotype[i].size(); j++) {
                if (taoxnGenotype[i].get(j)!=1) continue;
                taxonDerivedCount[i]=taxonDerivedCount[i]+1;
            }
        }
        return taxonDerivedCount;
    }

    public boolean retainAll(GenesDB.GeneRange geneRange){
        List<SNPGenotype> snpGenotypeList=this.snpList;
        Predicate<SNPGenotype> p=snp->snp.getPosition() < geneRange.start;
        boolean res=snpGenotypeList.removeIf(p);
        this.snpList=snpGenotypeList;
        return res;
    }
}
