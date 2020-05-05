package daxing.load.neutralSiteLoad;

import daxing.load.slidingWindow.GeneWindowDB;
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

    public TIntArrayList[] getAllTaxonNonsynGenotype(){
        TIntArrayList[] taoxnNonsynGenotype=new TIntArrayList[this.getTaxonNum()];
        for (int i = 0; i < taoxnNonsynGenotype.length; i++) {
            taoxnNonsynGenotype[i]=new TIntArrayList();
        }
        for (int i = 0; i < snpList.size(); i++) {
            if (!snpList.get(i).getSNPInfo().equals("NONSYNONYMOUS")) continue;
            for (int j = 0; j < this.getTaxonNum(); j++) {
                taoxnNonsynGenotype[j].add(this.getGenotype(i, j));
            }
        }
        return taoxnNonsynGenotype;
    }

    public int[][] countAllTaxonDerived(){
        TIntArrayList[] taoxnGenotype=this.getAllTaxonGenotype();
        int[][] taxonDerivedCount=new int[taoxnGenotype.length][];
        for (int i = 0; i < taxonDerivedCount.length; i++) {
            taxonDerivedCount[i]=new int[2];
        }
        for (int i = 0; i < taoxnGenotype.length; i++) {
            for (int j = 0; j < taoxnGenotype[i].size(); j++) {
                if (taoxnGenotype[i].get(j)==2 || taoxnGenotype[i].get(j) == 3) continue; //0: 0/0 1: 1/1 2: 0/1 3 ./.
                taxonDerivedCount[i][0]=taxonDerivedCount[i][0]+1;
                if (taoxnGenotype[i].get(j)==1){
                    taxonDerivedCount[i][1]=taxonDerivedCount[i][1]+1;
                }
            }
        }
        return taxonDerivedCount;
    }

    public int[][] countAllTaxonNonsynDerived(){
        TIntArrayList[] taoxnNonsynGenotype=this.getAllTaxonNonsynGenotype();
        int[][] taxonNonsynDerivedCount=new int[taoxnNonsynGenotype.length][];
        for (int i = 0; i < taxonNonsynDerivedCount.length; i++) {
            taxonNonsynDerivedCount[i]=new int[2];
        }
        for (int i = 0; i < taoxnNonsynGenotype.length; i++) {
            for (int j = 0; j < taoxnNonsynGenotype[i].size(); j++) {
                if (taoxnNonsynGenotype[i].get(j)==2 || taoxnNonsynGenotype[i].get(j) == 3) continue; //0: 0/0 1: 1/1 2: 0/1 3 ./.
                taxonNonsynDerivedCount[i][0]=taxonNonsynDerivedCount[i][0]+1;
                if (taoxnNonsynGenotype[i].get(j)==1){
                    taxonNonsynDerivedCount[i][1]=taxonNonsynDerivedCount[i][1]+1;
                }
            }
        }
        return taxonNonsynDerivedCount;
    }

    public boolean retainAll(GenesDB.GeneRange geneRange){
        List<SNPGenotype> snpGenotypeList=this.snpList;
        Predicate<SNPGenotype> p=snp->snp.getPosition() < geneRange.start;
        boolean res=snpGenotypeList.removeIf(p);
        this.snpList=snpGenotypeList;
        return res;
    }

    public boolean retainAll(GeneWindowDB.GeneWindow geneWindow){
        List<SNPGenotype> snpGenotypeList=this.snpList;
        Predicate<SNPGenotype> inGeneWindow=snp->geneWindow.containCDSPos(snp.getChromosome(), snp.getPosition());
        boolean res=snpGenotypeList.removeIf(inGeneWindow.negate());
        this.snpList=snpGenotypeList;
        return res;
    }
}
