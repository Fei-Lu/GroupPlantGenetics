package daxing.v2.localAncestryInfer.laidp;

import daxing.common.utiles.IOTool;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

public class TaxaGroup {

    List<String> taxaList;
    List<String> popList;

    List<Source> sourceList;

    public TaxaGroup(List<String> taxaList, List<String> popList, List<Source> sourceList){
        this.taxaList=taxaList;
        this.popList=popList;
        this.sourceList=sourceList;
    }

    public static TaxaGroup buildFrom(String taxaGroupFile){
        List<String> taxaList = new ArrayList<>();
        List<String> popList = new ArrayList<>();
        List<Source> sourceList = new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(taxaGroupFile)) {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                if (temp.get(2).equals("admixed")) continue;
                taxaList.add(temp.get(0));
                popList.add(temp.get(1));
                sourceList.add(Source.valueOf(temp.get(2)));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return new TaxaGroup(taxaList, popList, sourceList);
    }

    public int getTotalTaxaNumber(){
        return this.taxaList.size();
    }

    public List<String> getPop(){
        Set<String> popSet = new HashSet<>();
        for (String pop : this.popList){
            popSet.add(pop);
        }
        List<String> popList = new ArrayList<>(popSet);
        Collections.sort(popList);
        return popList;
    }

    public List<Source> getSourceList() {
        return sourceList;
    }

    public List<String> getTaxaOf(Source source){
        List<String> taxaList= new ArrayList<>();
        List<Source> sources = this.getSourceList();
        assert sources.contains(source) : "error, check your popName";
        int totalTaxaNum = this.getTotalTaxaNumber();
        for (int i = 0; i < totalTaxaNum; i++) {
            if (this.getSourceList().get(i).equals(source)){
                taxaList.add(this.taxaList.get(i));
            }
        }
        return taxaList;
    }

    public List<String> getTaxaOf_native_introgressed(){
        List<Source> native_introgressed_sources = this.get_Native_Introgressed_Source();
        List<String> taxaList = new ArrayList<>();
        for (Source source : native_introgressed_sources){
            taxaList.addAll(this.getTaxaOf(source));
        }
        return taxaList;
    }

    public List<String>[] getTaxaOf(List<Source> sources){
        List<String>[] res = new List[sources.size()];
        for (int i = 0; i < sources.size(); i++) {
            res[i] = this.getTaxaOf(sources.get(i));
        }
        return res;
    }

    public List<String>[] getIntrogressedPopTaxa(){
        Set<Source> sourceSet = new HashSet<>(this.sourceList);
        sourceSet.remove(Source.ADMIXED);
        sourceSet.remove(Source.NATIVE);
        List<Source> sources = new ArrayList<>(sourceSet);
        Collections.sort(sources);
        return this.getTaxaOf(sources);
    }

    public List<String>[] getNative_admixed_introgressed_Taxa(){
        int nativePopNum = 1;
        int admixedTaxaNum = this.getTaxaOf(Source.ADMIXED).size();
        int introgressedPopNum = this.getIntrogressedPopTaxa().length;
        List<String>[] res = new List[nativePopNum+admixedTaxaNum+introgressedPopNum];
        res[0] = this.getTaxaOf(Source.NATIVE);
        for (int i = 1; i < admixedTaxaNum + nativePopNum; i++) {
            res[i] = this.getTaxaOf(Source.ADMIXED).subList(i-1, i);
        }
        for (int i = 0; i < introgressedPopNum; i++) {
            res[i+nativePopNum+admixedTaxaNum]=this.getIntrogressedPopTaxa()[i];
        }
        return res;
    }

    public Map<String, Source> getTaxaSourceMap(List<String> taxaList){
        Map<String, Source> taxaSourceMap = new HashMap<>();
        int taxaIndex;
        for (String taxon : taxaList){
            taxaIndex = this.taxaList.indexOf(taxon);
            taxaSourceMap.put(taxon, this.sourceList.get(taxaIndex));
        }
        return taxaSourceMap;
    }

    public List<Source> get_Native_Introgressed_Source(){
        Set<Source> sourceSet = new HashSet<>(this.sourceList);
        sourceSet.remove(Source.ADMIXED);
        List<Source> sources = new ArrayList<>(sourceSet);
        Collections.sort(sources);
        return sources;
    }



}
