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

    public List<String> getTaxaOf(String popName){
        List<String> taxaList= new ArrayList<>();
        List<String> popList = this.getPop();
        assert popList.contains(popName) : "error, check your popName";
        int totalTaxaNum = this.getTotalTaxaNumber();
        for (int i = 0; i < totalTaxaNum; i++) {
            if (this.popList.get(i).equals(popName)){
                taxaList.add(this.taxaList.get(i));
            }
        }
        return taxaList;
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


}
