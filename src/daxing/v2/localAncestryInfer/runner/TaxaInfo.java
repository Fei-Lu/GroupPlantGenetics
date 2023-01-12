package daxing.v2.localAncestryInfer.runner;

import com.google.common.collect.Table;
import daxing.common.table.RowTableTool;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TaxaInfo {

    Table<String,String,String> taxon_Pop_PopName;
    List<String> popName;
    Map<String, String> taxaPopNameMap;
    Map<String, List<String>> popTaxaListMap;
    public TaxaInfo(String taxaInfo){
        this.taxon_Pop_PopName = RowTableTool.getTable(taxaInfo, 0 ,1,2);
        this.popName = new ArrayList<>();
        this.taxaPopNameMap = new HashMap<>();
        for (Table.Cell<String,String,String> cell : this.taxon_Pop_PopName.cellSet()){
            this.taxaPopNameMap.put(cell.getRowKey(), cell.getValue());
            this.popName.add(cell.getValue());
        }
        this.popTaxaListMap = new HashMap<>();
        for (String pop : this.popName){
            this.popTaxaListMap.put(pop, new ArrayList<>());
        }
        for (Table.Cell<String,String,String> cell : this.taxon_Pop_PopName.cellSet()){
            this.popTaxaListMap.get(cell.getValue()).add(cell.getRowKey());
        }
    }

    public List<String> getPopName() {
        return popName;
    }

    public Map<String, List<String>> getPopTaxaListMap() {
        return popTaxaListMap;
    }

    public Map<String, String> getTaxaPopNameMap() {
        return taxaPopNameMap;
    }

    public int getPopSampleSize(String popName){
        return this.getPopTaxaListMap().get(popName).size();
    }

    public List<String> getTaxaListOf(String popName){
        return this.popTaxaListMap.get(popName);
    }

    public String getTaxonPop(String taxon){
        return this.taxaPopNameMap.get(taxon);
    }
}
