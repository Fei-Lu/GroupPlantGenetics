package daxing.v2.localAncestryInfer;

import gnu.trove.list.array.TIntArrayList;

import java.util.List;
import java.util.Map;

public class LocalSolutions {

    String taxonID;
    List<TIntArrayList> solutions;
    List<String> srcTaxaList;
    int[] siteIndices;
    static Map<String, WindowSource.Source> taxaSourceMap;

    public LocalSolutions(String taxonID, List<TIntArrayList> solutions, List<String> srcTaxaList,
                          int[] siteIndices, Map<String, WindowSource.Source> taxaSourceMap){
        this.taxonID=taxonID;
        this.solutions=solutions;
        this.srcTaxaList=srcTaxaList;
        this.siteIndices=siteIndices;
        this.taxaSourceMap=taxaSourceMap;
    }


}
