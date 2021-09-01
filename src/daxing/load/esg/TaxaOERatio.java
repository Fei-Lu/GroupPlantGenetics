package daxing.load.esg;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import daxing.common.IOTool;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang.math.NumberUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;

public class TaxaOERatio {

    List<OERatio> oeRatioList;

    public TaxaOERatio(List<OERatio> oeRatioList){
        this.oeRatioList=oeRatioList;
        this.sortByTaxonName();
    }

    public List<OERatio> getOeRatioList() {
        return oeRatioList;
    }

    private void sortByTaxonName(){
        Collections.sort(this.oeRatioList);
    }

    public int getTaxaNum(){
        return this.getOeRatioList().size();
    }

    public void removeIf(Predicate<OERatio> p){
        this.oeRatioList.removeIf(p);
    }

    public int getTaxaIndex(String taxon){
        return Collections.binarySearch(this.oeRatioList, new OERatio(taxon));
    }

    public OERatio getOERatio(int taxonIndex){
        return this.getOeRatioList().get(taxonIndex);
    }

    public TaxaOERatio getSubsetTaxaOERatio(TIntArrayList taxonIndexes){
        List<OERatio> oeRatioList=new ArrayList<>();
        for (int i = 0; i < taxonIndexes.size(); i++) {
            oeRatioList.add(this.getOERatio(taxonIndexes.get(i)));
        }
        return new TaxaOERatio(oeRatioList);
    }

    /**
     *
     * @param nonNATaxaList
     * @return
     */
    public double[] getESGExpectedCount(List<String> nonNATaxaList){
        TIntArrayList taxonIndexes = new TIntArrayList();
        for (String taxon: nonNATaxaList){
            taxonIndexes.add(this.getTaxaIndex(taxon));
        }
        TaxaOERatio subTaxaOERatio = this.getSubsetTaxaOERatio(taxonIndexes);
        double[] expectedESGCount = new double[ESG.values().length];
        double[] pESG;
        for (int i = 0; i < subTaxaOERatio.getTaxaNum(); i++) {
            pESG = subTaxaOERatio.getOERatio(i).getpESG();
            for (int j = 0; j < expectedESGCount.length; j++) {
                expectedESGCount[j] += pESG[j];
            }
        }
        return expectedESGCount;
    }

    public static TaxaOERatio getInstance(String allRatioOEFile){
        long start=System.nanoTime();
        List<OERatio> oeRatioList=null;
        try (BufferedReader br = IOTool.getReader(allRatioOEFile)) {
            String line;
            List<String> temp;
            OERatio oeRatio;
            br.readLine();
            String taxa = null;
            IfPseudoHexaploid ifPseudoHexaploid = null;
            GroupPseudo groupPseudo = null;
            int nonNaTriadsNum=-1;
            double thresholdQuantileProbs = -1;
            double[] thresholdQuantile= null;  // A B D
            double[] actualQuantileProbs= null; // A B D
            List<ESG> esgList= new ArrayList<>();
            TDoubleArrayList pESGList=new TDoubleArrayList();
            TDoubleArrayList expectedList= new TDoubleArrayList();
            TDoubleArrayList observedList = new TDoubleArrayList();
            TDoubleArrayList ratioObservedExpectedList= new TDoubleArrayList();
            ESG[] esgs;
            double[] pESG, expected, observed, ratioObservedExpected;
            oeRatioList= new ArrayList<>();
            while ((line=br.readLine())!=null){
                if (esgList.size()==8){
                    esgs = esgList.stream().toArray(ESG[]::new);
                    pESG = pESGList.toArray();
                    expected = expectedList.toArray();
                    observed = observedList.toArray();
                    ratioObservedExpected = ratioObservedExpectedList.toArray();
                    oeRatio= new OERatio(taxa, ifPseudoHexaploid, groupPseudo, nonNaTriadsNum,
                            thresholdQuantileProbs, thresholdQuantile, actualQuantileProbs, esgs, pESG, expected,
                            observed, ratioObservedExpected);
                    oeRatioList.add(oeRatio);
                    esgList=new ArrayList<>();
                    pESGList= new TDoubleArrayList();
                    expectedList= new TDoubleArrayList();
                    observedList = new TDoubleArrayList();
                    ratioObservedExpectedList= new TDoubleArrayList();
                    temp = PStringUtils.fastSplit(line);
                    taxa = temp.get(0);
                    ifPseudoHexaploid = Integer.parseInt(temp.get(1)) == 1 ? IfPseudoHexaploid.PSEUDO_HEXAPLOID : IfPseudoHexaploid.HEXAPLOID;
                    groupPseudo = GroupPseudo.valueOf(temp.get(2));
                    nonNaTriadsNum=Integer.parseInt(temp.get(3));
                    thresholdQuantileProbs = Double.parseDouble(temp.get(4));
                    thresholdQuantile = new double[3];
                    thresholdQuantile[0]= NumberUtils.isNumber(temp.get(5)) ? Double.parseDouble(temp.get(5)) : Double.POSITIVE_INFINITY;
                    thresholdQuantile[1] = NumberUtils.isNumber(temp.get(6)) ? Double.parseDouble(temp.get(6)) : Double.POSITIVE_INFINITY;
                    thresholdQuantile[2] = NumberUtils.isNumber(temp.get(7)) ? Double.parseDouble(temp.get(7)) : Double.POSITIVE_INFINITY;
                    actualQuantileProbs = new double[3];
                    actualQuantileProbs[0] = NumberUtils.isNumber(temp.get(8)) ? Double.parseDouble(temp.get(8)) : Double.POSITIVE_INFINITY;
                    actualQuantileProbs[1] = NumberUtils.isNumber(temp.get(9)) ? Double.parseDouble(temp.get(9)) : Double.POSITIVE_INFINITY;
                    actualQuantileProbs[2] = NumberUtils.isNumber(temp.get(10)) ? Double.parseDouble(temp.get(10)) : Double.POSITIVE_INFINITY;
                    esgList.add(ESG.valueOf(temp.get(11)));
                    pESGList.add(Double.parseDouble(temp.get(12)));
                    expectedList.add(Double.parseDouble(temp.get(13)));
                    observedList.add(Double.parseDouble(temp.get(14)));
                    ratioObservedExpectedList.add(Double.parseDouble(temp.get(15)));
                }else {
                    temp = PStringUtils.fastSplit(line);
                    taxa = temp.get(0);
                    ifPseudoHexaploid = Integer.parseInt(temp.get(1)) == 1 ? IfPseudoHexaploid.PSEUDO_HEXAPLOID : IfPseudoHexaploid.HEXAPLOID;
                    groupPseudo = GroupPseudo.valueOf(temp.get(2));
                    nonNaTriadsNum=Integer.parseInt(temp.get(3));
                    thresholdQuantileProbs = Double.parseDouble(temp.get(4));
                    thresholdQuantile = new double[3];
                    thresholdQuantile[0]= NumberUtils.isNumber(temp.get(5)) ? Double.parseDouble(temp.get(5)) : Double.POSITIVE_INFINITY;
                    thresholdQuantile[1] = NumberUtils.isNumber(temp.get(6)) ? Double.parseDouble(temp.get(6)) : Double.POSITIVE_INFINITY;
                    thresholdQuantile[2] = NumberUtils.isNumber(temp.get(7)) ? Double.parseDouble(temp.get(7)) : Double.POSITIVE_INFINITY;
                    actualQuantileProbs = new double[3];
                    actualQuantileProbs[0] = NumberUtils.isNumber(temp.get(8)) ? Double.parseDouble(temp.get(8)) : Double.POSITIVE_INFINITY;
                    actualQuantileProbs[1] = NumberUtils.isNumber(temp.get(9)) ? Double.parseDouble(temp.get(9)) : Double.POSITIVE_INFINITY;
                    actualQuantileProbs[2] = NumberUtils.isNumber(temp.get(10)) ? Double.parseDouble(temp.get(10)) : Double.POSITIVE_INFINITY;
                    esgList.add(ESG.valueOf(temp.get(11)));
                    pESGList.add(Double.parseDouble(temp.get(12)));
                    expectedList.add(Double.parseDouble(temp.get(13)));
                    observedList.add(Double.parseDouble(temp.get(14)));
                    ratioObservedExpectedList.add(Double.parseDouble(temp.get(15)));
                }
            }
            esgs = esgList.stream().toArray(ESG[]::new);
            pESG = pESGList.toArray();
            expected = expectedList.toArray();
            observed = observedList.toArray();
            ratioObservedExpected = ratioObservedExpectedList.toArray();
            oeRatio= new OERatio(taxa, ifPseudoHexaploid, groupPseudo, nonNaTriadsNum,
                    thresholdQuantileProbs, thresholdQuantile, actualQuantileProbs, esgs, pESG, expected,
                    observed, ratioObservedExpected);
            oeRatioList.add(oeRatio);
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("new TaxaOERatio instance spend "+ Benchmark.getTimeSpanSeconds(start)+ " second");
        return new TaxaOERatio(oeRatioList);
    }

    public static class OERatio implements Comparable<OERatio> {

        String taxa;
        IfPseudoHexaploid ifPseudoHexaploid;
        GroupPseudo groupPseudo;
        int nonNaTriadsNum;
        double thresholdQuantileProbs;
        double[] thresholdQuantile;  // A B D
        double[] actualQuantileProbs; // A B D
        ESG[] esgs;
        double[] pESG;
        double[] expected;
        double[] observed;
        double[] ratioObservedExpected;

        public OERatio(String taxa, IfPseudoHexaploid ifPseudoHexaploid, GroupPseudo groupPseudo, int nonNaTriadsNum,
                       double thresholdQuantileProbs, double[] thresholdQuantile, double[] actualQuantileProbs,
                       ESG[] esgs, double[] pESG, double[] expected, double[] observed, double[] ratioObservedExpected){
            this.taxa=taxa;
            this.ifPseudoHexaploid=ifPseudoHexaploid;
            this.groupPseudo=groupPseudo;
            this.nonNaTriadsNum=nonNaTriadsNum;
            this.thresholdQuantileProbs=thresholdQuantileProbs;
            this.thresholdQuantile=thresholdQuantile;
            this.actualQuantileProbs=actualQuantileProbs;
            this.esgs=esgs;
            this.pESG=pESG;
            this.expected=expected;
            this.observed=observed;
            this.ratioObservedExpected=ratioObservedExpected;
            this.quickSort();
        }

        public OERatio(String taxa){
            this.taxa=taxa;
            this.ifPseudoHexaploid=null;
            this.groupPseudo=null;
            this.nonNaTriadsNum=-1;
            this.thresholdQuantileProbs=-1;
            this.thresholdQuantile=null;
            this.actualQuantileProbs=null;
            this.esgs=null;
            this.pESG=null;
            this.expected=null;
            this.observed=null;
            this.ratioObservedExpected=null;
        }

        public String getTaxa() {
            return taxa;
        }

        public IfPseudoHexaploid isIfPseudoHexaploid() {
            return ifPseudoHexaploid;
        }

        public double getThresholdQuantileProbs() {
            return thresholdQuantileProbs;
        }

        public double[] getActualQuantileProbs() {
            return actualQuantileProbs;
        }

        public double[] getExpected() {
            return expected;
        }

        public double[] getObserved() {
            return observed;
        }

        public double[] getpESG() {
            return pESG;
        }

        public double[] getRatioObservedExpected() {
            return ratioObservedExpected;
        }

        public double[] getThresholdQuantile() {
            return thresholdQuantile;
        }

        public ESG[] getEsgs() {
            return esgs;
        }

        public GroupPseudo getGroupPseudo() {
            return groupPseudo;
        }

        public int getNonNaTriadsNum() {
            return nonNaTriadsNum;
        }

        @Override
        public int compareTo(OERatio o) {
            return this.getTaxa().compareTo(o.getTaxa());
        }

        Swapper swapper = new Swapper() {
            @Override
            public void swap(int a, int b) {
                ESG esg = esgs[a];
                esgs[a] = esgs[b];
                esgs[b] = esg;
                double p = pESG[a];
                pESG[a] = pESG[b];
                pESG[b] = p;
                double e = expected[a];
                expected[a] = expected[b];
                expected[b] = e;
                double o = observed[a];
                observed[a] = observed[b];
                observed[b] = o;
                double r = ratioObservedExpected[a];
                ratioObservedExpected[a] = ratioObservedExpected[b];
                ratioObservedExpected[b] = r;
            }
        };

        IntComparator intComparator = new IntComparator() {
            @Override
            public int compare(int a, int b) {
                if (esgs[a].getValue() > esgs[b].getValue()){
                    return 1;
                }else if (esgs[a].getValue() < esgs[b].getValue()){
                    return -1;
                }
                return 0;
            }
        };

        private void quickSort(){
            GenericSorting.quickSort(0, this.esgs.length, intComparator, swapper);
        }

    }
}
