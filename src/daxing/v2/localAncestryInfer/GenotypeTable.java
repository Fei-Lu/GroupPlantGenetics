package daxing.v2.localAncestryInfer;

import daxing.common.bisnp.SNP;
import daxing.common.chrrange.ChrPos;
import daxing.common.utiles.IOTool;
import it.unimi.dsi.fastutil.ints.*;
import pgl.PGLConstraints;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import javax.swing.text.html.HTMLDocument;
import java.io.BufferedReader;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * 记得改BitSet loop
 */
public class GenotypeTable {

    String[] taxa = null;
    SNP[] snps = null;

    /**
     * Bit genotype by site, the first dimension is site; the second dimension is genotype(haploid), and missing
     */
    BitSet[][] genoSite;

    /**
     * Bit genotype by taxon, the first dimension is taxon; the second dimension is genotype(haploid), and missing
     */
    BitSet[][] genoTaxon;

    /**
     * for optium solution
     */
    static int iteration=0;

    public GenotypeTable(String haploidGenotypeFile){
        try {
            List<String> vcfAnnotationList = new ArrayList<>();
            BufferedReader br = IOTool.getReader(haploidGenotypeFile);
            String temp;
            while ((temp = br.readLine()).startsWith("##")) {
                vcfAnnotationList.add(temp);
            }
            List<String> l;
            l = PStringUtils.fastSplit(temp);
            this.taxa = new String[l.size()-9];
            for (int i = 9; i < l.size(); i++) {
                taxa[i-9] = l.get(i);
            }
            ExecutorService pool = Executors.newFixedThreadPool(PGLConstraints.parallelLevel);
            List<Future<BlockLines>> resultList = new ArrayList<>();
            int siteCount = 0;
            int startIndex = 0;
            List<String> lines = new ArrayList ();
            StringBuilder sb = new StringBuilder();
            NumberFormat numberFormat= NumberFormat.getNumberInstance();
            numberFormat.setGroupingUsed(true);
            while ((temp = br.readLine()) != null) {
                lines.add(temp);
                if (lines.size()%BlockLines.blockSize == 0) {
                    BlockLines gsb = new BlockLines(lines, startIndex);
                    Future<BlockLines> result = pool.submit(gsb);
                    resultList.add(result);
                    startIndex+=lines.size();
                    lines = new ArrayList<>();
                }
                siteCount++;
                if (siteCount%1000000 == 0) {
                    sb.setLength(0);
                    sb.append("Read in ").append(siteCount).append(" SNPs from ").append(haploidGenotypeFile);
                    System.out.println(sb);
                }
            }
            br.close();
            if (lines.size() != 0) {
                BlockLines gsb = new BlockLines(lines, startIndex);
                Future<BlockLines> result = pool.submit(gsb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.genoSite = new BitSet[siteCount][];
            this.snps = new SNP[siteCount];
            for (Future<BlockLines> blockVCFFuture : resultList) {
                BlockLines block = blockVCFFuture.get();
                for (int j = 0; j < block.getActBlockSize(); j++) {
                    this.snps[block.getStartIndex() + j] = block.getSNPBlock()[j];
                    this.genoSite[block.getStartIndex() + j] = block.getGenoSiteBlock()[j];
                }
            }
            sb.setLength(0);
            sb.append("A total of ").append(numberFormat.format(this.getSiteNumber())).append(" SNPs are in ").append(haploidGenotypeFile).append("\n");
            sb.append("Genotype table is successfully built");
            System.out.println(sb);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.transposeSiteToTaxon();
//        this.sortByTaxa();
    }

    public SNP[] getSnps() {
        return snps;
    }

    public BitSet[][] getGenoSite() {
        return genoSite;
    }

    public String[] getTaxa() {
        return taxa;
    }

    public int getTaxaNumber() {
        return this.taxa.length;
    }

    public int getSiteNumber() {
        return this.snps.length;
    }

    public int getPosition(int siteIndex) {
        return this.snps[siteIndex].getPos();
    }

    public int getTaxonIndex(String taxon) {
        return Arrays.binarySearch(taxa, taxon);
    }

    public int getSiteIndex(String chromosome, int position) {
        ChrPos query = new ChrPos(chromosome, position);
        return Arrays.binarySearch(this.snps, query);
    }

    private void transposeSiteToTaxon () {
        long start = System.nanoTime();
        genoTaxon = new BitSet[this.getTaxaNumber()][];
        for (int i = 0; i < genoTaxon.length; i++) {
            genoTaxon[i]=new BitSet[2];
            for (int j = 0; j < genoTaxon[i].length; j++) {
                genoTaxon[i][j]=new BitSet(this.getSiteNumber());
            }
        }
        for (int i = 0; i < genoSite.length; i++) {
            for (int j = 0; j < genoSite[i].length; j++) {
                for (int k = 0; k < genoSite[i][j].size(); k++) {
                    if (genoSite[i][j].get(k)){
                        genoTaxon[k][j].set(i);
                    }
                }
            }
        }
        System.out.println("Transpose genoSite to genoTaxon takes " + Benchmark.getTimeSpanSeconds(start) + " seconds.");
    }

    public void sortByTaxa() {
        System.out.println("Start sorting genotype table by taxon");
        long start = System.nanoTime();
        int[] indices = PArrayUtils.getIndicesByAscendingValue(this.taxa);
        Arrays.sort(this.taxa);
        BitSet[][] nGenoTaxon = new BitSet[this.getTaxaNumber()][];
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            nGenoTaxon[i] = this.genoTaxon[indices[i]];
        }
        this.genoTaxon = nGenoTaxon;
        System.out.println("Sorting finished in " + Benchmark.getTimeSpanSeconds(start) + " seconds.");
        this.transposeTaxonToSite();
    }

    public void transposeTaxonToSite () {
        long start = System.nanoTime();
        genoSite = new BitSet[this.getSiteNumber()][];
        for (int i = 0; i < genoSite.length; i++) {
            genoSite[i]=new BitSet[2];
            for (int j = 0; j <genoSite[i].length; j++) {
                genoSite[i][j]=new BitSet(this.getTaxaNumber());
            }
        }
        for (int i = 0; i < genoTaxon.length; i++) {
            for (int j = 0; j < genoTaxon[i].length; j++) {
                for (int k = 0; k < genoTaxon[i][j].size(); k++) {
                    if (genoTaxon[i][j].get(k)){
                        genoSite[k][j].set(i);
                    }
                }
            }
        }
        System.out.println("Transpose genoTaxon to genoSite takes " + Benchmark.getTimeSpanSeconds(start) + " seconds.");
    }

    /**
     *
     * @param taxaIndices taxon可能是不连续的, 也可能是未排序的
     * @param siteIndices 是连续的range, int[2] [start,end] [include,include]
     * @return double[][], the first dimension is taxa, the second dimension is genotype (missing = 0.5, alt=1, ref=0)
     */
    public double[][] getSrcGenotypeFrom(int[] taxaIndices, int[] siteIndices){
        BitSet[][] genoTaxa = new BitSet[taxaIndices.length][];
        for (int i = 0; i < taxaIndices.length; i++) {
            genoTaxa[i] = new BitSet[2];
            System.arraycopy(this.genoTaxon[taxaIndices[i]], 0, genoTaxa[i], 0, genoTaxa[0].length);
        }
        BitSet[][] genoTaxaSite=new BitSet[taxaIndices.length][2];
        for (int i = 0; i < genoTaxa.length; i++) {
            for (int j = 0; j < genoTaxa[i].length; j++) {
                genoTaxaSite[i][j]=genoTaxa[i][j].get(siteIndices[0],siteIndices[1]);
            }
        }
        double[][] genotype=new double[taxaIndices.length][];
        for (int i = 0; i < genotype.length; i++) {
            genotype[i]= new double[siteIndices[1]-siteIndices[0]+1];
            Arrays.fill(genotype[i], 0);
            for (int j = 0; j < genoTaxaSite[i][0].length(); j++) {
                if (genoTaxaSite[i][0].get(j)){
                    genotype[i][j]=1;
                }
            }
            for (int j = 0; j < genoTaxaSite[i][1].length(); j++) {
                if (genoTaxaSite[i][1].get(j)){
                    genotype[i][j]=0.5;
                }
            }
        }
        return genotype;
    }

    /**
     *
     * @param taxaIndices taxon可能是不连续的, 也可能是未排序的
     * @param siteIndices 是连续的range, int[2] [start,end] [include,include]
     * @return BitSet[], the first dimension is taxa, the bitset dimension is genotype
     */
    public BitSet[] getSrcGenotypeBitSetFrom(int[] taxaIndices, int[] siteIndices){
        BitSet[][] genoTaxa = new BitSet[taxaIndices.length][];
        for (int i = 0; i < taxaIndices.length; i++) {
            genoTaxa[i] = new BitSet[2];
            System.arraycopy(this.genoTaxon[taxaIndices[i]], 0, genoTaxa[i], 0, genoTaxa[0].length);
        }
        BitSet[] genoTaxaSite=new BitSet[taxaIndices.length];
        for (int i = 0; i < genoTaxa.length; i++) {
            genoTaxaSite[i] = genoTaxa[i][0].get(siteIndices[0],siteIndices[1]);
        }
        return genoTaxaSite;
    }

    /**
     *
     * @param taxonIndex taxon index
     * @param siteIndices 是连续的range, int[2] [start,end] [include,include]
     * @return double[], genotype (missing = 0.5, alt=1, ref=0)
     */
    public double[] getQueryGenotypeFrom(int taxonIndex, int[] siteIndices){
        BitSet[] genoTaxa = new BitSet[2];
        System.arraycopy(this.genoTaxon[taxonIndex], 0, genoTaxa, 0, genoTaxa.length);
        BitSet[] genoTaxaSite=new BitSet[2];
        for (int i = 0; i < genoTaxa.length; i++) {
            genoTaxaSite[i]=genoTaxa[i].get(siteIndices[0],siteIndices[1]);
        }
        double[] query=new double[siteIndices[1]-siteIndices[0]+1];
        Arrays.fill(query, 0);
        for (int i = 0; i < genoTaxaSite[0].length(); i++) {
            if (genoTaxaSite[0].get(i)){
                query[i] = 1;
            }
        }
        for (int i = 0; i < genoTaxaSite[1].length(); i++) {
            if (genoTaxaSite[1].get(i)){
                query[i]=0.5;
            }
        }
        return query;
    }

    /**
     *
     * @param taxonIndex taxon index
     * @param siteIndices 是连续的range, int[2] [start,end] [include,include]
     * @return double[], genotype (missing = 0.5, alt=1, ref=0)
     */
    public BitSet getQueryGenotypeBitSetFrom(int taxonIndex, int[] siteIndices){
        BitSet[] genoTaxa = new BitSet[2];
        System.arraycopy(this.genoTaxon[taxonIndex], 0, genoTaxa, 0, genoTaxa.length);
        BitSet[] genoTaxaSite=new BitSet[2];
        for (int i = 0; i < genoTaxa.length; i++) {
            genoTaxaSite[i]=genoTaxa[i].get(siteIndices[0],siteIndices[1]);
        }
        return genoTaxaSite[0];
    }

    public String[] getSrcTaxaFrom(int[] taxaIndices){
        String[] srcTaxa= new String[taxaIndices.length];
        String[] taxa=this.getTaxa();
        for (int i = 0; i < taxaIndices.length; i++) {
            srcTaxa[i]=taxa[taxaIndices[i]];
        }
        return srcTaxa;
    }

    /**
     *
     * @param srcGenotype the first dim is haplotype, the second dim is SNP position
     * @param queryGenotype
     * @param switchCostScore
     * @return mini cost score matrix, the first dim is haplotype, the second dim is SNP position
     */
    public static double[][] getMiniCostScore(double[][] srcGenotype, double[] queryGenotype,
                                              double switchCostScore){
        //        double switchCostScore= 1.5;
//        int[][] srcGenotype = {{0,1,0,1,0,1,0,0,0,0,1,1},
//                            {0,0,0,1,0,1,1,0,0,0,1,1},
//                            {0,0,1,0,1,0,0,0,1,0,1,1},
//                            {0,0,0,0,1,0,1,0,1,1,1,1},
//                            {1,1,0,0,0,0,1,1,1,1,0,0},
//                            {1,0,0,1,0,0,1,1,1,1,0,0}};
//        int[] queryGenotype =       {1,1,0,0,0,1,0,0,1,1,1,1};

//        long start = System.nanoTime();

        // distance
        double[][] distance = new double[srcGenotype.length][];
        for (int i = 0; i < distance.length; i++) {
            distance[i]= new double[srcGenotype[i].length];
            Arrays.fill(distance[i], -1);
        }
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                distance[i][j]=Math.abs(srcGenotype[i][j]- queryGenotype[j]);
            }
        }

        // initialize mini cost score
        double[][] miniCost = new double[distance.length][];
        for (int i = 0; i < miniCost.length; i++) {
            miniCost[i] = new double[distance[0].length];
            Arrays.fill(miniCost[i], -1);
            miniCost[i][0] = distance[i][0];
        }


        // i is SNP position
        // j is haplotype index of source population
        // miniCost
        for (int i = 1; i < distance[0].length; i++) {

            // j-1 SNP位置，单倍型路径发生switch对应的最小Cost
            double miniCostSwitch=Double.MAX_VALUE;
            for (int j = 0; j < distance.length; j++) {
                miniCostSwitch = miniCost[j][i-1] < miniCostSwitch ? miniCost[j][i-1] : miniCostSwitch;
            }

            for (int j = 0; j < distance.length; j++) {
                // 最小cost路径对应当前haplotype
                if (miniCost[j][i-1] < miniCostSwitch+switchCostScore){
                    miniCost[j][i] = miniCost[j][i-1] + distance[j][i];
                }else {
                    // 最小cost路径对应转换单倍型
                    miniCost[j][i] = miniCostSwitch+switchCostScore+distance[j][i];
                }
            }
        }
//        System.out.println("calculate mini cost matrix take "+Benchmark.getTimeSpanSeconds(start)+ " seconds");
        return miniCost;
    }

    /**
     * 考虑单倍型所属的群体, 进行罚分
     * @param srcGenotype the first dim is haplotype, the second dim is SNP position
     * @param queryGenotype genotype
     * @param switchCostScore λ
     * @param srcIndiList 单倍型的样本名称
     * @param taxaSourceMap 样本对应群体的映射
     * @return mini cost score matrix, the first dim is haplotype, the second dim is SNP position
     */
    public static double[][] getMiniCostScore(double[][] srcGenotype, double[] queryGenotype,
                                                            double switchCostScore,
                                                            List<String> srcIndiList,
                                                            Map<String, WindowSource.Source> taxaSourceMap){
//        long start = System.nanoTime();
        // distance
        double[][] distance = new double[srcGenotype.length][];
        for (int i = 0; i < distance.length; i++) {
            distance[i]= new double[srcGenotype[i].length];
            Arrays.fill(distance[i], -1);
        }
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                distance[i][j]=Math.abs(srcGenotype[i][j]- queryGenotype[j]);
            }
        }

        // initialize mini cost score
        double[][] miniCost = new double[distance.length][];
        for (int i = 0; i < miniCost.length; i++) {
            miniCost[i] = new double[distance[0].length];
            Arrays.fill(miniCost[i], -1);
            miniCost[i][0] = distance[i][0];
        }


        // i is SNP position
        // j is haplotype index of source population
        // miniCost

        IntSet miniCostSwitchIndexSet;
        for (int i = 1; i < distance[0].length; i++) {

            // j-1 SNP位置，单倍型路径发生switch对应的最小Cost
            double miniCostSwitch=Double.MAX_VALUE;
            for (int j = 0; j < distance.length; j++) {
                miniCostSwitch = miniCost[j][i-1] < miniCostSwitch ? miniCost[j][i-1] : miniCostSwitch;
            }

            // j-1 SNP位置 单倍型路径发生switch对应的最小Cost index
            miniCostSwitchIndexSet = new IntOpenHashSet();
            for (int j = 0; j < distance.length; j++) {
                if (miniCostSwitch==miniCost[j][i-1]){
                    miniCostSwitchIndexSet.add(j);
                }
            }

            for (int j = 0; j < distance.length; j++) {
                // 最小cost路径对应当前haplotype
                if (miniCost[j][i-1] < miniCostSwitch+switchCostScore){
                    miniCost[j][i] = miniCost[j][i-1] + distance[j][i];
                }else {
                    // 最小cost路径对应转换单倍型

                    // i-1 SNP和当前SNP所属的单倍型是否来自同一群体
                    // 转换单倍型中i-1 SNP和当前SNP是否有在同一群体的案例
                    boolean ifExistAnyPairBelongToSamePop=false;
                    for (int index:miniCostSwitchIndexSet){
                        if(taxaSourceMap.get(srcIndiList.get(index)).equals(taxaSourceMap.get(srcIndiList.get(j)))){
                            ifExistAnyPairBelongToSamePop = true;
                            break;
                        }
                    }

                    // || taxaSourceMap.get(srcIndiList.get(j)).equals(WindowSource.Source.NONE)
                    if (ifExistAnyPairBelongToSamePop){
                        // i-1 SNP和当前SNP i所属的单倍型有来自同一群体 λ
                        // λ
                        miniCost[j][i] = miniCostSwitch+switchCostScore+distance[j][i];
                    }else {
                        // i-1 SNP和当前SNP所属的单倍型来自不同群体
                        // λ+1
                        miniCost[j][i] = miniCostSwitch+(switchCostScore+1)+distance[j][i];
                    }

                }
            }
        }

//        System.out.println("calculate mini cost matrix take "+Benchmark.getTimeSpanSeconds(start)+ " seconds");

        return miniCost;
    }

    /**
     *
     * @param miniCost the first dim is haplotype, the second dim is SNP position
     * @param switchCostScore λ
     * @return candidate solution, the first dim is mini cost sum haplotype index, the second dim is SNP position
     */
    public static IntSet[][] getCandidateSolution(double[][] miniCost, double switchCostScore){

//        long start = System.nanoTime();
        int haplotypeLen=miniCost[0].length;

        // miniCost score indexList
        double miniCostScore= Double.MAX_VALUE;
        IntList miniCostScoreIndexList = new IntArrayList();
        for (int i = 0; i < miniCost.length; i++) {
            miniCostScore = miniCost[i][haplotypeLen-1] < miniCostScore ? miniCost[i][haplotypeLen-1] :miniCostScore;
        }
        for (int i = 0; i < miniCost.length; i++) {
            if (miniCostScore==miniCost[i][haplotypeLen-1]){
                miniCostScoreIndexList.add(i);
            }
        }

        // new solution, the first dim is miniCostScoreList, the second dim is SNP
        IntSet[][] solution = new IntSet[miniCostScoreIndexList.size()][];
        for (int i = 0; i < solution.length; i++) {
            solution[i] = new IntSet[haplotypeLen];
            for (int j = 0; j < solution[i].length; j++) {
                solution[i][j] = new IntOpenHashSet();
            }
        }
        for (int i = 0; i < solution.length; i++) {
            solution[i][haplotypeLen-1].add(miniCostScoreIndexList.getInt(i));
        }


        // find all solution
        IntSet currentIndexSet;
        IntIterator tIntIterator;
        int index;
        for (int i = 0; i < solution.length; i++) {
            index = miniCostScoreIndexList.getInt(i);
            currentIndexSet = new IntOpenHashSet();
            currentIndexSet.add(index);
            for (int j = haplotypeLen - 1; j > 0; j--) {
                tIntIterator = currentIndexSet.iterator();
                while (tIntIterator.hasNext()){
                    index = tIntIterator.nextInt();

                    // 当前单倍型
                    if (miniCost[index][j-1] <= miniCost[index][j]){
                        solution[i][j-1].add(index);
                    }

                    // 转换单倍型
                    for (int k = 0; k < miniCost.length; k++) {
                        if (k==index) continue;
                        if ((miniCost[k][j-1]+switchCostScore) <= miniCost[index][j]){
                            solution[i][j-1].add(k);
                        }
                    }
                }
                currentIndexSet = solution[i][j-1];
            }
        }

//        System.out.println("calculate candidate solution take "+Benchmark.getTimeSpanSeconds(start)+" seconds");
        return solution;
    }


    /**
     *
     * @param miniCost the first dim is haplotype, the second dim is SNP position
     * @param switchCostScore λ
     * @param srcIndiList 单倍型的样本名称
     * @param taxaSourceMap 样本对应群体的映射
     * @return candidate source solution, the first dim is mini cost sum haplotype index, the second dim is SNP position
     */
    public static List<WindowSource.Source>[][] getCandidateSourceSolution(double[][] miniCost, double switchCostScore,
                                                  List<String> srcIndiList,
                                                  Map<String, WindowSource.Source> taxaSourceMap){
//        long start = System.nanoTime();
        int haplotypeLen=miniCost[0].length;

        IntSet[][] solution = GenotypeTable.getCandidateSolution(miniCost, switchCostScore);

        EnumSet<WindowSource.Source>[][] solutionSource= new EnumSet[solution.length][];
        for (int i = 0; i < solutionSource.length; i++) {
            solutionSource[i] = new EnumSet[haplotypeLen];
            for (int j = 0; j < solutionSource[i].length; j++) {
                solutionSource[i][j] = EnumSet.noneOf(WindowSource.Source.class);
            }
        }

        for (int i = 0; i < solutionSource.length; i++) {
            for (int j = 0; j < solutionSource[i].length; j++) {
                for (int ele:solution[i][j]){
                    solutionSource[i][j].add(taxaSourceMap.get(srcIndiList.get(ele)));
                }
            }
        }

        List<WindowSource.Source>[][] solutionSourceList = new List[solution.length][];
        for (int i = 0; i < solutionSourceList.length; i++) {
            solutionSourceList[i] = new List[haplotypeLen];
            for (int j = 0; j < solutionSource[i].length; j++) {
                solutionSourceList[i][j]=new ArrayList<>(solutionSource[i][j]);
                Collections.sort(solutionSourceList[i][j]);
            }
        }


        // 从单倍型solution到Source solution的转变过程中, 会有重复的solution
        Set<List<List<WindowSource.Source>>> set= new HashSet<>();
        List<List<WindowSource.Source>> list;
        for (int i = 0; i < solutionSourceList.length; i++) {
            list = new ArrayList<>();
            for (int j = 0; j < solutionSourceList[i].length; j++) {
                list.add(solutionSourceList[i][j]);
            }
            set.add(list);
        }

        List<List<List<WindowSource.Source>>> listRemovedDup= new ArrayList<>(set);
        List<WindowSource.Source>[][] solutionSourceListRemovedDup= new List[listRemovedDup.size()][];
        for (int i = 0; i < solutionSourceListRemovedDup.length; i++) {
            solutionSourceListRemovedDup[i] = new List[haplotypeLen];
            for (int j = 0; j < solutionSourceListRemovedDup[i].length; j++) {
                solutionSourceListRemovedDup[i][j] = listRemovedDup.get(i).get(j);
                Collections.sort(solutionSourceListRemovedDup[i][j]);
            }
        }

//        System.out.println("calculate candidate solution take "+Benchmark.getTimeSpanSeconds(start)+" seconds");
        return solutionSourceListRemovedDup;
    }

    private static int getSolutionSize(List<WindowSource.Source>[] solution, int count){
        int size=solution[0].size();
        for (int i = 1; i < solution.length; i++) {
            if (solution[i].size() == 1) continue;
            if (solution[i-1].size() == solution[i].size()) continue;
            size *=solution[i].size();
            if (size > (Integer.MAX_VALUE/count)){
                return Integer.MAX_VALUE;
            }
        }
        return size;
    }

    private static int getSolutionSize(List<WindowSource.Source>[][] solution){
        int res=GenotypeTable.getSolutionSize(solution[0], solution.length);
        if (res==Integer.MAX_VALUE) return Integer.MAX_VALUE;
        int size;
        for (int i = 1; i < solution.length; i++) {
            size= GenotypeTable.getSolutionSize(solution[i], solution.length);
            if (size==Integer.MAX_VALUE){
                return Integer.MAX_VALUE;
            }
            res+=size;
            if (res > (Integer.MAX_VALUE/4)){
                return Integer.MAX_VALUE;
            }
        }

        return res;
    }

    private static int[] getSolutionSizeArray(List<WindowSource.Source>[][] solution){
        int[] res = new int[solution.length];
        for (int i = 0; i < solution.length; i++) {
            res[i] = getSolutionSize(solution[i], solution.length);
        }
        return res;
    }

    /**
     *
     * @param srcGenotype the first dimension is haplotype; the second dimension is SNP
     * @param queryGenotype
     * @return
     */
    public static List<WindowSource.Source[]> getMiniPath(double[][] srcGenotype, double[] queryGenotype,
                                                          double switchCostScore,
                                                          List<String> srcIndiList,
                                                          Map<String, WindowSource.Source> taxaSourceMap,
                                                          int maxSolutionCount, double maxSwitchCostScore){

        iteration++;
        // distance
        double[][] miniCost = GenotypeTable.getMiniCostScore(srcGenotype, queryGenotype, switchCostScore);
        List<WindowSource.Source>[][] solutionSource =
                GenotypeTable.getCandidateSourceSolution(miniCost, switchCostScore, srcIndiList, taxaSourceMap);

        // transform solution to array
        List<WindowSource.Source[]> optiumSolutionList = new ArrayList<>();

        WindowSource.Source[] subSolution;
        int cumSize=0;
        for (int m = 0; m < solutionSource.length; m++) {

            int currentSolutionSize, multiplySolutionSize;
            currentSolutionSize = solutionSource[m][0].size();
            multiplySolutionSize = currentSolutionSize;

            for (int i = 0; i < currentSolutionSize; i++) {
                subSolution = new WindowSource.Source[queryGenotype.length];
                Arrays.fill(subSolution, null);
                optiumSolutionList.add(subSolution);
            }

            for (int i = 0; i < currentSolutionSize; i++) {
                optiumSolutionList.get(cumSize+i)[0]=solutionSource[m][0].get(i);
            }

            for (int i = 1; i < solutionSource[m].length; i++) {
                currentSolutionSize = solutionSource[m][i].size();

                // {WE}
                if (currentSolutionSize == 1){
                    for (int j = cumSize; j < optiumSolutionList.size(); j++) {
                        optiumSolutionList.get(j)[i] = solutionSource[m][i].get(0);
                    }

                    // {WE,DE}, {WE,DE}
                }else if (solutionSource[m][i].equals(solutionSource[m][i-1])){
                    for (int j = cumSize; j < optiumSolutionList.size(); j++) {
                        optiumSolutionList.get(j)[i]=optiumSolutionList.get(j)[i-1];
                    }
                }//            else if (CollectionTool.hasIntersection(solution[i], solution[i-1])){
//
//                // intersection
//                IntSet intersectionSet = new IntOpenHashSet(solution[i]);
//                intersectionSet.retainAll(solution[i-1]);
//                for (int ele: intersectionSet){
//                    int eleIndexIMinus1 = Collections.binarySearch(solutionList[i-1], ele);
//                    for (int j = 0; j < optiumSolutionList.size(); j++) {
//                        if (optiumSolutionList.get(j)[i-1]!=ele) continue;
//                        optiumSolutionList.get(j)[i]= solutionList[i-1].getInt(eleIndexIMinus1);
//                    }
//                }
//
//                // removedAll
//                IntSet removedAllSet = new IntOpenHashSet(solution[i]);
//                removedAllSet.removeAll(solution[i-1]);
//                for (int j = 0; j < (multiplySolutionSize*removedAllSet.size()+intersectionSet.size()-multiplySolutionSize); j++) {
//                    subSolution = new int[queryGenotype.length];
//                    Arrays.fill(subSolution, -1);
//                    optiumSolutionList.add(subSolution);
//                }
//                for (int j = 0; j < multiplySolutionSize; j++) {
//                    for (int k = 0; k < (currentSolutionSize-1); k++) {
//                        System.arraycopy(optiumSolutionList.get(j),0,optiumSolutionList.get(multiplySolutionSize*(k+1)+j),
//                                0, i);
//                    }
//                }
//
//                for (int j = 0; j < currentSolutionSize; j++) {
//                    for (int k = 0; k < multiplySolutionSize; k++) {
//                        if ((k+j*multiplySolutionSize)<(intersectionSet.size()*multiplySolutionSize)) continue;
//                        optiumSolutionList.get(k+j*multiplySolutionSize)[i]=solutionList[i].getInt(j);
//                    }
//                }
//                multiplySolutionSize=optiumSolutionList.size();
//            }
                else if (!solutionSource[m][i].equals(solutionSource[m][i-1])){
                    // new
                    for (int j = 0; j < (multiplySolutionSize*currentSolutionSize-multiplySolutionSize); j++) {
                        subSolution = new WindowSource.Source[queryGenotype.length];
                        Arrays.fill(subSolution, null);
                        optiumSolutionList.add(subSolution);
                    }

                    // 递归调用
                    if (optiumSolutionList.size() > maxSolutionCount && switchCostScore < maxSwitchCostScore){
                        System.out.println("iteration "+iteration);
                        System.out.println("Switch cost score is "+switchCostScore);
                        System.out.println();
                        return GenotypeTable.getMiniPath(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
                                taxaSourceMap, maxSolutionCount, maxSwitchCostScore);
                    }

                    // i-1 SNP 赋值
                    for (int j = 0; j < multiplySolutionSize; j++) {
                        for (int k = 0; k < (currentSolutionSize-1); k++) {
                            System.arraycopy(optiumSolutionList.get(cumSize+j),0,
                                    optiumSolutionList.get(multiplySolutionSize*(k+1)+j+cumSize),0, i);
                        }
                    }

                    // i SNP 赋值
                    for (int j = 0; j < currentSolutionSize; j++) {
                        for (int k = 0; k < multiplySolutionSize; k++) {
                            optiumSolutionList.get(k+j*multiplySolutionSize+cumSize)[i]=solutionSource[m][i].get(j);
                        }
                    }
                    multiplySolutionSize *=currentSolutionSize;
                }
            }
            cumSize = optiumSolutionList.size();
        }

        System.out.println("iteration "+iteration);
        System.out.println("Switch cost score is "+switchCostScore);
        System.out.println();
        System.out.println("optium switch cost score is "+switchCostScore+", solution size is "+optiumSolutionList.size());
        iteration=0;
        return optiumSolutionList;
    }

    /**
     *
     * @param srcGenotype the first dimension is haplotype; the second dimension is SNP
     * @param queryGenotype
     * @return
     */
    public static List<List<SolutionElement>> getMiniPath2(double[][] srcGenotype,
                                                                          double[] queryGenotype,
                                                            double switchCostScore,
                                                            List<String> srcIndiList,
                                                            Map<String, WindowSource.Source> taxaSourceMap,
                                                           int maxSolutionCount){

        iteration++;
        double[][] miniCostCurrent, miniCostNext;
        List<WindowSource.Source>[][] solutionSourceCurrent, solutionSourceNext;

        miniCostCurrent = GenotypeTable.getMiniCostScore(srcGenotype, queryGenotype, switchCostScore);
        solutionSourceCurrent = GenotypeTable.getCandidateSourceSolution(miniCostCurrent, switchCostScore, srcIndiList,
                taxaSourceMap);

        miniCostNext =  GenotypeTable.getMiniCostScore(srcGenotype, queryGenotype, switchCostScore+1);
        solutionSourceNext = GenotypeTable.getCandidateSourceSolution(miniCostNext, switchCostScore+1, srcIndiList,
                taxaSourceMap);

        int totalSolutionSizeCurrent = GenotypeTable.getSolutionSize(solutionSourceCurrent);
        int totalSolutionSizeNext = GenotypeTable.getSolutionSize(solutionSourceNext);


        if (totalSolutionSizeNext < totalSolutionSizeCurrent/2){
            System.out.println();
            System.out.println("iteration "+iteration );
            System.out.println("Switch cost score is "+switchCostScore);
            System.out.println("Total solution size is "+totalSolutionSizeCurrent);
            System.out.println();
            return GenotypeTable.getMiniPath2(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
                    taxaSourceMap, maxSolutionCount);
        }

        if (totalSolutionSizeCurrent > maxSolutionCount){
            System.out.println();
            System.out.println("iteration "+iteration);
            System.out.println("Switch cost score is "+switchCostScore);
            System.out.println("Total solution size is "+totalSolutionSizeCurrent);
            System.out.println("Total solution size greater than maxSolutionCount "+maxSolutionCount);
            System.out.println("do not calculate");
            System.out.println();
            iteration=0;
            return null;
        }

        System.out.println();
        System.out.println("iteration "+iteration);
        System.out.println("Switch cost score is "+switchCostScore);
        System.out.println("Total solution size is "+totalSolutionSizeCurrent);


        // transform solution to array
        List<List<SolutionElement>> optiumSolutionList = new ArrayList<>();

        List<SolutionElement> subSolution;
        int cumSize=0;
        WindowSource.Source currentSource;
        SolutionElement solutionElement;
        int start=0, end=1;
        Iterator<SolutionElement> iterator;
        for (int m = 0; m < solutionSourceCurrent.length; m++) {

            int currentSolutionSize, multiplySolutionSize;
            currentSolutionSize = solutionSourceCurrent[m][0].size();
            multiplySolutionSize = currentSolutionSize;


            for (int i = 0; i < currentSolutionSize; i++) {
                subSolution = new ArrayList<>();
                optiumSolutionList.add(subSolution);
            }

            for (int i = 0; i < currentSolutionSize; i++) {
                currentSource = solutionSourceCurrent[m][0].get(i);
                solutionElement = new SolutionElement(currentSource, start, end);
                optiumSolutionList.get(cumSize+i).add(solutionElement);
            }

            for (int i = 1; i < solutionSourceCurrent[m].length; i++) {
                currentSolutionSize = solutionSourceCurrent[m][i].size();

                // {WE}
                if (currentSolutionSize == 1){
                    for (int j = cumSize; j < optiumSolutionList.size(); j++) {
                        subSolution = optiumSolutionList.get(j);

                        // 最后一个solutionElement进行延伸
                        solutionElement = subSolution.get(subSolution.size()-1);
                        solutionElement.extend();
                    }

                    // {WE,DE}, {WE,DE}
                }else if (solutionSourceCurrent[m][i].equals(solutionSourceCurrent[m][i-1])){
                    for (int j = cumSize; j < optiumSolutionList.size(); j++) {
                        subSolution = optiumSolutionList.get(j);
                        // 最后一个solutionElement进行延伸
                        solutionElement = subSolution.get(subSolution.size()-1);
                        solutionElement.extend();
                    }
                }
                else if (!solutionSourceCurrent[m][i].equals(solutionSourceCurrent[m][i-1])){
                    // new
                    for (int j = 0; j < (multiplySolutionSize*currentSolutionSize-multiplySolutionSize); j++) {
                        subSolution = new ArrayList<>();
                        optiumSolutionList.add(subSolution);
                    }

//                     // 递归调用
//                    if (optiumSolutionList.size() > 10000 && switchCostScore < maxSwitchCostScore){
//                        System.out.println();
//                        System.out.println("iteration "+iteration);
//                        System.out.println("Switch cost score is "+switchCostScore);
//                        System.out.println();
//                        return GenotypeTable.getMiniPath2(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
//                                taxaSourceMap, maxSolutionCount, maxSwitchCostScore);
//                    }

                    // i-1 SNP 赋值
                    for (int j = 0; j < multiplySolutionSize; j++) {
                        for (int k = 0; k < (currentSolutionSize-1); k++) {
                            iterator = optiumSolutionList.get(cumSize+j).iterator();
                            while (iterator.hasNext()){
                                optiumSolutionList.get(multiplySolutionSize*(k+1)+j+cumSize).add(iterator.next().clone());
                            }
                        }
                    }

//                    // i SNP 赋值
                    for (int j = 0; j < currentSolutionSize; j++) {
                        for (int k = 0; k < multiplySolutionSize; k++) {
                            subSolution = optiumSolutionList.get(k+j*multiplySolutionSize+cumSize);
                            currentSource = solutionSourceCurrent[m][i].get(j);
                            if (subSolution.get(subSolution.size()-1).source==currentSource){
                                subSolution.get(subSolution.size()-1).extend();
                            }else {
                                solutionElement = new SolutionElement(currentSource, i, i+1);
                                subSolution.add(solutionElement);
                            }
                        }
                    }
                    multiplySolutionSize *=currentSolutionSize;
                }
            }
            cumSize = optiumSolutionList.size();
        }

        // 由solutionSource生成的optiumSolutionList会存在相同的路径
        Set<List<SolutionElement>> set = new HashSet<>(optiumSolutionList);
        List<List<SolutionElement>> list = new ArrayList<>(set);

//        if (list.size() > maxSolutionCount && switchCostScore < maxSwitchCostScore){
//            System.out.println();
//            System.out.println("iteration "+iteration);
//            System.out.println("Switch cost score is "+switchCostScore);
//            System.out.println();
//            return GenotypeTable.getMiniPath2(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
//                    taxaSourceMap, maxSolutionCount, maxSwitchCostScore);
//        }

        System.out.println();
        System.out.println("optium switch cost score is "+switchCostScore+", solution count after removing redundancy" +
                " is "+list.size());
        iteration =0;
        return list;
    }

}
