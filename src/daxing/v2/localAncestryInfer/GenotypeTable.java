package daxing.v2.localAncestryInfer;

import daxing.common.bisnp.SNP;
import daxing.common.chrrange.ChrPos;
import daxing.common.utiles.IOTool;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import it.unimi.dsi.fastutil.ints.*;
import pgl.PGLConstraints;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
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
     * @param srcGenotype the first dimension is haplotype; the second dimension is SNP
     * @param queryGenotype
     * @return
     */
    public static List<TIntList> getMiniPath(double[][] srcGenotype, double[] queryGenotype, double switchCostScore){
//        double switchCostScore= 1.5;
//        int[][] srcGenotype = {{0,1,0,1,0,1,0,0,0,0,1,1},
//                            {0,0,0,1,0,1,1,0,0,0,1,1},
//                            {0,0,1,0,1,0,0,0,1,0,1,1},
//                            {0,0,0,0,1,0,1,0,1,1,1,1},
//                            {1,1,0,0,0,0,1,1,1,1,0,0},
//                            {1,0,0,1,0,0,1,1,1,1,0,0}};
//        int[] queryGenotype =       {1,1,0,0,0,1,0,0,1,1,1,1};


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


        // initialize solution of path
        TIntList initializeSolution;
        List<TIntList>[] solutionCurrent = new List[srcGenotype.length];
        List<TIntList>[] solutionPast = new List[srcGenotype.length];
        for (int i = 0; i < solutionCurrent.length; i++) {
            solutionCurrent[i] = new ArrayList<>();
            initializeSolution = new TIntArrayList();
            initializeSolution.add(i);
            solutionCurrent[i].add(initializeSolution);
        }

        for (int i = 0; i < solutionPast.length; i++) {
            solutionPast[i] = new ArrayList<>();
            initializeSolution = new TIntArrayList();
            initializeSolution.add(i);
            solutionPast[i].add(initializeSolution);
        }


        // i is SNP position
        // j is haplotype index of source population
        TIntList switchIndexList;
        TIntList switchHaplotype;
        for (int i = 1; i < distance[0].length; i++) {

            // j-1 SNP位置，单倍型路径发生switch对应的最小Cost
            double miniCostSwitch=Double.MAX_VALUE;
            switchIndexList = new TIntArrayList();
            for (int j = 0; j < distance.length; j++) {
                miniCostSwitch = miniCost[j][i-1] < miniCostSwitch ? miniCost[j][i-1] : miniCostSwitch;
            }
            for (int j = 0; j < distance.length; j++) {
                if (miniCost[j][i-1] == miniCostSwitch){
                    switchIndexList.add(j);
                }
            }

            for (int j = 0; j < distance.length; j++) {
                // 最小cost路径对应当前haplotype
                if (miniCost[j][i-1] < miniCostSwitch+switchCostScore){
                    miniCost[j][i] = miniCost[j][i-1] + distance[j][i];
                    for (int k = 0; k < solutionCurrent[j].size(); k++) {
                        solutionCurrent[j].get(k).add(j);
                    }
                }else {
                    // 最小cost路径对应转换单倍型
                    miniCost[j][i] = miniCostSwitch+switchCostScore+distance[j][i];
                    solutionCurrent[j] = new ArrayList<>();
                    for (int k = 0; k < switchIndexList.size(); k++) {
                        if (switchIndexList.get(k)==j) continue;
                        for (int l = 0; l < solutionPast[switchIndexList.get(k)].size(); l++) {
                            switchHaplotype = new TIntArrayList();
                            switchHaplotype.addAll(solutionPast[switchIndexList.get(k)].get(l));
                            solutionCurrent[j].add(switchHaplotype);
                        }
                    }

                    // 添加当前单倍型所在的索引
                    for (int k = 0; k < solutionCurrent[j].size(); k++) {
                        solutionCurrent[j].get(k).add(j);
                    }

                    // 如果转换单倍型和当前单倍型最小Cost相同，则两种路径同时保留
                    if (miniCost[j][i-1] == miniCostSwitch+switchCostScore){
                        for (int k = 0; k < solutionPast[j].size(); k++) {
                            switchHaplotype = new TIntArrayList();
                            switchHaplotype.addAll(solutionPast[j].get(k));
                            switchHaplotype.add(j);
                            solutionCurrent[j].add(switchHaplotype);
                        }
                    }
                }
            }

            // assign solutionPast to solutionCurrent
            for (int j = 0; j < solutionCurrent.length; j++) {
                solutionPast[j] = new ArrayList<>();
                for (int k = 0; k < solutionCurrent[j].size(); k++) {
                    solutionPast[j].add(new TIntArrayList());
                }
                for (int k = 0; k < solutionCurrent[j].size(); k++) {
                    solutionPast[j].get(k).addAll(solutionCurrent[j].get(k));
                }
            }
        }

        // miniCost index
        double mini = Double.MAX_VALUE;
        for (int i = 0; i < miniCost.length; i++) {
            mini = miniCost[i][miniCost[0].length-1] < mini ? miniCost[i][miniCost[0].length-1]: mini;
        }
        TIntList solutionIndexList = new TIntArrayList();
        for (int i = 0; i < miniCost.length; i++) {
            if (miniCost[i][miniCost[0].length-1]==mini){
                solutionIndexList.add(i);
            }
        }

        // optium solution
        List<TIntList> optium=new ArrayList<>();
        for (int i = 0; i < solutionIndexList.size(); i++) {
            for (int j = 0; j < solutionCurrent[solutionIndexList.get(i)].size(); j++) {
                optium.add(solutionCurrent[solutionIndexList.get(i)].get(j));
            }
        }

        return optium;
    }

    /**
     *
     * @param scrGenotype the first dimension is haplotype; the BitSet dimension is SNP
     * @param queryGenotype
     * @param seqLength
     * @param switchCostScore
     * @return
     */
    public static List<WindowSource.Source[]> getMiniCostPath(BitSet[] scrGenotype, BitSet queryGenotype, int seqLength,
                                                 double switchCostScore, List<String> srcIndiList,
                                                 Map<String, WindowSource.Source> taxaSourceMap,
                                                 int maxSolutionCount, double maxSwitchCostScore) {
//        double switchCostScore = 1.5;
        BitSet[] distance = new BitSet[scrGenotype.length];
        for (int i = 0; i < scrGenotype.length; i++) {
            scrGenotype[i].xor(queryGenotype);
            distance[i] = (BitSet) scrGenotype[i].clone();
        }

        // initialize mini cost score
        double[][] miniCost = new double[distance.length][];
        for (int i = 0; i < miniCost.length; i++) {
            miniCost[i] = new double[seqLength];
            Arrays.fill(miniCost[i], -1);
            miniCost[i][0] = distance[i].get(0) == true ? 1 : 0;
        }

        // i is SNP position
        // j is haplotype index of source population
        // miniCost
        for (int i = 1; i < seqLength; i++) {

            // j-1 SNP位置，单倍型路径发生switch对应的最小Cost
            double miniCostSwitch=Double.MAX_VALUE;
            for (int j = 0; j < distance.length; j++) {
                miniCostSwitch = miniCost[j][i-1] < miniCostSwitch ? miniCost[j][i-1] : miniCostSwitch;
            }

            for (int j = 0; j < distance.length; j++) {
                // 最小cost路径对应当前haplotype
                if (miniCost[j][i-1] < miniCostSwitch+switchCostScore){
                    miniCost[j][i] = miniCost[j][i-1] + (distance[j].get(i)==true ? 1:0);
                }else {
                    // 最小cost路径对应转换单倍型
                    miniCost[j][i] = miniCostSwitch+switchCostScore+(distance[j].get(i)==true ? 1:0);
                }
            }
        }

        // new solution
        IntSet[] solution = new IntSet[miniCost[0].length];
        for (int i = 0; i < solution.length; i++) {
            solution[i] = new IntOpenHashSet();
        }

        // initialize solution
        double currentHaplotypeMiniValue= Double.MAX_VALUE;
        IntSet currentMiniValueIndexSet = new IntOpenHashSet();
        for (int i = 0; i < miniCost.length; i++) {
            currentHaplotypeMiniValue = miniCost[i][seqLength-1] < currentHaplotypeMiniValue ? miniCost[i][seqLength-1] :currentHaplotypeMiniValue;
        }
        for (int i = 0; i < miniCost.length; i++) {
            if (currentHaplotypeMiniValue==miniCost[i][seqLength-1]){
                currentMiniValueIndexSet.add(i);
            }
        }
        solution[seqLength-1].addAll(currentMiniValueIndexSet);

        // find all solution

        IntIterator tIntIterator;
        int index;
        for (int i = miniCost[0].length-1; i > 0; i--) {
            tIntIterator = currentMiniValueIndexSet.iterator();
            while (tIntIterator.hasNext()){
                index = tIntIterator.nextInt();

                // 当前单倍型
                if (miniCost[index][i-1] <= miniCost[index][i]){
                    solution[i-1].add(index);
                }

                // 转换单倍型
                for (int k = 0; k < miniCost.length; k++) {
                    if (k==index) continue;
                    if ((miniCost[k][i-1]+switchCostScore) <= miniCost[index][i]){
                        solution[i-1].add(k);
                    }
                }
            }
            currentMiniValueIndexSet = solution[i-1];
        }
        EnumSet<WindowSource.Source>[] solutionSource= new EnumSet[solution.length];
        for (int i = 0; i < solutionSource.length; i++) {
            solutionSource[i] = EnumSet.noneOf(WindowSource.Source.class);
        }
        for (int i = 0; i < solutionSource.length; i++) {
            for (int ele:solution[i]){
                solutionSource[i].add(taxaSourceMap.get(srcIndiList.get(ele)));
            }
        }
        List<WindowSource.Source>[] solutionSourceList = new List[solutionSource.length];
        for (int i = 0; i < solutionSourceList.length; i++) {
            solutionSourceList[i]=new ArrayList<>(solutionSource[i]);
            Collections.sort(solutionSourceList[i]);
        }

        // transform solution to array
        List<WindowSource.Source[]> optiumSolutionList = new ArrayList<>();
        WindowSource.Source[] subSolution = new WindowSource.Source[seqLength];
        Arrays.fill(subSolution, null);
        optiumSolutionList.add(subSolution);
        int currentSolutionSize, multiplySolutionSize;
        currentSolutionSize = solutionSourceList[0].size();
        multiplySolutionSize = currentSolutionSize;
        for (int i = 0; i < currentSolutionSize-1; i++) {
            subSolution = new WindowSource.Source[seqLength];
            Arrays.fill(subSolution, null);
            optiumSolutionList.add(subSolution);
        }
        for (int i = 0; i < currentSolutionSize; i++) {
            optiumSolutionList.get(i)[0]=solutionSourceList[0].get(i);
        }

        for (int i = 1; i < solutionSourceList.length; i++) {
            currentSolutionSize = solutionSourceList[i].size();

            // {WE}
            if (currentSolutionSize == 1){
                for (int j = 0; j < optiumSolutionList.size(); j++) {
                    optiumSolutionList.get(j)[i] = solutionSourceList[i].get(0);
                }

                // {WE,DE}, {WE,DE}
            }else if (solution[i].equals(solution[i-1])){
                for (int j = 0; j < optiumSolutionList.size(); j++) {
                    optiumSolutionList.get(j)[i]=optiumSolutionList.get(j)[i-1];
                }
            }
//            else if (CollectionTool.hasIntersection(solution[i], solution[i-1])){
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
            else if (!solution[i].equals(solution[i-1])){
                // new
                for (int j = 0; j < (multiplySolutionSize*currentSolutionSize-multiplySolutionSize); j++) {
                    subSolution = new WindowSource.Source[seqLength];
                    Arrays.fill(subSolution, null);
                    optiumSolutionList.add(subSolution);
                }

                // 递归调用
                if (optiumSolutionList.size() > maxSolutionCount && switchCostScore < maxSwitchCostScore){
                    System.out.println("iteration "+iteration);
                    System.out.println("Switch cost score is "+switchCostScore);
                    System.out.println();
                    return GenotypeTable.getMiniCostPath(scrGenotype, queryGenotype, seqLength,switchCostScore+1,
                            srcIndiList,taxaSourceMap, maxSolutionCount, maxSwitchCostScore);
                }

                // i-1 SNP 赋值
                for (int j = 0; j < multiplySolutionSize; j++) {
                    for (int k = 0; k < (currentSolutionSize-1); k++) {
                        System.arraycopy(optiumSolutionList.get(j),0,optiumSolutionList.get(multiplySolutionSize*(k+1)+j),
                                0, i);
                    }
                }

                // i SNP 赋值
                for (int j = 0; j < currentSolutionSize; j++) {
                    for (int k = 0; k < multiplySolutionSize; k++) {
                        optiumSolutionList.get(k+j*multiplySolutionSize)[i]=solutionSourceList[i].get(j);
                    }
                }
                multiplySolutionSize *=currentSolutionSize;
            }
        }
        System.out.println("iteration "+iteration);
        System.out.println("Switch cost score is "+switchCostScore);
        System.out.println();

        System.out.println("optium switch cost score is "+switchCostScore+" solution size is "+optiumSolutionList.size());
        iteration=0;
        return optiumSolutionList;
    }

    /**
     *
     * @param srcGenotype the first dimension is haplotype; the second dimension is SNP
     * @param queryGenotype
     * @return
     */
    public static List<WindowSource.Source[]> getMiniPath2(double[][] srcGenotype, double[] queryGenotype,
                                                           double switchCostScore,
                                                           List<String> srcIndiList,
                                                           Map<String, WindowSource.Source> taxaSourceMap,
                                                           int maxSolutionCount, double maxSwitchCostScore){
//        double switchCostScore= 1.5;
//        int[][] srcGenotype = {{0,1,0,1,0,1,0,0,0,0,1,1},
//                            {0,0,0,1,0,1,1,0,0,0,1,1},
//                            {0,0,1,0,1,0,0,0,1,0,1,1},
//                            {0,0,0,0,1,0,1,0,1,1,1,1},
//                            {1,1,0,0,0,0,1,1,1,1,0,0},
//                            {1,0,0,1,0,0,1,1,1,1,0,0}};
//        int[] queryGenotype =       {1,1,0,0,0,1,0,0,1,1,1,1};

        iteration++;
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

        // new solution
        IntSet[] solution = new IntSet[miniCost[0].length];
        for (int i = 0; i < solution.length; i++) {
            solution[i] = new IntOpenHashSet();
        }

        // initialize solution
        double currentHaplotypeMiniValue= Double.MAX_VALUE;
        IntSet currentMiniValueIndexSet = new IntOpenHashSet();
        for (int i = 0; i < miniCost.length; i++) {
            currentHaplotypeMiniValue = miniCost[i][queryGenotype.length-1] < currentHaplotypeMiniValue ? miniCost[i][queryGenotype.length-1] :currentHaplotypeMiniValue;
        }
        for (int i = 0; i < miniCost.length; i++) {
            if (currentHaplotypeMiniValue==miniCost[i][queryGenotype.length-1]){
                currentMiniValueIndexSet.add(i);
            }
        }
        solution[queryGenotype.length-1].addAll(currentMiniValueIndexSet);

        // find all solution

        IntIterator tIntIterator;
        int index;
        for (int i = miniCost[0].length-1; i > 0; i--) {
            tIntIterator = currentMiniValueIndexSet.iterator();
            while (tIntIterator.hasNext()){
                index = tIntIterator.nextInt();

                // 当前单倍型
                if (miniCost[index][i-1] <= miniCost[index][i]){
                    solution[i-1].add(index);
                }

                // 转换单倍型
                for (int k = 0; k < miniCost.length; k++) {
                    if (k==index) continue;
                    if ((miniCost[k][i-1]+switchCostScore) <= miniCost[index][i]){
                        solution[i-1].add(k);
                    }
                }
            }
            currentMiniValueIndexSet = solution[i-1];
        }
        EnumSet<WindowSource.Source>[] solutionSource= new EnumSet[solution.length];
        for (int i = 0; i < solutionSource.length; i++) {
            solutionSource[i] = EnumSet.noneOf(WindowSource.Source.class);
        }
        for (int i = 0; i < solutionSource.length; i++) {
            for (int ele:solution[i]){
                solutionSource[i].add(taxaSourceMap.get(srcIndiList.get(ele)));
            }
        }
        List<WindowSource.Source>[] solutionSourceList = new List[solutionSource.length];
        for (int i = 0; i < solutionSourceList.length; i++) {
            solutionSourceList[i]=new ArrayList<>(solutionSource[i]);
            Collections.sort(solutionSourceList[i]);
        }

        // transform solution to array
        List<WindowSource.Source[]> optiumSolutionList = new ArrayList<>();
        WindowSource.Source[] subSolution = new WindowSource.Source[queryGenotype.length];
        Arrays.fill(subSolution, null);
        optiumSolutionList.add(subSolution);
        int currentSolutionSize, multiplySolutionSize;
        currentSolutionSize = solutionSourceList[0].size();
        multiplySolutionSize = currentSolutionSize;
        for (int i = 0; i < currentSolutionSize-1; i++) {
            subSolution = new WindowSource.Source[queryGenotype.length];
            Arrays.fill(subSolution, null);
            optiumSolutionList.add(subSolution);
        }
        for (int i = 0; i < currentSolutionSize; i++) {
            optiumSolutionList.get(i)[0]=solutionSourceList[0].get(i);
        }

        for (int i = 1; i < solutionSourceList.length; i++) {
            currentSolutionSize = solutionSourceList[i].size();

            // {WE}
            if (currentSolutionSize == 1){
                for (int j = 0; j < optiumSolutionList.size(); j++) {
                    optiumSolutionList.get(j)[i] = solutionSourceList[i].get(0);
                }

                // {WE,DE}, {WE,DE}
            }else if (solution[i].equals(solution[i-1])){
                for (int j = 0; j < optiumSolutionList.size(); j++) {
                    optiumSolutionList.get(j)[i]=optiumSolutionList.get(j)[i-1];
                }
            }
//            else if (CollectionTool.hasIntersection(solution[i], solution[i-1])){
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
            else if (!solution[i].equals(solution[i-1])){
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
                    return GenotypeTable.getMiniPath2(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
                            taxaSourceMap, maxSolutionCount, maxSwitchCostScore);
                }

                // i-1 SNP 赋值
                for (int j = 0; j < multiplySolutionSize; j++) {
                    for (int k = 0; k < (currentSolutionSize-1); k++) {
                        System.arraycopy(optiumSolutionList.get(j),0,optiumSolutionList.get(multiplySolutionSize*(k+1)+j),
                                0, i);
                    }
                }

                // i SNP 赋值
                for (int j = 0; j < currentSolutionSize; j++) {
                    for (int k = 0; k < multiplySolutionSize; k++) {
                        optiumSolutionList.get(k+j*multiplySolutionSize)[i]=solutionSourceList[i].get(j);
                    }
                }
                multiplySolutionSize *=currentSolutionSize;
            }
        }
        System.out.println("iteration "+iteration);
        System.out.println("Switch cost score is "+switchCostScore);
        System.out.println();

        System.out.println("optium switch cost score is "+switchCostScore+" solution size is "+optiumSolutionList.size());
        iteration=0;
        return optiumSolutionList;
    }


}
