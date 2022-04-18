package daxing.individualIntrogression.pop2Indi;

import daxing.common.chrrange.ChrRange;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import daxing.individualIntrogression.individual.Donor;
import gnu.trove.list.array.TIntArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import java.io.File;
import java.util.*;


/**
 * this class used to storage all individuals fd(pop2Indi), dxy(pop2Indi), Donor, Timing
 * ChrRange Window: 200 SNP, Step: 100 SNP, be storage in array, can not add or delete ChrRange
 * Taxon info such as sitesUsed:RelativeTiming_AC is storage in List, can add or delete taxon
 */
public class IndividualsFdDxyDonor {

    // dim1: ChrRange
    ChrRange[] chrRanges;
    int[] mid;
    int[] sites;

    // dim1: taxon
    List<String> taxa;
    Map<String, Integer> taxonIndexMap;

    // list: taxon, array: ChrRange
    IntList[] siteUsed;
    DoubleList[] abba;
    DoubleList[] baba;
    DoubleList[] d;
    DoubleList[] fd;
    DoubleList[] fdM;
    DoubleList[] fd_mean;
    IntList[] ifIntrogressionByfd;
    List<Donor>[] donorPredictedByMindxy;
    List<Donor>[] donorPredictedByMindxy_fd;
    DoubleList[] rndMin;
    DoubleList[] relativeTiming_AC;

    public IndividualsFdDxyDonor(String individualFdDxyDonorDir){
        List<File> files = IOTool.getFileListInDirEndsWith(individualFdDxyDonorDir, ".txt.gz");
        RowTableTool<String> table = new RowTableTool<>(files.get(0).getAbsolutePath());
        this.initializeChrRange(table);
        int chrRangNum = this.getChrRangeNum();
        this.siteUsed =  new IntList[chrRangNum];
        this.abba = new DoubleList[chrRangNum];
        this.baba = new DoubleList[chrRangNum];
        this.d = new DoubleList[chrRangNum];
        this.fd = new DoubleList[chrRangNum];
        this.fdM = new DoubleList[chrRangNum];
        this.fd_mean = new DoubleList[chrRangNum];
        this.ifIntrogressionByfd = new IntList[chrRangNum];
        this.donorPredictedByMindxy = new ArrayList[chrRangNum];
        this.donorPredictedByMindxy_fd = new ArrayList[chrRangNum];
        this.rndMin = new DoubleList[chrRangNum];
        this.relativeTiming_AC = new DoubleList[chrRangNum];
        for (int i = 0; i < chrRangNum; i++) {
            this.siteUsed[i] = new IntArrayList();
            this.abba[i] = new DoubleArrayList();
            this.baba[i] = new DoubleArrayList();
            this.d[i] = new DoubleArrayList();
            this.fd[i] = new DoubleArrayList();
            this.fdM[i] = new DoubleArrayList();
            this.fd_mean[i] = new DoubleArrayList();
            this.ifIntrogressionByfd[i] = new IntArrayList();
            this.donorPredictedByMindxy[i] = new ArrayList<>();
            this.donorPredictedByMindxy_fd[i] = new ArrayList<>();
            this.rndMin[i] = new DoubleArrayList();
            this.relativeTiming_AC[i] = new DoubleArrayList();
        }
        List<RowTableTool<String>> tables = new ArrayList<>();
        this.taxa = new ArrayList<>();
        this.taxonIndexMap = new HashMap<>();
        String taxon;
        int index = 0;
        for (File file : files){
            tables.add(new RowTableTool<>(file.getAbsolutePath()));
            taxon = file.getName().substring(13).replaceAll("_vmap2.1.txt.gz","");
            this.taxa.add(taxon);
            this.taxonIndexMap.put(taxon, index);
            index++;
        }
        this.initializeTaxa(tables);
    }

    /**
     *
     * @param chrRangesTable The first five columns must contain: scaffold, start, end, mid, sites
     */
    private void initializeChrRange(RowTableTool<String> chrRangesTable){
       List<String> chr = chrRangesTable.getColumn(0);
       int[] start = chrRangesTable.getColumnAsIntArray("start");
       int[] end = chrRangesTable.getColumnAsIntArray("end");
       int[] mid = chrRangesTable.getColumnAsIntArray("mid");
       int[] sites = chrRangesTable.getColumnAsIntArray("sites");
       ChrRange[] chrRanges = new ChrRange[chr.size()];
        for (int i = 0; i < start.length; i++) {
            chrRanges[i] = new ChrRange(chr.get(i), start[i], end[i]);
        }
        Arrays.sort(chrRanges);
        this.chrRanges=chrRanges;
        this.mid=mid;
        this.sites=sites;
    }

    /**
     *
     * @param individualTables see follow
     * scaffold	start	end	mid	sites	sitesUsed	ABBA	BABA	D	fd	fdM	fd_mean	IfIntrogressionByfd	DonorPredictedByMindxy	RNDMin	RelativeTiming_AC
     */
    private void initializeTaxa(List<RowTableTool<String>> individualTables){
        List<String> sitesUsed, abba, baba, d, fd, fdM, fd_mean, ifIntrogressionByfd, rndMin, relativeTiming_AC;
        List<String> donorPredictedByMindxy, donorPredictedByMindxy_fd;
        int chrRangeNum = this.getChrRangeNum();
        List<String> headList = this.getHeader();
        Comparator<List<String>> comparatorPos= Comparator.comparing(l->l.get(0));
        comparatorPos = comparatorPos.thenComparing(l->Integer.parseInt(l.get(1)));
        comparatorPos = comparatorPos.thenComparing(l->Integer.parseInt(l.get(2)));
        Donor donor;
        int ifIntrogressionByFd;
        for (RowTableTool<String> tableTool : individualTables){
            assert tableTool.getRowNumber() == chrRangeNum : "individualFdDxyDonor must contain same num row";
            assert tableTool.getHeader().equals(headList) : "Column names must meet method requirements";
            tableTool.sortBy(comparatorPos);
            sitesUsed = tableTool.getColumn(5);
            abba = tableTool.getColumn(6);
            baba = tableTool.getColumn(7);
            d = tableTool.getColumn(8);
            fd = tableTool.getColumn(9);
            fdM = tableTool.getColumn(10);
            fd_mean = tableTool.getColumn(11);
            ifIntrogressionByfd = tableTool.getColumn(12);
            donorPredictedByMindxy = tableTool.getColumn(13);
            rndMin = tableTool.getColumn(14);
            relativeTiming_AC = tableTool.getColumn(15);
            for (int i = 0; i < this.getChrRangeNum(); i++) {
                this.siteUsed[i].add(sitesUsed.get(i).equals("NA") ? -1 : Integer.parseInt(sitesUsed.get(i)));
                this.abba[i].add(abba.get(i).equals("NA") ? -1 : Double.parseDouble(abba.get(i)));
                this.baba[i].add(baba.get(i).equals("NA") ? -1 : Double.parseDouble(baba.get(i)));
                this.d[i].add(d.get(i).equals("NA") ? -1 : Double.parseDouble(d.get(i)));
                this.fd[i].add(fd.get(i).equals("NA") ? -1 : Double.parseDouble(fd.get(i)));
                this.fdM[i].add(fdM.get(i).equals("NA") ? -1 : Double.parseDouble(fdM.get(i)));
                this.fd_mean[i].add(fd_mean.get(i).equals("NA") ? -1 : Double.parseDouble(fd_mean.get(i)));
                ifIntrogressionByFd = ifIntrogressionByfd.get(i).equals("NA") ? -1 : Integer.parseInt(ifIntrogressionByfd.get(i));
                this.ifIntrogressionByfd[i].add(ifIntrogressionByFd);
                donor = donorPredictedByMindxy.get(i).equals("NA") ? Donor.NONE : Donor.valueOf(donorPredictedByMindxy.get(i));
                this.donorPredictedByMindxy[i].add(donor);
                this.donorPredictedByMindxy_fd[i].add(ifIntrogressionByFd==1 ? donor : Donor.NONE);
                this.rndMin[i].add(rndMin.get(i).equals("NA") ? -1 : Double.parseDouble(rndMin.get(i)));
                this.relativeTiming_AC[i].add(relativeTiming_AC.get(i).equals("NA") ? -1 : Double.parseDouble(relativeTiming_AC.get(i)));
            }
        }
    }

    private List<String> getHeader(){
        String[] head={"scaffold","start","end","mid","sites","sitesUsed","ABBA","BABA","D","fd","fdM","fd_mean",
                "IfIntrogressionByfd","DonorPredictedByMindxy","RNDMin","RelativeTiming_AC"};
        List<String> headList = new ArrayList<>();
        for (String str:head){
            headList.add(str);
        }
        return headList;
    }

    public int getChrRangeNum(){
        return this.chrRanges.length;
    }

    public int getTaxaNum(){
        return this.getTaxa().size();
    }

    public List<String> getTaxa() {
        return taxa;
    }

    public ChrRange[] getChrRanges() {
        return chrRanges;
    }

    public List<Donor>[] getDonorPredictedByMindxy() {
        return donorPredictedByMindxy;
    }

    public IntList[] getIfIntrogressionByfd() {
        return ifIntrogressionByfd;
    }

    public List<Donor>[] getDonorPredictedByMindxy_fd() {
        return donorPredictedByMindxy_fd;
    }

    public Map<String, Integer> getTaxonIndexMap() {
        return taxonIndexMap;
    }

    public int getTaxonIndex(String taxon){
        return this.getTaxonIndexMap().get(taxon);
    }

    /**
     *
     * @param refChr ref chr
     * @param refPos ref pos
     * @return new int[0] if do not contain
     */
    public int[] getChrRangeIndicesContainsPosition (String refChr, int refPos) {
        TIntArrayList indexList = new TIntArrayList();
        ChrRange query = new ChrRange(refChr, refPos, refPos+1);
        int hit = Arrays.binarySearch(this.getChrRanges(), query);
        hit = hit < -1 ? -hit-2 : hit;
        if (hit < 0) return new int[0];
        for (int i = hit; i > -1; i--) {
            if (this.getChrRanges()[i].contain(refChr, refPos)) indexList.add(i);
        }
        indexList.sort();
        return indexList.toArray();
    }

    public boolean containTaxon(String taxon){
        List<String> taxaList = this.getTaxa();
        return taxaList.contains(taxon);
    }


    /**
     * if no introgression, return blank set
     * @param taxon
     * @param chr
     * @param pos
     * @return
     */
    public EnumSet<Donor> getDonorFrom(String taxon, String chr, int pos){
        EnumSet<Donor> donorEnumSet= EnumSet.noneOf(Donor.class);
        if (!this.containTaxon(taxon)) return donorEnumSet;
        int taxonIndex = this.getTaxonIndex(taxon);
        int[] indices=this.getChrRangeIndicesContainsPosition(chr, pos);
        for (int index : indices){
            donorEnumSet.add(this.getDonorPredictedByMindxy_fd()[index].get(taxonIndex));
        }
        if (donorEnumSet.size() <=1) return donorEnumSet;
        donorEnumSet.remove(Donor.NONE);
        return donorEnumSet;
    }
}
