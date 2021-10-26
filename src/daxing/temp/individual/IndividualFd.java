package daxing.temp.individual;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import daxing.common.chrrange.ChrRange;
import daxing.common.chrrange.ChrRanges;
import daxing.common.utiles.IOTool;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

public class IndividualFd implements Comparable<IndividualFd> {

    ChrRanges chrRanges;
    Donor[] donors;
    double[] maxFdM;
    String taxon;

    static final double thresholdIntrogression=1;

    IntComparator intComparator = (int a, int b) -> chrRanges.getChrRange(a).compareTo(chrRanges.getChrRange(b));
    Swapper swapper = new Swapper() {
        @Override
        public void swap(int a, int b) {
            ChrRange chrRange= chrRanges.getChrRange(a);
            chrRanges.setChrRange(a, chrRanges.getChrRange(b));
            chrRanges.setChrRange(b, chrRange);
            Donor donor = donors[a];
            donors[a]=donors[b];
            donors[b]=donor;
            double fdM = maxFdM[a];
            maxFdM[a]=maxFdM[b];
            maxFdM[b]=fdM;
        }
    };

    public IndividualFd(String individualFdFile, String taxonName){
        List<ChrRange> chrRangeList= new ArrayList<>();
        List<Donor> donorList= new ArrayList<>();
        DoubleList maxFdMList=new DoubleArrayList();
        try (BufferedReader br = IOTool.getReader(individualFdFile)) {
            String line, refChr;
            List<String> temp;
            br.readLine();
            int start,end;
            ChrRange chrRange;
            double maxFdM;
            Donor donor;
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                refChr = temp.get(0);
                start=Integer.parseInt(temp.get(1));
                end=Integer.parseInt(temp.get(2));
                chrRange=new ChrRange(refChr, start, end);
                donor= temp.get(3).equals("NA") ? Donor.NONE : Donor.valueOf(temp.get(3));
                maxFdM = Double.parseDouble(temp.get(4));
                donor = maxFdM < thresholdIntrogression ? Donor.NONE : donor;
                maxFdMList.add(maxFdM);
                chrRangeList.add(chrRange);
                donorList.add(donor);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        Donor[] donors=new Donor[donorList.size()];
        for (int i = 0; i < donors.length; i++) {
            donors[i]=donorList.get(i);
        }
        this.chrRanges= new ChrRanges(chrRangeList);
        this.donors=donors;
        this.maxFdM=maxFdMList.toDoubleArray();
        this.taxon=taxonName;
        GenericSorting.quickSort(0, this.chrRanges.getRangeNumber(), intComparator, swapper);
        this.chrRanges.setSortType(ChrRanges.SortType.POSITION);
    }

    public IndividualFd(String taxon){
        this.taxon=taxon;
        this.chrRanges=null;
        this.maxFdM=null;
        this.donors=null;
    }

    public Donor[] getDonors() {
        return donors;
    }

    private int[] getIndicesContainPosition(String chr, int pos){
        return this.chrRanges.getIndicesContainsPosition(chr, pos);
    }

    public EnumSet<Donor> getDonorFrom(String chr, int pos){
        EnumSet<Donor> donorEnumSet= EnumSet.noneOf(Donor.class);
        int[] indices=getIndicesContainPosition(chr, pos);
        for (int index : indices){
            donorEnumSet.add(this.getDonors()[index]);
        }
        if (donorEnumSet.size() <=1) return donorEnumSet;
        donorEnumSet.remove(Donor.NONE);
        return donorEnumSet;
    }

    @Override
    public int compareTo(IndividualFd o) {
        return this.taxon.compareTo(o.taxon);
    }
}
