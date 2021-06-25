package daxing.hybrid.dection;

import daxing.common.IOTool;
import gnu.trove.list.array.TIntArrayList;

import java.io.BufferedWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.List;

public class IndividualIntervalHeterozygosity {

    List<String> taxonList;
    SlidingWindow[][] taxonHomoHeteroCount; //dim1: taxon, dim2: homo or hetero

    IndividualIntervalHeterozygosity(List<String> taxonList, int chrLen, int windowSize, int windowStep){
        this.taxonList=taxonList;
        taxonHomoHeteroCount=new SlidingWindow[taxonList.size()][];
        for (int i = 0; i < taxonHomoHeteroCount.length; i++) {
            taxonHomoHeteroCount[i]=new SlidingWindow[2];
            for (int j = 0; j < taxonHomoHeteroCount[i].length; j++) {
                taxonHomoHeteroCount[i][j]=new SlidingWindow(chrLen, windowSize, windowStep);
            }
        }
    }

    public List<String> getTaxonList() {
        return taxonList;
    }

    public SlidingWindow[][] getTaxonHomoHeteroCount() {
        return taxonHomoHeteroCount;
    }

    public int getWindowNum(){
        return this.getTaxonHomoHeteroCount()[0][0].starts.length;
    }

    public int[] getWindowStart(){
        return this.getTaxonHomoHeteroCount()[0][0].getWindowStarts();
    }

    public int[] getWindowEnd(){
        return this.getTaxonHomoHeteroCount()[0][0].getWindowEnds();
    }

    public int getTaxonIndex(String taxon){
        return taxonList.indexOf(taxon);
    }

    public void addTaxonHomoHeteroCount(List<IndividualHomoHeteroPosition> taxaHomoHeteroPositionList){
        String taxon;
        TIntArrayList homoSNPPosList, heteroSNPPosList;
        int taxonIndex;
        for (IndividualHomoHeteroPosition individualHomoHeteroPosition: taxaHomoHeteroPositionList){
            taxon=individualHomoHeteroPosition.getTaxon();
            homoSNPPosList=individualHomoHeteroPosition.getHomoSNPPositionList();
            heteroSNPPosList=individualHomoHeteroPosition.getHeteroSNPPositionList();
            taxonIndex=this.getTaxonIndex(taxon);
            if (taxonIndex < 0){
                System.out.println("check your taxon list, it had duplicated taxon name, program will quit");
                System.exit(1);
            }
            this.taxonHomoHeteroCount[taxonIndex][0].addPositionCount(homoSNPPosList.toArray());
            this.taxonHomoHeteroCount[taxonIndex][1].addPositionCount(heteroSNPPosList.toArray());
        }
    }

    public void write(String outFile, String chr){
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb = new StringBuilder();
            List<String> taxonList=this.getTaxonList();
            sb.append("Chr").append("\t").append("Start").append("\t").append("End").append("\t");
            sb.append(String.join("\t", taxonList));
            bw.write(sb.toString());
            bw.newLine();
            double individualIntervalHeterozygosity;
            double[] homoWindowValueArray, heteroWindowValueArray;
            int[] windowStart=this.getWindowStart();
            int[] windowEnd= this.getWindowEnd();
            NumberFormat numberFormat = NumberFormat.getInstance();
            numberFormat.setGroupingUsed(false);
            numberFormat.setMaximumFractionDigits(5);
            for (int i = 0; i < this.getWindowNum(); i++) {
                sb.setLength(0);
                sb.append(chr).append("\t").append(windowStart[i]).append("\t");
                sb.append(windowEnd[i]).append("\t");
                for (int j = 0; j < taxonList.size(); j++) {
                    homoWindowValueArray= this.getTaxonHomoHeteroCount()[j][0].getWindowValuesDouble();
                    heteroWindowValueArray= this.getTaxonHomoHeteroCount()[j][1].getWindowValuesDouble();
                    individualIntervalHeterozygosity = heteroWindowValueArray[i]/(heteroWindowValueArray[i]+homoWindowValueArray[i]);
                    if (Double.isNaN(individualIntervalHeterozygosity)){
                        sb.append("NA").append("\t");
                    }else {
                        sb.append(numberFormat.format(individualIntervalHeterozygosity)).append("\t");
                    }
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static class IndividualHomoHeteroPosition{

        String taxon;
        // position: 1-based
        TIntArrayList homoSNPPositionList;
        TIntArrayList heteroSNPPositionList;

        IndividualHomoHeteroPosition(String taxon, TIntArrayList homoSNPPositionList, TIntArrayList heteroSNPPositionList){
            this.taxon=taxon;
            this.homoSNPPositionList= homoSNPPositionList;
            this.heteroSNPPositionList= heteroSNPPositionList;
        }

        IndividualHomoHeteroPosition(String taxon){
            this.taxon=taxon;
            this.homoSNPPositionList= new TIntArrayList();
            this.heteroSNPPositionList= new TIntArrayList();
        }

        public String getTaxon() {
            return taxon;
        }

        public TIntArrayList getHeteroSNPPositionList() {
            return heteroSNPPositionList;
        }

        public TIntArrayList getHomoSNPPositionList() {
            return homoSNPPositionList;
        }

        /**
         *
         * @param position 1-based
         */
        public void addHomoPos(int position){
            this.homoSNPPositionList.add(position);
        }

        /**
         *
         * @param position 1-based
         */
        public void addHeteroPos(int position){
            this.heteroSNPPositionList.add(position);
        }
    }



}
