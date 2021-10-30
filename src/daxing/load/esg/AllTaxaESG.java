package daxing.load.esg;

import daxing.common.utiles.IOTool;
import daxing.common.factors.Ploidy;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang.ArrayUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AllTaxaESG {

    String[] taxa;
    String[] traidBlockID;
    // -1: NA
    int[][][] cells; // dim1: triads, dim2: taxa, dim3: subgenome

    public AllTaxaESG(String[] taxa, String[] traidBlockID, int[][][] cells){
        this.taxa =taxa;
        this.traidBlockID = traidBlockID;
        this.cells=cells;
    }

    public static AllTaxaESG getInstance(String allESGFile){
        long start=System.nanoTime();
        String[] taxa=null, triadsBlock=null;
        int[][][] cells=null;
        try (BufferedReader br = IOTool.getReader(allESGFile)) {
            String line;
            List<String> temp, tem;
            line = br.readLine();
            temp = PStringUtils.fastSplit(line);
            taxa = temp.subList(1, temp.size()).stream().toArray(String[]::new);
            List<String> triadsBlockIDList= new ArrayList<>();
            int[][] esgTaxa;
            List<int[][]> esgTaxaTriadsList = new ArrayList<>();
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                triadsBlockIDList.add(temp.get(0));
                esgTaxa= new int[taxa.length][Ploidy.HEXAPLOID.getSubgenomewNum()];
                for (int i = 0; i < esgTaxa.length; i++) {
                    esgTaxa[i] = new int[Ploidy.HEXAPLOID.getSubgenomewNum()];
                    Arrays.fill(esgTaxa[i], -1);
                }
                for (int i = 1; i < temp.size(); i++) {
                    if (temp.get(i).equals("NA")) continue;
                    tem = PStringUtils.fastSplit(temp.get(i),";");
                    for (int j = 0; j < tem.size(); j++) {
                        esgTaxa[i-1][j] = Integer.parseInt(tem.get(j));
                    }
                }
                esgTaxaTriadsList.add(esgTaxa);
            }
            triadsBlock= triadsBlockIDList.stream().toArray(String[]::new);
            cells = new int[esgTaxaTriadsList.size()][taxa.length][Ploidy.HEXAPLOID.getSubgenomewNum()];
            for (int i = 0; i < esgTaxaTriadsList.size(); i++) {
                cells[i] = esgTaxaTriadsList.get(i);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("new AllTaxaESG instance spend "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
        return new AllTaxaESG(taxa, triadsBlock, cells);
    }

    public int[][][] getCells() {
        return cells;
    }

    public String[] getTaxa() {
        return taxa;
    }

    public String[] getTraidBlockID() {
        return traidBlockID;
    }

    public int getTriadsNum(){
        return this.getTraidBlockID().length;
    }

    public int getTaxaNum(){
        return this.getTaxa().length;
    }

    public TIntArrayList getTaxaIndex(IfPseudoHexaploid ifPseudoHexaploid){
        String[] taxa=this.getTaxa();
        TIntArrayList indexes=new TIntArrayList();
        for (int i = 0; i < taxa.length; i++) {
            switch (ifPseudoHexaploid){
                case HEXAPLOID:
                    if (taxa[i].contains("_D_")) continue;
                    indexes.add(i);
                    break;
                case PSEUDO_HEXAPLOID:
                    if (!taxa[i].contains("_D_")) continue;
                    indexes.add(i);
                    break;
            }
        }
        return indexes;
    }

    public int[] getCell(String triadsBlockID, String taxon){
        int triadsBlockIDIndex = ArrayUtils.indexOf(this.getTraidBlockID(), triadsBlockID);
        int taxonIndex = ArrayUtils.indexOf(this.getTaxa(), taxon);
        return getCell(triadsBlockIDIndex, taxonIndex);
    }

    public int[] getCell(int triadsBlockIDIndex, int taxonIndex){
        return this.getCells()[triadsBlockIDIndex][taxonIndex];
    }

    public List<String> getNonNATaxaList(int triadsBlockIDIndex, IfPseudoHexaploid ifPseudoHexaploid){
        List<String> nonNATaxaList= new ArrayList<>();
        TIntArrayList taxaIndexes=this.getTaxaIndex(ifPseudoHexaploid);
        int[][] cells = this.getCells()[triadsBlockIDIndex];
        String[] taxa=this.getTaxa();
        for (int i = 0; i < taxaIndexes.size(); i++) {
            for (int j = 0; j < cells[taxaIndexes.get(i)].length; j++) {
                if (cells[taxaIndexes.get(i)][j] < 0) continue;
                nonNATaxaList.add(taxa[taxaIndexes.get(i)]);
                break;
            }
        }
        return nonNATaxaList;
    }

    public int getNonNATaxaNum(int triadsBlockIDIndex, IfPseudoHexaploid ifPseudoHexaploid){
        return this.getNonNATaxaList(triadsBlockIDIndex, ifPseudoHexaploid).size();
    }

    public int[] getESGObservedCount(int triadsBlockIDIndex, IfPseudoHexaploid ifPseudoHexaploid){
        int[] esgObservedCount=new int[ESG.values().length];
        int[][] cells= this.getCells()[triadsBlockIDIndex];
        TIntArrayList taxaIndexes=this.getTaxaIndex(ifPseudoHexaploid);
        String esgInfo;
        StringBuilder sb=new StringBuilder();
        ESG esg;
        for (int i = 0; i < taxaIndexes.size(); i++) {
            if (!(ArrayUtils.indexOf(cells[taxaIndexes.get(i)], -1) < 0)) continue;
            sb.setLength(0);
            for (int j = 0; j < cells[taxaIndexes.get(i)].length; j++) {
                assert cells[taxaIndexes.get(i)][j] > -1 : " exist -1";
                sb.append(cells[taxaIndexes.get(i)][j]);
            }
            esgInfo = sb.toString();
            esg = ESG.newInstanceFrom(esgInfo);
            esgObservedCount[esg.getValue()]++;
        }
        return esgObservedCount;
    }

    public String getTriadsBlockID(int triadsBlockIDIndex){
        return this.getTraidBlockID()[triadsBlockIDIndex];
    }
}
