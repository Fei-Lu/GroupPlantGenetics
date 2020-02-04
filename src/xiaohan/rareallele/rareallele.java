package xiaohan.rareallele;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class rareallele {
    public rareallele(){
        this.rankGenes();
    }

    public void getup5kbSNP() throws IOException {
        //String SNPpositioninfile = "";
        String PhenotypeInfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/output/phenotypes.bed copy.txt";
        HashMap<String,Integer> GeneStartPointMap = new HashMap<String,Integer>();
        BufferedReader  br = IOUtils.getTextReader(PhenotypeInfile);
        String temp = null;
        String[] temps = null;
        String geneName = null;
        String geneStartPoint =null;
        Set<String> geneSet = new HashSet<String>();
        while((temp = br.readLine())!= null){
            if(temp.startsWith("chr")) continue;
            temps = temp.split("\t");
            geneName = temps[3];
            geneSet.add(geneName);
            geneStartPoint = temps[1];
            GeneStartPointMap.put(geneName,Integer.valueOf(geneStartPoint));
            System.out.println(GeneStartPointMap.get(geneName));




        }




    }
    public void rankGenes(){

    }

    public static void main (String[] args){
        new rareallele();
    }
}
