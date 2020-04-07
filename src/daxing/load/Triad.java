package daxing.load;

import daxing.common.IOTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Triad {

    List<String> triadID;
    List<String[]> triad;
    TIntArrayList ifSyntenic;  //0 is non-syntenic and 1 is syntenic
    TIntArrayList ifExpressed;  // 0 is false and 1 is true

    public Triad(String triadFile){
        triadID=new ArrayList<>(19000);
        triad=new ArrayList<>(19000);
        ifSyntenic=new TIntArrayList();
        ifExpressed=new TIntArrayList();
        try (BufferedReader br = IOTool.getReader(triadFile)) {
            String line;
            List<String> temp;
            br.readLine();
            String[] abdGenes;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                triadID.add(temp.get(0));
                abdGenes=new String[3];
                abdGenes[0]=temp.get(1);
                abdGenes[1]=temp.get(2);
                abdGenes[2]=temp.get(3);
                triad.add(abdGenes);
                if (temp.get(4).equals("syntenic")){
                    ifSyntenic.add(1);
                }else {
                    ifSyntenic.add(0);
                }
                if (temp.get(5).equals("TRUE")){
                    ifExpressed.add(1);
                }else {
                    ifExpressed.add(0);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public List<String[]> getTriad() {
        return triad;
    }

    public TIntArrayList getIfExpressed() {
        return ifExpressed;
    }

    public TIntArrayList getIfSyntenic() {
        return ifSyntenic;
    }

    public List<String> getTriadID() {
        return triadID;
    }

    public WheatLineage getSubgenome(String geneName){
        List<String> temp=PStringUtils.fastSplit(geneName, "G");
        String subgenome=temp.get(0).substring(8,9);
        return WheatLineage.valueOf(subgenome);
    }

    public List<String> getSubgenomeGenes(WheatLineage subgenome){
        String[] subgenomes={"A","B","D"};
        int index= Arrays.binarySearch(subgenomes, subgenome.name());
        if (index < 0) {
            System.out.println("subgenome must be A or B or D");
            System.exit(1);
        }
        List<String[]> genes=this.getTriad();
        List<String> subgenomeGenes=new ArrayList<>();
        for (int i = 0; i < genes.size(); i++) {
            subgenomeGenes.add(genes.get(i)[index]);
        }
        return subgenomeGenes;
    }

    public List<String> getSubgenomeGens(String geneName){
        WheatLineage subgenome=this.getSubgenome(geneName);
        return this.getSubgenomeGenes(subgenome);
    }

    public List<String> getAllGenes(){
        List<String> genes=new ArrayList<>(57000);
        List<String[]> data=this.getTriad();
        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data.get(i).length; j++) {
                genes.add(data.get(i)[j]);
            }
        }
        return genes;
    }

    public String[] getTriad(int geneIndex){
        return this.getTriad().get(geneIndex);
    }

    public String[] getTriad(String geneName){
        int geneIndex=this.getGeneIndex(geneName);
        return getTriad(geneIndex);
    }

    public String[] getTraidGenes(String triadID){
        List<String> triadNames=this.triadID;
        int index=triadNames.indexOf(triadID);
        return this.getTriad().get(index);
    }

    public String getTraidID(String geneName){
        int geneIndex=this.getGeneIndex(geneName);
        return this.getTriadID().get(geneIndex);
    }

    public int getGeneIndex(String geneName){
        String subgenome=this.getSubgenome(geneName).name();
        List<String> subgenomeGenes=this.getSubgenomeGenes(WheatLineage.valueOf(subgenome));
        return subgenomeGenes.indexOf(geneName);
    }

    public boolean contain(String geneName){
        String subgenome=this.getSubgenome(geneName).name();
        List<String> subgenomeGenes=this.getSubgenomeGenes(WheatLineage.valueOf(subgenome));
        return subgenomeGenes.contains(geneName);
    }

}
