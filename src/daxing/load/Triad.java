package daxing.load;

import daxing.common.IOTool;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Triad {

    List<String>[] triad;
    TIntArrayList ifSyntenic;  //0 is non-syntenic and 1 is syntenic
    TIntArrayList ifExpressed;  // 0 is false and 1 is true

    public Triad(String triadFile){
        triad=new List[3];
        ifSyntenic=new TIntArrayList();
        ifExpressed=new TIntArrayList();
        for (int i = 0; i < triad.length; i++) {
            triad[i]=new ArrayList<>(19000);
        }
        try (BufferedReader br = IOTool.getReader(triadFile)) {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                triad[0].add(temp.get(0).replaceFirst("01", "02"));
                triad[1].add(temp.get(1).replaceFirst("01","02"));
                triad[2].add(temp.get(2).replaceFirst("01", "02"));
                if (temp.get(3).equals("syntenic")){
                    ifSyntenic.add(1);
                }else {
                    ifSyntenic.add(0);
                }
                if (temp.get(4).equals("TRUE")){
                    ifExpressed.add(1);
                }else {
                    ifExpressed.add(0);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public Triad(List<String> aTriad, List<String> bTriad, List<String> dTriad){
        triad=new List[3];
        for (int i = 0; i < triad.length; i++) {
            triad[i]=new ArrayList<>();
        }
        triad[0]=aTriad;
        triad[1]=bTriad;
        triad[2]=dTriad;
    }

    public List<String>[] getTriad() {
        return triad;
    }

    public TIntArrayList getIfExpressed() {
        return ifExpressed;
    }

    public TIntArrayList getIfSyntenic() {
        return ifSyntenic;
    }

    public List<String> getSubgenomeGenes(String subgenome){
        String[] subgenomes={"A","B","D"};
        int index= Arrays.binarySearch(subgenomes, subgenome);
        if (index < 0) {
            System.out.println("subgenome must be A or B or D");
            System.exit(1);
        }
        return this.getTriad()[index];
    }

    public String[] getTriad(int geneIndex){
        List<String>[] triad=this.getTriad();
        String[] res=new String[3];
        res[0]=triad[0].get(geneIndex);
        res[1]=triad[1].get(geneIndex);
        res[2]=triad[2].get(geneIndex);
        return res;
    }

    public String[] getTriad(String geneName){
        int geneIndex=this.getGeneIndex(geneName);
        return getTriad(geneIndex);
    }

    public int getGeneIndex(String geneName){
        String[] subgenomes={"A","B","D"};
        List<String> temp=PStringUtils.fastSplit(geneName, "G");
        String subgenome=temp.get(0).substring(8,9);
        int subgenomeIndex=Arrays.binarySearch(subgenomes, subgenome);
        int geneIndex=this.getTriad()[subgenomeIndex].indexOf(geneName);
        return geneIndex;
    }

    public boolean contain(String geneName){
        String[] subgenomes={"A","B","D"};
        List<String> temp=PStringUtils.fastSplit(geneName, "G");
        String subgenome=temp.get(0).substring(8,9);
        int index=Arrays.binarySearch(subgenomes, subgenome);
        List<String> subgenomeTriad=this.triad[index];
        return subgenomeTriad.contains(geneName);
    }

}
