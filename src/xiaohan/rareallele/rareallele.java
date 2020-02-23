package xiaohan.rareallele;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.*;

public class rareallele {
    public rareallele(){
        //this.getExampleVCF();
        //this.rankGenes();
        //this.getTransNumber();
        this.getVCFposition();
    }

    public void getVCFposition(){
        String SampleName = "AT18488\n" +
                "AT18493\n" +
                "AT18494\n" +
                "AT18496\n" +
                "AT18497\n" +
                "AT18500\n" +
                "AT18502\n" +
                "AT18504\n" +
                "AT18509\n" +
                "AT18510\n" +
                "AT18514\n" +
                "AT18517\n" +
                "AT18520\n" +
                "AT18523\n" +
                "AT18528\n" +
                "AT18530\n" +
                "AT18531\n" +
                "AT18534\n" +
                "AT18535\n" +
                "AT18546\n" +
                "AT18549\n" +
                "AT18554\n" +
                "AT18556\n" +
                "AT18558\n" +
                "AT18564\n" +
                "AT18565\n" +
                "AT18570\n" +
                "AT18572\n" +
                "AT18580\n" +
                "AT18582\n" +
                "AT18597\n" +
                "AT18984\n" +
                "AT18606\n" +
                "AT18608\n" +
                "AT18615\n" +
                "AT18619\n" +
                "AT18620\n" +
                "AT18622\n" +
                "AT18632\n" +
                "AT18646\n" +
                "AT18650\n" +
                "AT18659\n" +
                "AT18663\n" +
                "AT18664\n" +
                "AT18667\n" +
                "AT18678\n" +
                "AT18681\n" +
                "AT18682\n" +
                "AT18683\n" +
                "AT18688\n" +
                "AT18692\n" +
                "AT18693\n" +
                "AT18696\n" +
                "SYR-L1\n" +
                "AT18837\n" +
                "AT18710\n" +
                "AT18713\n" +
                "AT18714\n" +
                "AT18719\n" +
                "AT18721\n" +
                "IRN-L2\n" +
                "AT18727\n" +
                "AT18728\n" +
                "AT18729\n" +
                "AT18732\n" +
                "TJK-L1\n" +
                "AT18741\n" +
                "AT18743\n" +
                "AT18746\n" +
                "AT18747\n" +
                "AT18749\n" +
                "AT18752\n" +
                "AFG-L1\n" +
                "AT18761\n" +
                "AT18838\n" +
                "AT18839\n" +
                "AT18765\n" +
                "AT18773\n" +
                "AT18774\n" +
                "AT18779\n" +
                "AT18780\n" +
                "AT18786\n" +
                "AT18790\n" +
                "AT18792\n" +
                "AT18794\n" +
                "AT18798\n" +
                "UZB-L1\n" +
                "AT18805\n" +
                "AT18808\n" +
                "AT18810\n" +
                "AT18819\n" +
                "AT18820\n" +
                "AT18842\n" +
                "AT18843\n" +
                "AT18828\n" +
                "AT18829\n";
        String[] Sample = null;
        Sample = SampleName.split("\n");
        Set<String> tempS = new HashSet<>();
        String infileS = "/Users/yxh/Documents/RareAllele/003rawdata/header.vcf";
        BufferedReader br =IOUtils.getTextReader(infileS);
        String temp =null;
        String[] temps = null;
        ArrayList<String> NameList = new ArrayList<>();
        try{

            while((temp = br.readLine())!=null) {
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].startsWith("AT")) {
                            temps[i] = temps[i].toString().substring(0,7);
                        }
                    }
                }
                for (int j = 0; j < Sample.length; j++) {
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].equals(Sample[j])) {
                            System.out.println(i);
                        }
                    }
                }
            }
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }


    public void getTransNumber(){
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/origin/genotypetrans.txt";
        String infor ="/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/origin/genotype1.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedReader br1 =IOUtils.getTextReader(infor);
        String temp1 =null;
        String temp =null;
        String[] temps =null;
        String[] temps1 = null;
        Set<String> Sample = new HashSet<String>();
        HashMap<String ,String > transMap = new HashMap<>();
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                transMap.put(temps[0], temps[1]);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        try{
            while((temp1 =br1.readLine()) != null) {
                temps1 = temp1.split("\t");
                String SampleName = temps1[1].toString();
                if (!Sample.contains(SampleName)) {
                    Sample.add(SampleName);
                }
            }
            String[] SampleS = Sample.toArray(new String[Sample.size()]);
            Arrays.sort(SampleS);
            for(int i = 0; i<SampleS.length;i++){
                System.out.print(transMap.get(SampleS[i])+"\t");
            }
        }
        catch (Exception e){
            e.printStackTrace();
        }

    }

    public void getExampleVCF(){
        String inforS ="/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/origin/genotypetrans.txt";
        BufferedReader br1 = IOUtils.getTextReader(inforS);
        String temp1 =null;
        String[] temps1 = null;
        HashMap<String,String> NameMap = new HashMap<>();
        try {
            while((temp1 = br1.readLine())!= null){
                temps1 = temp1.split("\t");
                String Name = temps1[0];
                String transName = temps1[1];
                NameMap.put(Name,transName);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String infileS = "/Users/yxh/Documents/RareAllele/003rawdata/header.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        String temp = null ;
        String[] temps = null;
        try{
            while((temp = br.readLine())!= null){
                if(temp.startsWith("##"))continue;
                if(temp.startsWith("#")){
                    temps = temp.split("\t");
                    for(int i =0;i<temps.length;i++){
                        if(NameMap.containsValue(temps[i])){
                            System.out.println(i+"\t"+temps[i]);
                        }
                    }
                }
                }


        }
        catch (Exception e){
            e.printStackTrace();
        }
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
