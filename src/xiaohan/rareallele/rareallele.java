package xiaohan.rareallele;

import smile.stat.Stat;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

public class rareallele {
    public rareallele(){
        //this.getExampleVCF();
        //this.rankGenes();
        //this.getTransNumber();
        //this.getVCFposition();
        this.getsubVCF();
    }

    public void getsubVCF(){
        String infileS = "/data2/junxu/geneSNP/all.vcf.gz";
        String outputDir = "/data1/home/xiaohan/rareallele/subvcf";
        String index ="12\n" +
                "17\n" +
                "18\n" +
                "20\n" +
                "21\n" +
                "24\n" +
                "26\n" +
                "28\n" +
                "33\n" +
                "34\n" +
                "38\n" +
                "41\n" +
                "44\n" +
                "47\n" +
                "52\n" +
                "54\n" +
                "55\n" +
                "58\n" +
                "59\n" +
                "70\n" +
                "73\n" +
                "78\n" +
                "80\n" +
                "82\n" +
                "88\n" +
                "89\n" +
                "94\n" +
                "96\n" +
                "104\n" +
                "106\n" +
                "119\n" +
                "335\n" +
                "126\n" +
                "128\n" +
                "134\n" +
                "138\n" +
                "139\n" +
                "141\n" +
                "149\n" +
                "163\n" +
                "166\n" +
                "173\n" +
                "177\n" +
                "178\n" +
                "181\n" +
                "192\n" +
                "195\n" +
                "196\n" +
                "197\n" +
                "201\n" +
                "205\n" +
                "206\n" +
                "209\n" +
                "351\n" +
                "328\n" +
                "221\n" +
                "224\n" +
                "225\n" +
                "229\n" +
                "231\n" +
                "342\n" +
                "236\n" +
                "237\n" +
                "238\n" +
                "241\n" +
                "352\n" +
                "249\n" +
                "251\n" +
                "254\n" +
                "255\n" +
                "257\n" +
                "259\n" +
                "9\n" +
                "266\n" +
                "329\n" +
                "330\n" +
                "268\n" +
                "276\n" +
                "277\n" +
                "281\n" +
                "282\n" +
                "288\n" +
                "292\n" +
                "294\n" +
                "296\n" +
                "300\n" +
                "357\n" +
                "306\n" +
                "309\n" +
                "311\n" +
                "317\n" +
                "318\n" +
                "331\n" +
                "332\n" +
                "322\n" +
                "323\n";
        String[] indexes = null;
        indexes = index.split("\n");
        BufferedReader br = IOUtils.getTextGzipReader(infileS);
        String temp = null;
        String[] temps = null;
        BufferedWriter bw = IOUtils.getTextGzipWriter(new File(outputDir,"colVCF.vcf.gz").getAbsolutePath());
        try{
            while ((temp = br.readLine()) != null){
                if(temp.startsWith("##")){
                    bw.write(temp);bw.newLine();
                    continue;
                }
                else {
                    temps = temp.split("\t");
                    bw.write(temps[0] + "\t" + temps[1] + "\t" + temps[2] + "\t" + temps[3] + "\t" + temps[4] + "\t" + temps[5] + "\t" +
                            temps[6] + "\t" + temps[7] + "\t" + temps[8] + "\t");
                    for (int i = 0; i < indexes.length; i++) {
                        bw.write(temps[Integer.parseInt(indexes[i])]+"\t");
                    }
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }

        catch (Exception e){
            e.printStackTrace();
        }

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
