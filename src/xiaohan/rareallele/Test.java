package xiaohan.rareallele;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class Test {

    public Test() throws IOException {
        //this.findTaxon();
        //this.findSNPnumber();
        this.findFalseSample();
    }

    public void findFalseSample(){
        String infileS = "/Users/yxh/Documents/eQTL/003experiment/浓度汇总/RNAconcentration.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        String[] temps = null ;
        String temp = null ;
        HashSet<String> NameSet = new HashSet<String>() ;
        HashMap<String, Double> A280Map = new HashMap<>();
        HashMap<String, Double> A230Map = new HashMap<>();
        HashMap<String, Double> ConcentrationMap = new HashMap<>();
        try{
            while((temp = br.readLine()) != null){
                if(temp.startsWith("#")){
                    continue;
                }
                temps = temp.split("\t");
                String Name = temps[0];
                Double A280 = Double.valueOf(temps[1]);
                Double A230 = Double.valueOf(temps[2]);
                Double Concentration = Double.valueOf(temps[3]);
                NameSet.add(Name);
                A280Map.put(Name,A280);
                A230Map.put(Name,A230);
                ConcentrationMap.put(Name,Concentration);
                //A280不合格，A230合格
                /*if((A280 > 2.3 || A280 < 1.6) && (A230<2.5 && A230 > 1.6)){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }
                //A280合格，A230不合格
                if((A280 < 2.3 && A280 > 1.6) && (A230 > 2.5 || A230 < 1.6)){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }*/
                //A280，A230都不合格
                /*if((A280 > 2.3 || A280 < 1.6) && (A230 > 2.5 || A230 < 1.6)){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }*/
                //浓度不合格
                /*if(Concentration < 30){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }*/
                if((A280 > 2.3 || A280 < 1.6) || (A230 > 2.5 || A230 < 1.6) || (Concentration < 30)){
                    System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }
            }
            
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    public void findTaxon(){
    }

    public void findSNPnumber() throws IOException {
        String infile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/output/genotype.vcf";
        BufferedReader br = IOUtils.getTextReader(infile);
        String temp = null;
        int count = 0;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                } else {
                    count = count++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            {
            }
        }
    }

    public static void main (String[] args) throws IOException {
        new Test();
    }
}

