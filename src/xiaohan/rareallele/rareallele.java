package xiaohan.rareallele;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import pgl.infra.table.RowTable;
import smile.sort.Sort;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

import static java.lang.Integer.parseInt;

public class rareallele {
    public rareallele() {
        String infileS = "all.non-filtering";
        String inputDir ="/data2/xiaohan/nonfiltering/";
        String outputDir = "/data1/home/xiaohan/rareallele/fastQTL/eGenes/nonfiltering-col/";
        //this.getExampleVCF();
        //this.rankGenes();
        //this.getTransNumber();
        this.getVCFposition(infileS, inputDir);
        this.getsubVCF(infileS,inputDir,outputDir);
        //this.checklines();
        //this.countlines();
        //this.test();
        //this.changeSampleName();
        //this.getGTvcf(infileS,outputDir);
        //this.changeName();
        //this.get5kSNPcount();
//        this.SNPcount();
//        this.SNPTable();
//        this.rankcorrelation();
    }

    public void rankcorrelation(){
        String exprfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/exprfile.txt";
        String snpfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/SNPcountTable.txt";
        String outputDirS ="/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/";
        BufferedReader brexpr = IOUtils.getTextReader(exprfile);
        BufferedReader brsnp = IOUtils.getTextReader(snpfile);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS,"rank.txt").getAbsolutePath());
        String tempexpr = null;
        String tempsnp = null;
        String[] tempe = null;
        String[] temps = null;
        String Sample = "B18-E002\tB18-E007\tB18-E008\tB18-E010\tB18-E011\tB18-E014\tB18-E016\tB18-E018\tB18-E023\tB18-E024\tB18-E029\tB18-E032\tB18-E035\tB18-E038\tB18-E043\tB18-E045\tB18-E046\tB18-E049\tB18-E051\tB18-E062\tB18-E065\tB18-E070\tB18-E072\tB18-E074\tB18-E081\tB18-E082\tB18-E087\tB18-E089\tB18-E097\tB18-E099\tB18-E115\tB18-E118\tB18-E124\tB18-E127\tB18-E134\tB18-E138\tB18-E139\tB18-E141\tB18-E152\tB18-E166\tB18-E170\tB18-E180\tB18-E184\tB18-E185\tB18-E188\tB18-E199\tB18-E203\tB18-E204\tB18-E205\tB18-E210\tB18-E214\tB18-E215\tB18-E218\tB18-E219\tB18-E228\tB18-E233\tB18-E236\tB18-E237\tB18-E242\tB18-E244\tB18-E245\tB18-E251\tB18-E252\tB18-E253\tB18-E256\tB18-E262\tB18-E265\tB18-E267\tB18-E270\tB18-E271\tB18-E273\tB18-E277\tB18-E280\tB18-E286\tB18-E288\tB18-E289\tB18-E290\tB18-E298\tB18-E299\tB18-E305\tB18-E306\tB18-E312\tB18-E316\tB18-E318\tB18-E320\tB18-E324\tB18-E330\tB18-E332\tB18-E335\tB18-E337\tB18-E346\tB18-E347\tB18-E348\tB18-E355\tB18-E356\tB18-E357";
        String[] SampleName = Sample.split("\t");
        HashMap<String,Integer> exprRank = new HashMap<>();
        try{
            while((tempsnp=brsnp.readLine())!=null) {
                if (tempsnp.startsWith("#")) {
                    continue;
                }
                temps = tempsnp.split("\t");
                bw.write(temps[0]+"\t");
                for(int i = 1 ;i<temps.length;i++){
                    exprRank.put(SampleName[i-1],Integer.parseInt(temps[i]));
                }
                tempexpr = brexpr.readLine();
                tempe = tempexpr.split("\t");
                for(int j = 1;j<tempe.length;j++){
                    bw.write(exprRank.get(tempe[j])+"\t");
                }
                bw.newLine();
                }
            bw.flush();bw.close();
            brexpr.close();
            brsnp.close();
        }
        catch (Exception e ){
            e.printStackTrace();
        }
    }

    public void SNPTable(){
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/SNPcount.txt";
        String Table = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/Top_512__expressed_genes_median_counts.txt";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/";
        Set<String> SNP = new HashSet<>();
        String geneName = null;
        String temp = null;
        String[] SNPtemp = null;
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir,"SNPcountTable.txt").getAbsolutePath());
        BufferedReader br = IOUtils.getTextReader(infileS);
        try{
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp.toString());
                    bw.newLine();
                    continue;
                }
                SNP.add(temp);
            }
            SNPtemp = SNP.toArray(new String[SNP.size()]);
            RowTable<String> t =new RowTable<>(Table);
            for(int i = 0;i<t.getRowNumber();i++){
                geneName = t.getCell(i,0);
                System.out.println(geneName);
                for(int j = 0;j<SNPtemp.length;j++){
                    if(SNPtemp[j].startsWith(geneName)){
                        bw.write(SNPtemp[j].toString());
                        bw.newLine();
                        continue;
                    }
                }
            }
            bw.flush();bw.close();
            br.close();
        }
        catch (Exception e ){
            e.printStackTrace();
        }

    }

    public void SNPcount() {
        String VCFfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/5kSNP.vcf";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/";
        String expressionfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/phenotypes36.txt";
        BufferedReader brVCF = IOUtils.getTextReader(VCFfile);
        BufferedReader brexpr = IOUtils.getTextReader(expressionfile);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir,"SNPcount.txt").getAbsolutePath());
        String tempVCF = null;
        String[] VCF = null;
        String Sample = "B18-E002\tB18-E007\tB18-E008\tB18-E010\tB18-E011\tB18-E014\tB18-E016\tB18-E018\tB18-E023\tB18-E024\tB18-E029\tB18-E032\tB18-E035\tB18-E038\tB18-E043\tB18-E045\tB18-E046\tB18-E049\tB18-E051\tB18-E062\tB18-E065\tB18-E070\tB18-E072\tB18-E074\tB18-E081\tB18-E082\tB18-E087\tB18-E089\tB18-E097\tB18-E099\tB18-E115\tB18-E118\tB18-E124\tB18-E127\tB18-E134\tB18-E138\tB18-E139\tB18-E141\tB18-E152\tB18-E166\tB18-E170\tB18-E180\tB18-E184\tB18-E185\tB18-E188\tB18-E199\tB18-E203\tB18-E204\tB18-E205\tB18-E210\tB18-E214\tB18-E215\tB18-E218\tB18-E219\tB18-E228\tB18-E233\tB18-E236\tB18-E237\tB18-E242\tB18-E244\tB18-E245\tB18-E251\tB18-E252\tB18-E253\tB18-E256\tB18-E262\tB18-E265\tB18-E267\tB18-E270\tB18-E271\tB18-E273\tB18-E277\tB18-E280\tB18-E286\tB18-E288\tB18-E289\tB18-E290\tB18-E298\tB18-E299\tB18-E305\tB18-E306\tB18-E312\tB18-E316\tB18-E318\tB18-E320\tB18-E324\tB18-E330\tB18-E332\tB18-E335\tB18-E337\tB18-E346\tB18-E347\tB18-E348\tB18-E355\tB18-E356\tB18-E357";
        String[] SampleName = null;
        SampleName = Sample.split("\t");
        String tempexpr = null;
        String[] expr = null;
        HashMap<String, String> TSSMap = new HashMap<>();
        Set<String> TSSset = new HashSet<>();
        String[] TSS = null;
        try {
            while ((tempexpr = brexpr.readLine()) != null) {
                if (tempexpr.startsWith("#")) continue;
                expr = tempexpr.split("\t");
//                System.out.println(expr[1]);
                TSSset.add(expr[1]);
                TSSMap.put(expr[1],expr[3]);
            }
            TSS = TSSset.toArray(new String[TSSset.size()]);
            Arrays.sort(TSS);
            int[][] SNPcount = new int[TSS.length][SampleName.length];

            while((tempVCF = brVCF.readLine())!=null) {
                if (tempVCF.startsWith("##")) continue;
                if (tempVCF.startsWith("#")) {
                    bw.write("SNPcount" + "\t");
                    bw.write(Sample);
                    bw.newLine();
                    continue;
                }
                VCF = tempVCF.split("\t");
                for (int i = 3; i < SampleName.length + 3; i++) {
                    for (int j = 0; j < TSS.length; j++) {
                        int StartPoint = Integer.parseInt(TSS[j]);
                        int distance =  StartPoint - Integer.parseInt(VCF[1]);
//                            if(distance < 5120 && distance > 0 && VCF[i].endsWith("1")){
//                                SNPcount[j][i - 3] = SNPcount[j][i-3] + 1 ;
//                            }
                            if(distance < 5120 && distance > 0 && VCF[i].equals("0/1")){
                                SNPcount[j][i - 3] = SNPcount[j][i-3] + 1 ;
                            }
                            if(distance < 5120 && distance > 0 && VCF[i].equals("1/1")){
                                SNPcount[j][i - 3] = SNPcount[j][i-3] + 2;
                            }
                    }
                }
            }
            for(int m =0;m<TSS.length;m++){
                bw.write(TSSMap.get(TSS[m])+"\t");
                for(int n = 0;n<SampleName.length;n++) {
                    bw.write(SNPcount[m][n] + "\t");
                }
                bw.newLine();
            }
            bw.flush();bw.close();
            brVCF.close();brexpr.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void get5kSNPcount(){
        String VCFfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/genotypes36.vcf";
        String expressionfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/phenotypes36.txt";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/";
        String[] TSS = null;
        BufferedReader brVCF = IOUtils.getTextReader(VCFfile);
        BufferedReader brexpr = IOUtils.getTextReader(expressionfile);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir,"5kSNP.vcf").getAbsolutePath());
        String tempVCF = null;
        String[] VCF = null;
        String tempexpr = null;
        String[] expr = null;
        HashMap<String,String> TSSMap = new HashMap<>();
        Set<String> TSSset = new HashSet<>();
        try{
            while((tempexpr = brexpr.readLine()) != null){
                if(tempexpr.startsWith("#"))continue;
                expr = tempexpr.split("\t");
                TSSset.add(expr[1]);
            }
            TSS = TSSset.toArray(new String[TSSset.size()]);
            Arrays.sort(TSS);
            BufferedWriter [] bw1 = new BufferedWriter[TSS.length];
//            for (int i = 0; i < TSS.length; i++) {
//                bw1[i]= pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "").getAbsolutePath());
//            }
            while((tempVCF = brVCF.readLine())!=null){
                if(tempVCF.startsWith("#")){
                    bw.write(tempVCF);
                    bw.newLine();
                    continue;
                }
                VCF = tempVCF.split("\t");
                int snpsite = parseInt(VCF[1]);
                for(int i = 0;i<TSS.length;i++){
                    int startsite = Integer.parseInt(TSS[i]);
                    int distance = (int) (startsite - snpsite);
                    if(distance > 0 && distance < 5120){
                        for(int m = 0;m<3;m++){
                            bw.write(VCF[m]+"\t");
                        }
                        for(int j = 9;j<VCF.length;j++) {
                            bw.write(VCF[j]+"\t");
                        }
                        bw.newLine();
                    }
                }
            }
            bw.flush();bw.close();
            brexpr.close();brVCF.close();
        }

        catch (Exception e){
            e.printStackTrace();
        }
    }
    public void changeName() {
        String infileS = "/data1/home/xiaohan/rareallele/fastQTL/chr36.vcf";
        String outputDir = "/data1/home/xiaohan/rareallele/fastQTL";
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "Chr36.vcf").getAbsolutePath());
        String temp = null;
        String[] temps = null;
        String[] tems = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < 1; i++) {
                        tems = temps[i].split("]");
                        bw.write("#" + tems[1] + "\t");
                    }
                    for (int i = 1; i < temps.length; i++) {
                        tems = temps[i].split("]");
                        bw.write(tems[1] + "\t");

                    }
                    bw.newLine();
                } else {
                    bw.write(temp);
                    bw.newLine();

                }

            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getGTvcf(String infileS,String outputDir) {
        BufferedReader br = IOUtils.getTextGzipReader(outputDir + infileS+".new.vcf.gz");
        BufferedWriter bw = IOUtils.getTextGzipWriter(new File(outputDir, "genotypes"+infileS+".vcf.gz").getAbsolutePath());
        String temp = null;
        String[] temps = null;
        String[] tems = null;
        String tem = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                    continue;
                } else {
                    temps = temp.split("\t");
                    for (int j = 0; j < 8; j++) {
                        bw.write(temps[j] + "\t");
                    }
                    for (int i = 8; i < temps.length; i++) {
                        tems = temps[i].split(":");
                        if (tems[0].startsWith("0/2") || tems[0].startsWith("0/3")) {
                            tems[0] = "0/1";
                        }
                        if (tems[0].startsWith("1/2") || tems[0].startsWith("1/3") || tems[0].startsWith("2/2") || tems[0].startsWith("2/3") || tems[0].startsWith("3/3")) {
                            tems[0] = "1/1";
                        }
                        bw.write(tems[0]+"\t");
                    }
//                    for(int i = 8;i<9;i++){
//                        tems = temps[i].split(":");
//                        bw.write(tems[0]+"\t");
//                    }
//                    for(int m = 9;m<temps.length;m++) {
//                        bw.write(temps[m]+"\t");
//                    }
//                    for(int i = 9;i<temps.length;i++){
//                        tems = temps[i].split(":");
//                        //bw.write(tems[0]+tems[1]+"\t");
//                        bw.write(tems[0]+"\t");
//                    }
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void changeSampleName() {
        String infileS = "/data2/xiaohan/nonfiltering/getsub/all.non-filtering.new.vcf";
        String outputDir = "/data2/xiaohan/nonfiltering/getsub/";
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "all.nonfilter.vcf").getAbsolutePath());
        String temp = null;
        String[] temps = null;
        String name = "B18-E002,B18-E007,B18-E008,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
        String[] names = name.split(",");
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) {
                    continue;
                } else if(temp.startsWith("#C")){
                    temps = temp.split("\t");
                    for(int i =0;i<9;i++) {
                        bw.write(temps[i]+"\t");
                    }
                    for(int j = 0;j<names.length-1;j++){
                        bw.write(names[j] + "\t");
                    }
                    bw.write(names[names.length-1]);
                    bw.newLine();
                    continue;
                }
                else{
                    bw.write(temp);
                    bw.newLine();
                    continue;
                }
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void test() {
        for (int i = 0; i < 43; i++) {
            System.out.print(i + ".snp.vcf.gz" + "\t");
        }
    }

    public void countlines() {
        String infileS1 = "/data2/junxu/geneSNP/all.vcf.gz";
        String outfileDir = "/data2/xiaohan/test";
        BufferedReader br1 = IOUtils.getTextGzipReader(infileS1);
        int count1 = 0;
        String temp1 = null;
        try {
            while ((temp1 = br1.readLine()) != null) {
                count1++;
                System.out.println(count1);
            }
            br1.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


//    public void countlines(){
//        String infileS1 = "/data2/junxu/geneSNP/0002.N349.snp.vcf.gz";
//        String outfileDir = "/data2/xiaohan/test";
//        BufferedReader br1 = IOUtils.getTextGzipReader(infileS1);
//        BufferedWriter bw = IOUtils.getTextWriter(new File(outfileDir,"result1.txt").getAbsolutePath());
//        int count1 = 0;
//        String temp1 = null;
//        try {
//            if ((temp1 = br1.readLine()) != null) {
//                count1++;
//                bw.write(count1);
//                bw.newLine();
//            }
//            else {
//            bw.flush();bw.close();
//            br1.close();
//            }
//        }
//        catch (Exception e){
//                e.printStackTrace();
//            }
//    }

    public void checklines() {
        String infileS1 = "/data2/junxu/geneSNP/all.vcf.gz";
        String infileS2 = "/data1/home/xiaohan/rareallele/subvcf/colVCF.vcf";
        BufferedReader br1 = IOUtils.getTextGzipReader(infileS1);
        BufferedReader br2 = IOUtils.getTextReader(infileS2);
        BufferedWriter bw = IOUtils.getTextWriter(new File(infileS2, "ifTrue.txt").getAbsolutePath());
        int count1 = 0;
        int count2 = 0;
        String temp1 = null;
        String temp2 = null;
        try {
            while ((temp1 = br1.readLine()) != null) {
                count1++;
            }
            while ((temp2 = br2.readLine()) != null) {
                count2++;
            }
            if (count1 == count2) {
                System.out.println("True");
            } else {
                System.out.println("False");
            }
            bw.flush();
            bw.close();
            br1.close();
            br2.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getsubVCF(String infileS, String inputDir,String outputDir) {
        String index = getVCFposition(infileS, inputDir);
//        String index ="12\n" +
//                "17\n" +
//                "18\n" +
//                "20\n" +
//                "21\n" +
//                "24\n" +
//                "26\n" +
//                "28\n" +
//                "33\n" +
//                "34\n" +
//                "38\n" +
//                "41\n" +
//                "44\n" +
//                "47\n" +
//                "52\n" +
//                "54\n" +
//                "55\n" +
//                "58\n" +
//                "59\n" +
//                "70\n" +
//                "73\n" +
//                "78\n" +
//                "80\n" +
//                "82\n" +
//                "88\n" +
//                "89\n" +
//                "94\n" +
//                "96\n" +
//                "104\n" +
//                "106\n" +
//                "119\n" +
//                "335\n" +
//                "126\n" +
//                "128\n" +
//                "134\n" +
//                "138\n" +
//                "139\n" +
//                "141\n" +
//                "149\n" +
//                "163\n" +
//                "166\n" +
//                "173\n" +
//                "177\n" +
//                "178\n" +
//                "181\n" +
//                "192\n" +
//                "195\n" +
//                "196\n" +
//                "197\n" +
//                "201\n" +
//                "205\n" +
//                "206\n" +
//                "209\n" +
//                "351\n" +
//                "328\n" +
//                "221\n" +
//                "224\n" +
//                "225\n" +
//                "229\n" +
//                "231\n" +
//                "342\n" +
//                "236\n" +
//                "237\n" +
//                "238\n" +
//                "241\n" +
//                "352\n" +
//                "249\n" +
//                "251\n" +
//                "254\n" +
//                "255\n" +
//                "257\n" +
//                "259\n" +
//                "9\n" +
//                "266\n" +
//                "329\n" +
//                "330\n" +
//                "268\n" +
//                "276\n" +
//                "277\n" +
//                "281\n" +
//                "282\n" +
//                "288\n" +
//                "292\n" +
//                "294\n" +
//                "296\n" +
//                "300\n" +
//                "357\n" +
//                "306\n" +
//                "309\n" +
//                "311\n" +
//                "317\n" +
//                "318\n" +
//                "331\n" +
//                "332\n" +
//                "322\n" +
//                "323\n";
        String[] indexes = null;
        indexes = index.split("\t");
        BufferedReader br = IOUtils.getTextGzipReader(inputDir + infileS + ".vcf.gz");
        String temp = null;
        String[] temps = null;
        String name = "B18-E002,B18-E007,B18-E008,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
        String[] names = name.split(",");
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, infileS + ".new.vcf").getAbsolutePath());
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) {
                    continue;
                   
                } else if(temp.startsWith("#C")) {
                     temps = temp.split("\t");
                    bw.write(temps[0]);
                    for (int j = 1; j < 9; j++) {
                        bw.write("\t"+temps[j]);
                    }
                    for(int j = 0;j<names.length;j++){
                        bw.write("\t"+names[j] );
                    }
                    bw.newLine();
                    continue;
                }
                else {
                    temps = temp.split("\t");
                    bw.write(temps[0]+"\t"+temps[1]);
                    bw.write("\t"+"snp_"+temps[1]+"\t");
                    for(int i = 3;i<9;i++){
                        bw.write("\t"+temps[i]);
                    }
                    for (int i = 0; i < indexes.length; i++) {
                        bw.write("\t"+temps[parseInt(indexes[i])]);
                    }
                    bw.newLine();
                    continue;
                }
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public String getVCFposition(String infileS, String inputDir) {
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
        BufferedReader br = IOUtils.getTextGzipReader(inputDir + infileS + ".vcf.gz");
        String temp = null;
        String[] temps = null;
        String[] tempsOrigin = null;
        StringBuilder sb = new StringBuilder();
        ArrayList<String> NameList = new ArrayList<>();
        try {

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#")) {
                    temps = temp.split("\t");
                    tempsOrigin = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].startsWith("AT")) {
                            temps[i] = temps[i].substring(0, 7);
                        }
                    }
                    for (int j = 0; j < Sample.length; j++) {
                        for (int i = 0; i < temps.length; i++) {
                            if (temps[i].equals(Sample[j])) {
                                sb.append(i+"\t");
                                //System.out.println(i);
                                //System.out.println(tempsOrigin[i]+"\t"+Sample[j]);
                            }
                        }
                    }
                }
                if (!temp.startsWith("#")) {
                    break;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return sb.toString();
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
                System.out.println(SampleS[i]+"\t"+transMap.get(SampleS[i])+"\t");
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
