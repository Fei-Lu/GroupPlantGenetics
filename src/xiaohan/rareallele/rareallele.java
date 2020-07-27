package xiaohan.rareallele;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

import static java.lang.Integer.parseInt;

public class rareallele {
    public rareallele() {
        this.getGTvcf();//将VCF转换成GT格式
//        this.getupstreamSNPcount();//根据上游的不同位置call出rare allele count
//        this.SplitPhenoBychr();
//        this.extractTop();//获取高表达的SNP文件
//        this.extractHomologousGene();//获取同源基因相关的文件
//        this.extractsubgenome();//获取亚基因组的文件
//        this.extractTopPheno();//获取高表达的表达文件
    }

    public void extractTopPheno() {
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3";
        String infileS = "countResult3_DESeq2.txt";
//        String outfileS = "Top5000_S3.txt";
        String outfileS = "Top5000-10000_S3.txt";
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7";
//        String infileS = "countResult7_DESeq2.txt";
//        String outfileS = "Top5000_S7.txt";
//        String outfileS = "Top5000-10000_S7.txt";
//        String infor = "Top_5000_median_exp_genes.txt";
        String infor = "Top_5001_to_10000_median_exp_genes.txt";
        try {
            String temp = null;
            String[] temps = null;
            HashSet GeneSet = new HashSet();
            HashMap<String, Integer> geneMap = new HashMap();
            BufferedReader br = IOUtils.getTextReader(new File(inputDir, infor).getAbsolutePath());
            System.out.println(new File(inputDir, infor).getAbsolutePath());
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("x")) continue;
                temps = temp.split("\t");
                System.out.println(temps[1]);
                GeneSet.add(temps[1]);
                geneMap.put(temps[1], Integer.parseInt(temps[0]));
                continue;
            }
            String temp1 = null;
            String[] temps1 = null;
            String[] rank = new String[5000];
            BufferedReader br1 = IOUtils.getTextReader(new File(inputDir, infileS).getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(new File(inputDir, outfileS).getAbsolutePath());
            while ((temp1 = br1.readLine()) != null) {
                temps1 = temp1.split("\t");
                if (temp1.startsWith("g")) {
                    bw.write(temp1);
                    bw.newLine();
                }
                if (GeneSet.contains(temps1[0])) {
                    String gene = temps1[0];
                    int index = geneMap.get(gene);
                    rank[index - 1] = temp1;
                }
            }
            for (int i = 0; i < rank.length; i++) {
                bw.write(rank[i]);
                bw.newLine();
            }
            br1.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void extractTop() {
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/snp_count/major_SNP_site";
//        String inforDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/expressionTable";
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/snp_count/major_SNP_site";
        String inforDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionTable";
        String infor = "Top_5000_median_exp_genes.txt";
//        String infor = "Top_5001_to_10000_median_exp_genes.txt";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("_")[1] + "_" + fs[i].getName().split("_")[2] ;
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String temp = null;
                String[] temps = null;
                HashSet GeneSet = new HashSet();
                HashMap<String, Integer> geneMap = new HashMap();
                BufferedReader br = IOUtils.getTextReader(new File(inforDir, infor).getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("x")) continue;
                    temps = temp.split("\t");
                    GeneSet.add(temps[1]);
                    geneMap.put(temps[1], Integer.parseInt(temps[0]));
                    continue;
                }
                String temp1 = null;
                String[] temps1 = null;
                String[] rank = new String[5000];
                    BufferedReader br1 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
                    BufferedReader br2 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
//                    BufferedWriter bw = IOUtils.getTextWriter(new File(inputDir, "Top5000_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw = IOUtils.getTextWriter(new File(inputDir, "Top5000-10000_" +f+ "_count.txt").getAbsolutePath());
                    while ((temp1 = br2.readLine()) != null) {
                        temps1 = temp1.split("\t");
                        if (temps1[0].startsWith("gene")) {
                            bw.write(temp1);
                            bw.newLine();
                        }
                    }
                    while ((temp1 = br1.readLine()) != null) {
                        temps1 = temp1.split("\t");
                        if (GeneSet.contains(temps1[0])) {
                            String gene = temps1[0];
                            int index = geneMap.get(gene);
                            rank[index - 1] = temp1;
                        }
                    }
                    for (int i = 0; i < rank.length; i++) {
                        bw.write(rank[i]);
                        bw.newLine();
                    }
                    br1.close();
                    br2.close();
                    bw.flush();
                    bw.close();
                    br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
    }

    public void extractHomologousGene() {
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/snp_count/major_SNP_site";
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/snp_count/major_SNP_site";
        String inforDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data";
        String infor = "TheABD.txt";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("_")[1] + "_" + fs[i].getName().split("_")[2] ;
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String temp = null;
                String[] temps = null;
                HashSet GeneSetA = new HashSet();
                HashSet GeneSetB = new HashSet();
                HashSet GeneSetD = new HashSet();
                HashMap<String, Integer> genesubAMap = new HashMap();
                HashMap<String, Integer> genesubBMap = new HashMap();
                HashMap<String, Integer> genesubDMap = new HashMap();
                BufferedReader br = IOUtils.getTextReader(new File(inforDir, infor).getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    GeneSetA.add(temps[1]);
                    GeneSetB.add(temps[2]);
                    GeneSetD.add(temps[3]);
                    genesubAMap.put(temps[1], Integer.parseInt(temps[0]));
                    genesubBMap.put(temps[2], Integer.parseInt(temps[0]));
                    genesubDMap.put(temps[3], Integer.parseInt(temps[0]));
                    continue;
                }
                String temp1 = null;
                String[] temps1 = null;
                String[] rankA = new String[17108];
                String[] rankB = new String[17108];
                String[] rankD = new String[17108];
                    BufferedReader br1 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
                    BufferedReader br2 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw1 = IOUtils.getTextWriter(new File(inputDir, "Ahomo_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw2 = IOUtils.getTextWriter(new File(inputDir, "Bhomo_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw3 = IOUtils.getTextWriter(new File(inputDir, "Dhomo_" +f+ "_count.txt").getAbsolutePath());
                    while ((temp1 = br2.readLine()) != null) {
                        temps1 = temp1.split("\t");
                        if (temps1[0].startsWith("gene")) {
                            bw1.write(temp1);
                            bw1.newLine();
                            bw2.write(temp1);
                            bw2.newLine();
                            bw3.write(temp1);
                            bw3.newLine();
                        }
                    }
                    while ((temp1 = br1.readLine()) != null) {
                        temps1 = temp1.split("\t");
                        if (GeneSetA.contains(temps1[0])) {
                            String gene = temps1[0];
                            int indexA = genesubAMap.get(gene);
                            rankA[indexA - 1] = temp1;
                        }
                        if (GeneSetB.contains(temps1[0])) {
                            String gene = temps1[0];
                            int indexB = genesubBMap.get(gene);
                            rankB[indexB - 1] = temp1;
                        }
                        if (GeneSetD.contains(temps1[0])) {
                            String gene = temps1[0];
                            int indexD = genesubDMap.get(gene);
                            rankD[indexD - 1] = temp1;
                        }
                    }
                    for (int i = 0; i < rankA.length; i++) {
                        bw1.write(rankA[i]);
                        bw1.newLine();
                    }
                    for (int i = 0; i < rankB.length; i++) {
                        bw2.write(rankB[i]);
                        bw2.newLine();
                    }
                    for (int i = 0; i < rankD.length; i++) {
                        bw3.write(rankD[i]);
                        bw3.newLine();
                    }
                    br1.close();
                    br2.close();
                    bw1.flush();bw2.flush();bw3.flush();
                    bw1.close();bw2.close();bw3.close();
                    br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
    
    public void extractsubgenome() {
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/snp_count/major_SNP_site";
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/snp_count/major_SNP_site";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("_")[1] + "_" + fs[i].getName().split("_")[2] ;
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String temp = null;
                String[] temps = null;
                String[] rank = new String[5000];
                    BufferedReader br1 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw1 = IOUtils.getTextWriter(new File(inputDir, "A_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw2 = IOUtils.getTextWriter(new File(inputDir, "B_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw3 = IOUtils.getTextWriter(new File(inputDir, "D_" +f+ "_count.txt").getAbsolutePath());
                    temp = br1.readLine();
                    bw1.write(temp);
                    bw2.write(temp);
                    bw3.write(temp);
                    bw1.newLine();
                    bw2.newLine();
                    bw3.newLine();
                    while((temp = br1.readLine())!= null){
                        if(!temp.startsWith("T"))continue;
                        temps = temp.split("\t");
                        String geneName = temps[0];
                        String subgenome = geneName.toString().substring(8, 9);
                        if(subgenome.equals("A")){
                            bw1.write(temp);
                            bw1.newLine();
                        }
                        if(subgenome.equals("B")){
                            bw2.write(temp);
                            bw2.newLine();
                        }
                        if(subgenome.equals("D")){
                            bw3.write(temp);
                            bw3.newLine();
                        }
                    }
                    br1.close();
                    bw1.flush();bw1.close();
                    bw2.flush();bw2.close();
                    bw3.flush();bw3.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void SplitPhenoBychr() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/expressionboxcox3.txt";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/splitexpression";
//        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionboxcox7.txt";
//        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/splitexpression";
        try {
            BufferedReader br = pgl.infra.utils.IOUtils.getTextReader(infileS);
            BufferedWriter[] bw = new BufferedWriter[43];
            BufferedWriter[] bw1 = new BufferedWriter[43];
            for (int i = 0; i < 43; i++) {
//                bw[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S7expression" + i + ".bed").getAbsolutePath());
//                bw1[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S7expression" + i + ".txt").getAbsolutePath());
                bw[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S3expression" + i + ".bed").getAbsolutePath());
                bw1[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S3expression" + i + ".txt").getAbsolutePath());
            }
            String temp = null;
            String[] temps = null;
            temp = br.readLine();
            for (int i = 0; i < 43; i++) {
                bw[i].write(temp);
                bw[i].newLine();
                bw1[i].write(temp);
                bw1[i].newLine();
            }
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                int count = Integer.parseInt(temps[0]);
                bw[count].write(temp);
                bw[count].newLine();
                bw1[count].write(temp);
                bw1[count].newLine();
            }
            for (int i = 0; i < 43; i++) {
                bw[i].flush();
                bw[i].close();
                bw1[i].flush();
                bw1[i].close();
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    public void getupstreamSNPcount() {
        //GT
        String VCFfileDir = "/data2/xiaohan/GT/S3";
        String expressionfileDir = "/data1/home/xiaohan/rareallele/fastQTL/expression/S3";
//        String VCFfileDir = "/data2/xiaohan/GT/S7";
//        String expressionfileDir = "/data1/home/xiaohan/rareallele/fastQTL/expression/S7";
        String outputDir = "/data1/home/xiaohan/rareallele/rankcorrelation/S3/minor";
//        String outputDir = "/data1/home/xiaohan/rareallele/rankcorrelation/S7/minor";
        File[] fs = new File(VCFfileDir).listFiles();
        fs = pgl.infra.utils.IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        HashSet<String> nameSet = new HashSet<String>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String name = fs[i].getName().split("\\.")[0].split("r")[1].split("G")[0];
            nameSet.add(name);
        }
        nameSet.stream().forEach((String p) -> {
            try {
                System.out.println("this is running :" + p);
                String[] TSS = null;
                BufferedReader brVCF = IOUtils.getTextGzipReader(new File(VCFfileDir, "chr" + p + "GT.vcf.gz").getAbsolutePath());
                BufferedReader brexpr = IOUtils.getTextGzipReader(new File(expressionfileDir, "S3expression" + p + ".bed.gz").getAbsolutePath());
//                BufferedReader brexpr = IOUtils.getTextGzipReader(new File(expressionfileDir, "S7expression" + p + ".bed.gz").getAbsolutePath());
                BufferedWriter[] bwS = new BufferedWriter[12];
                bwS[0] = IOUtils.getTextWriter(new File(outputDir, p + "_5k_all_count.txt").getAbsolutePath());
                bwS[11] = IOUtils.getTextWriter(new File(outputDir, p + "_10k_all_count.txt").getAbsolutePath());
                for (int i = 1; i < 11; i++) {
                    int count = i - 1;
                    bwS[i] = IOUtils.getTextWriter(new File(outputDir, p + "_" + count + "k_" + i + "k_count.txt").getAbsolutePath());
                }
                String tempVCF = null;
                String[] VCF = null;
                String tempexpr = null;
                String[] expr = null;
                HashMap<String, String> TSSMap = new HashMap<>();
                Set<String> TSSset = new HashSet<>();
//                S3
                String name = "B18-E007,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E029,B18-E032,B18-E035,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E348,B18-E355,B18-E356,B18-E357";
//                S7
//                String name = "B18-E002,B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E355,B18-E356,B18-E357";
                String[] names = name.split(",");
                System.out.println("this is reading expression file:" + p);
                while ((tempexpr = brexpr.readLine()) != null) {
                    if (tempexpr.startsWith("#")) continue;
                    expr = tempexpr.split("\t");
                    TSSMap.put(expr[3], expr[1]);
                    TSSset.add(expr[3]);
                }
                TSS = TSSset.toArray(new String[TSSset.size()]);
                Arrays.sort(TSS);
                int[][] count5k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count5k[i][j] = 0;
                    }
                }
                int[][] count10k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count10k[i][j] = 0;
                    }
                }
                int[][] count0_1k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count0_1k[i][j] = 0;
                    }
                }
                int[][] count1_2k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count1_2k[i][j] = 0;
                    }
                }
                int[][] count2_3k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count2_3k[i][j] = 0;
                    }
                }
                int[][] count3_4k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count3_4k[i][j] = 0;
                    }
                }
                int[][] count4_5k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count4_5k[i][j] = 0;
                    }
                }
                int[][] count5_6k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count5_6k[i][j] = 0;
                    }
                }
                int[][] count6_7k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count6_7k[i][j] = 0;
                    }
                }
                int[][] count7_8k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count7_8k[i][j] = 0;
                    }
                }
                int[][] count8_9k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count8_9k[i][j] = 0;
                    }
                }
                int[][] count9_10k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count9_10k[i][j] = 0;
                    }
                }
                System.out.println("complete create int array: " + p);
                while ((tempVCF = brVCF.readLine()) != null) {
                    if (tempVCF.startsWith("#")) {
                        VCF = tempVCF.split("\t");
                        for (int i = 0; i < 12; i++) {
                            bwS[i].write(VCF[0] + "\t");
                            for (int m = 9; m < names.length+9; m++) {
                                bwS[i].write(VCF[m] + "\t");
                            }
                            bwS[i].newLine();
                        }
                        continue;
                    }
                    VCF = tempVCF.split("\t");
                    int snpsite = parseInt(VCF[1]);
                    for (int i = 0; i < TSS.length; i++) {
                        int startsite = Integer.parseInt(TSSMap.get(TSS[i]));
                        int distance = (int) (startsite - snpsite);
                        if (distance >= 0 && distance < 5000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count5k[i][m - 9]++;
                                        count5k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count5k[i][m - 9]++;
                                        count5k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 0 && distance < 10000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count10k[i][m - 9]++;
                                        count10k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count10k[i][m - 9]++;
                                        count10k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 0 && distance < 1000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count0_1k[i][m - 9]++;
                                        count0_1k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count0_1k[i][m - 9]++;
                                        count0_1k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 1000 && distance < 2000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count1_2k[i][m - 9]++;
                                        count1_2k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count1_2k[i][m - 9]++;
                                        count1_2k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 2000 && distance < 3000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count2_3k[i][m - 9]++;
                                        count2_3k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count2_3k[i][m - 9]++;
                                        count2_3k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 3000 && distance < 4000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count3_4k[i][m - 9]++;
                                        count3_4k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count3_4k[i][m - 9]++;
                                        count3_4k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 4000 && distance < 5000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count4_5k[i][m - 9]++;
                                        count4_5k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count4_5k[i][m - 9]++;
                                        count4_5k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 5000 && distance < 6000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count5_6k[i][m - 9]++;
                                        count5_6k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count5_6k[i][m - 9]++;
                                        count5_6k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 6000 && distance < 7000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count6_7k[i][m - 9]++;
                                        count6_7k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count6_7k[i][m - 9]++;
                                        count6_7k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 7000 && distance < 8000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count7_8k[i][m - 9]++;
                                        count7_8k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count7_8k[i][m - 9]++;
                                        count7_8k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 8000 && distance < 9000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count8_9k[i][m - 9]++;
                                        count8_9k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count8_9k[i][m - 9]++;
                                        count8_9k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 9000 && distance < 10000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count9_10k[i][m - 9]++;
                                        count9_10k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count9_10k[i][m - 9]++;
                                        count9_10k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        continue;
                    }
                }
                System.out.println("complete making int array : " + p);
                for (int i = 0; i < TSS.length; i++) {
                    bwS[0].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[0].write(count5k[i][j] + "\t");
                    }
                    bwS[0].newLine();
                }
                System.out.println("complete writing file : " + p + "_5k_all_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[11].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[11].write(count10k[i][j] + "\t");
                    }
                    bwS[11].newLine();
                }
                System.out.println("complete writing file : " + p + "_10k_all_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[1].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[1].write(count0_1k[i][j] + "\t");
                    }
                    bwS[1].newLine();

                }
                System.out.println("complete writing file : " + p + "_0k_1k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[2].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[2].write(count1_2k[i][j] + "\t");
                    }
                    bwS[2].newLine();
                }
                System.out.println("complete writing file : " + p + "_1k_2k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[3].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[3].write(count2_3k[i][j] + "\t");
                    }
                    bwS[3].newLine();
                }
                System.out.println("complete writing file : " + p + "_2k_3k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[4].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[4].write(count3_4k[i][j] + "\t");
                    }
                    bwS[4].newLine();
                }
                System.out.println("complete writing file : " + p + "_3k_4k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[5].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[5].write(count4_5k[i][j] + "\t");
                    }
                    bwS[5].newLine();
                }
                System.out.println("complete writing file : " + p + "_4k_5k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[6].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[6].write(count5_6k[i][j] + "\t");
                    }
                    bwS[6].newLine();
                }
                System.out.println("complete writing file : " + p + "_5k_6k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[7].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[7].write(count6_7k[i][j] + "\t");
                    }
                    bwS[7].newLine();
                }
                System.out.println("complete writing file : " + p + "_6k_7k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[8].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[8].write(count7_8k[i][j] + "\t");
                    }
                    bwS[8].newLine();
                }
                System.out.println("complete writing file : " + p + "_7k_8k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[9].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[9].write(count8_9k[i][j] + "\t");
                    }
                    bwS[9].newLine();
                }
                System.out.println("complete writing file : " + p + "_8k_9k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[10].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[10].write(count9_10k[i][j] + "\t");
                    }
                    bwS[10].newLine();
                }
                System.out.println("complete writing file : " + p + "_9k_10k_count.txt");
                System.out.println("complete writing file : " + p);
                for (int i = 0; i < 12; i++) {
                    bwS[i].flush();
                    bwS[i].close();
                }
                brexpr.close();
                brVCF.close();
            } catch (
                    Exception e) {
                e.printStackTrace();
            }
        });
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

    public void getGTvcf() {
        String inputDir = "/data2/xiaohan/sub7/";
        String outputDir = "/data2/xiaohan/GT/S7";
        HashSet<String> nameSet = new HashSet<String>();
        File[] fs = new File(inputDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "vcf.gz");
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
//            String[] SampleName = fs[i].getName().split("\\.");
//            String Name = SampleName[0].replace(SampleName[0].substring(0, 3), "");
            String Name = fs[i].getName().split("\\.")[0].split("p")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach((String p) -> {
            try {
                System.out.println("Start reading file snp" + p + "vcf.gz");
                BufferedReader br = IOUtils.getTextGzipReader(new File(inputDir, "snp" + p + ".vcf.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(new File(outputDir, "chr" + p + "GT.vcf.gz").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String[] tems = null;
//                S3
//                String name = "B18-E007,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E029,B18-E032,B18-E035,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E348,B18-E355,B18-E356,B18-E357";
//                S7
                String name = "B18-E002,B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E355,B18-E356,B18-E357";
                String[] names = name.split(",");
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##"))
                        continue;
                    if (temp.startsWith("#C")) {
                        temps = temp.split("\t");
                        for (int i = 0; i < 9; i++) {
                            bw.write(temps[i] + "\t");
                        }
                        for (int j = 0; j < names.length - 1; j++) {
                            bw.write(names[j] + "\t");
                        }
                        bw.write(names[names.length - 1]);
                        bw.newLine();
                        continue;
                    }
                    temps = temp.split("\t");
                    for (int i = 0; i < 2; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    //ID
                    bw.write("snp_" + temps[0] + "_" + temps[1] + "\t");
                    //REF ALT
                    for (int i = 3; i < 5; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    //QUAL FILTER
                    for (int i = 5; i < 7; i++) {
                        bw.write("." + "\t");
                    }
                    bw.write("INFO" + "\t" + "GT" + "\t");
                    for (int i = 9; i < temps.length; i++) {
                        tems = temps[i].split(":");
                        bw.write(tems[0] + "\t");
                    }
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public static void main(String[] args) {
        new rareallele();
    }
}
