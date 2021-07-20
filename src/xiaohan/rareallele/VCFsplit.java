package xiaohan.rareallele;

import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class VCFsplit {

    public VCFsplit() {
//        this.split();
//        this.splitGerp();
//        this.splitgff3();
        this.splitphyloP();
    }

    public void splitgff3(){

        String positionFileS =
                        "1	0	471304005	chr1A	0	471304005\n" +
                        "3	0	438720154	chr1B	0	438720154\n" +
                        "5	0	452179604	chr1D	0	452179604\n" +
                        "7	0	462376173	chr2A	0	462376173\n" +
                        "9	0	453218924	chr2B	0	453218924\n" +
                        "11	0	462216879	chr2D	0	462216879\n" +
                        "13	0	454103970	chr3A	0	454103970\n" +
                        "15	0	448155269	chr3B	0	448155269\n" +
                        "17	0	476235359	chr3D	0	476235359\n" +
                        "19	0	452555092	chr4A	0	452555092\n" +
                        "21	0	451014251	chr4B	0	451014251\n" +
                        "23	0	451004620	chr4D	0	451004620\n" +
                        "25	0	453230519	chr5A	0	453230519\n" +
                        "27	0	451372872	chr5B	0	451372872\n" +
                        "29	0	451901030	chr5D	0	451901030\n" +
                        "31	0	452440856	chr6A	0	452440856\n" +
                        "33	0	452077197	chr6B	0	452077197\n" +
                        "35	0	450509124	chr6D	0	450509124\n" +
                        "37	0	450046986	chr7A	0	450046986\n" +
                        "39	0	453822637	chr7B	0	453822637\n" +
                        "41	0	453812268	chr7D	0	453812268\n";
        HashMap<String, Integer> ChromosomeMap = new HashMap<>();
        HashMap<String, Integer> ChromosomeLength = new HashMap<>();
        String[] temps = positionFileS.split("\n");
        String[] temp = null;
        for (int i = 0; i < temps.length; i++) {
            temp = temps[i].split("\t");
            ChromosomeMap.put(temp[3], Integer.parseInt(temp[0]));
            ChromosomeLength.put(temp[3], Integer.parseInt(temp[5]));
        }
        String infile = "/data2/xiaohan/Transposon/iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3";
        String outputDir = "/data2/xiaohan/Transposon/chr";
            try {
                BufferedReader br = IOUtils.getTextReader(infile);
                String temp1 = null;
                String[] temps1 = null;
                BufferedWriter[] bw = new BufferedWriter[43];
                for (int i = 0; i < 42; i++) {
                    int chr = i +1;
                    bw[i] = IOUtils.getTextWriter(new File(outputDir,"chr"+chr+"_Transposon.gff3").getAbsolutePath());
                }
                bw[42] = IOUtils.getTextWriter(new File(outputDir,"chrUn_Transposon.gff3").getAbsolutePath());
                while ((temp1 = br.readLine()) != null) {
                    if(temp1.startsWith("#")){
                        for(int i = 0;i<42;i++){
                            bw[i].write(temp1);
                            bw[i].newLine();
                        }
                        continue;
                    }
                    if(temp1.startsWith("chrUn")){
                        bw[42].write(temp1);
                        bw[42].newLine();
                        continue;
                    }
                    temps1 = temp1.split("\t");
                    String chrABD = temps1[0];
                    int chrNumber = ChromosomeMap.get(chrABD);
                    int positionstart = Integer.parseInt(temps1[3]);
                    int positionend = Integer.parseInt(temps1[4]);
                    if (ChromosomeLength.get(chrABD) >= positionstart) {
                        bw[chrNumber - 1].write(chrNumber+"\t");
                        for (int i = 1; i < temps1.length; i++) {
                            bw[chrNumber - 1].write(temps1[i] +"\t");
                        }
                        bw[chrNumber - 1].newLine();
                    } else if (ChromosomeLength.get(chrABD) <= positionstart) {
                        int posstart = positionstart - ChromosomeLength.get(chrABD);
                        int posend = positionend - ChromosomeLength.get(chrABD);
                        bw[chrNumber].write(chrNumber+1 + "\t");
                        for (int i = 1; i < 3; i++) {
                            bw[chrNumber].write(temps1[i] + "\t");
                        }
                        bw[chrNumber].write(posstart +"\t"+posend+"\t");
                        for (int i = 5; i < temps1.length; i++) {
                            bw[chrNumber].write(temps1[i] + "\t");
                        }
                        bw[chrNumber].newLine();
                    }
                }
                br.close();
                for (int i = 0; i < 43; i++) {
                    bw[i].flush();bw[i].close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
//        });

    }

    public void splitphyloP(){
        String positionFileS =
                "1	0	471304005	chr1A	0	471304005\n" +
                        "3	0	438720154	chr1B	0	438720154\n" +
                        "5	0	452179604	chr1D	0	452179604\n" +
                        "7	0	462376173	chr2A	0	462376173\n" +
                        "9	0	453218924	chr2B	0	453218924\n" +
                        "11	0	462216879	chr2D	0	462216879\n" +
                        "13	0	454103970	chr3A	0	454103970\n" +
                        "15	0	448155269	chr3B	0	448155269\n" +
                        "17	0	476235359	chr3D	0	476235359\n" +
                        "19	0	452555092	chr4A	0	452555092\n" +
                        "21	0	451014251	chr4B	0	451014251\n" +
                        "23	0	451004620	chr4D	0	451004620\n" +
                        "25	0	453230519	chr5A	0	453230519\n" +
                        "27	0	451372872	chr5B	0	451372872\n" +
                        "29	0	451901030	chr5D	0	451901030\n" +
                        "31	0	452440856	chr6A	0	452440856\n" +
                        "33	0	452077197	chr6B	0	452077197\n" +
                        "35	0	450509124	chr6D	0	450509124\n" +
                        "37	0	450046986	chr7A	0	450046986\n" +
                        "39	0	453822637	chr7B	0	453822637\n" +
                        "41	0	453812268	chr7D	0	453812268\n";
        HashMap<String, Integer> ChromosomeMap = new HashMap<>();
        HashMap<String, Integer> ChromosomeLength = new HashMap<>();
        String[] temps = positionFileS.split("\n");
        String[] temp = null;
        for (int i = 0; i < temps.length; i++) {
            temp = temps[i].split("\t");
            ChromosomeMap.put(temp[3], Integer.parseInt(temp[0]));
            ChromosomeLength.put(temp[3], Integer.parseInt(temp[5]));
        }
        String infileDir = "/data1/home/lipeng/data/05_conservation/91way/phyloP/";
        String outputDir = "/data1/home/xiaohan/annotation/phyloP";
        File[] fs = new File(infileDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "monocots_phyloP.bed.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("_")[0];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextGzipReader(new File(infileDir, f + "_monocots_phyloP.bed.gz").getAbsolutePath());
                String temp1 = null;
                String[] temps1 = null;
                BufferedWriter[] bw = new BufferedWriter[2];
                String name = "chr"+f.toString();
                System.out.println(name);
                int number = ChromosomeMap.get(name);
                int number1 = number +1;
                System.out.println(number);
                bw[0] = IOUtils.getTextWriter(new File(outputDir, "chr" + number + ".bed").getAbsolutePath());
                bw[1] = IOUtils.getTextWriter(new File(outputDir, "chr" + number1 + ".bed").getAbsolutePath());
                while ((temp1 = br.readLine()) != null) {
                    temps1 = temp1.split("\t");
                    String chrABD = temps1[0];
                    int chrNumber = ChromosomeMap.get(chrABD);
                    int position = Integer.parseInt(temps1[1]);
                    if (ChromosomeLength.get(chrABD) >= position) {
                        bw[0].write(number+"\t");
                        bw[0].write(temps1[1]+"\t"+temps1[2]+"\t"+temps1[3]+"\t"+temps1[4]);
                        bw[0].newLine();
                    } else if (ChromosomeLength.get(chrABD) <= position) {
                        int pos = position - ChromosomeLength.get(chrABD);
//                        int pos = position;
                        int pos1 = pos +1;
                        bw[1].write(number1+"\t"+pos+"\t"+pos1+"\t"+temps1[3]+"\t"+temps1[4]);
                        bw[1].newLine();
                    }
                }
                br.close();
                bw[0].flush();bw[0].close();
                bw[1].flush();bw[1].close();
                System.out.println("Finished split file "+ f);
                StringBuilder sb = new StringBuilder();
                sb.append("bgzip chr"+ number +".bed"+" && bgzip chr"+ number1 +".bed");
                String command = sb.toString();
                try {
                    File dir = new File(new File(outputDir).getAbsolutePath());
                    String[] cmdarry = {"/bin/bash", "-c", command};
                    Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                    p.waitFor();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
    }

    public void splitGerp() {
        String positionFileS =
                "1	0	471304005	chr1A	0	471304005\n" +
                        "3	0	438720154	chr1B	0	438720154\n" +
                        "5	0	452179604	chr1D	0	452179604\n" +
                        "7	0	462376173	chr2A	0	462376173\n" +
                        "9	0	453218924	chr2B	0	453218924\n" +
                        "11	0	462216879	chr2D	0	462216879\n" +
                        "13	0	454103970	chr3A	0	454103970\n" +
                        "15	0	448155269	chr3B	0	448155269\n" +
                        "17	0	476235359	chr3D	0	476235359\n" +
                        "19	0	452555092	chr4A	0	452555092\n" +
                        "21	0	451014251	chr4B	0	451014251\n" +
                        "23	0	451004620	chr4D	0	451004620\n" +
                        "25	0	453230519	chr5A	0	453230519\n" +
                        "27	0	451372872	chr5B	0	451372872\n" +
                        "29	0	451901030	chr5D	0	451901030\n" +
                        "31	0	452440856	chr6A	0	452440856\n" +
                        "33	0	452077197	chr6B	0	452077197\n" +
                        "35	0	450509124	chr6D	0	450509124\n" +
                        "37	0	450046986	chr7A	0	450046986\n" +
                        "39	0	453822637	chr7B	0	453822637\n" +
                        "41	0	453812268	chr7D	0	453812268\n";
        HashMap<String, Integer> ChromosomeMap = new HashMap<>();
        HashMap<String, Integer> ChromosomeLength = new HashMap<>();
        String[] temps = positionFileS.split("\n");
        String[] temp = null;
        for (int i = 0; i < temps.length; i++) {
            temp = temps[i].split("\t");
            ChromosomeMap.put(temp[3], Integer.parseInt(temp[0]));
            ChromosomeLength.put(temp[3], Integer.parseInt(temp[5]));
        }
        String infileDir = "/data2/xiaohan/GerpOrigin";
        String outputDir = "/data2/xiaohan/GerpOrigin/Chr_constitute";
        File[] fs = new File(infileDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("\\.")[0];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextGzipReader(new File(infileDir, f + ".bed.gz").getAbsolutePath());
                String temp1 = null;
                String[] temps1 = null;
                BufferedWriter[] bw = new BufferedWriter[2];
                String name = "chr"+f.toString();
                System.out.println(name);
                int number = ChromosomeMap.get(name);
                int number1 = number +1;
                System.out.println(number);
                bw[0] = IOUtils.getTextWriter(new File(outputDir, "chr" + number + ".bed").getAbsolutePath());
                bw[1] = IOUtils.getTextWriter(new File(outputDir, "chr" + number1 + ".bed").getAbsolutePath());
                while ((temp1 = br.readLine()) != null) {
                    temps1 = temp1.split("\t");
                    String chrABD = "chr" + temps1[0];
                    int chrNumber = ChromosomeMap.get(chrABD);
                    int position = Integer.parseInt(temps1[1]);
                    if (ChromosomeLength.get(chrABD) >= position) {
                        bw[0].write(number+"\t");
                        bw[0].write(temps1[1]+"\t"+temps1[2]+"\t"+temps1[3]+"\t"+temps1[4]);
                        bw[0].newLine();
                    } else if (ChromosomeLength.get(chrABD) <= position) {
                        int pos = position - ChromosomeLength.get(chrABD);
//                        int pos = position;
                        int pos1 = pos +1;
                        bw[1].write(number1+"\t"+pos+"\t"+pos1+"\t"+temps1[3]+"\t"+temps1[4]);
                        bw[1].newLine();
                    }
                }
                br.close();
                bw[0].flush();bw[0].close();
                bw[1].flush();bw[1].close();
                System.out.println("Finished split file "+ f);
                StringBuilder sb = new StringBuilder();
                sb.append("bgzip chr"+ number +".bed"+" && bgzip chr"+ number1 +".bed");
                String command = sb.toString();
                try {
                    File dir = new File(new File(outputDir).getAbsolutePath());
                    String[] cmdarry = {"/bin/bash", "-c", command};
                    Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                    p.waitFor();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        });

    }


    public void split() {
        String positionFileS =
                "1	0	471304005	chr1A	0	471304005\n" +
                        "2	0	122798051	chr1A	471304005	594102056\n" +
                        "3	0	438720154	chr1B	0	438720154\n" +
                        "4	0	251131716	chr1B	438720154	689851870\n" +
                        "5	0	452179604	chr1D	0	452179604\n" +
                        "6	0	43273582	chr1D	452179604	495453186\n" +
                        "7	0	462376173	chr2A	0	462376173\n" +
                        "8	0	318422384	chr2A	462376173	780798557\n" +
                        "9	0	453218924	chr2B	0	453218924\n" +
                        "10	0	348037791	chr2B	453218924	801256715\n" +
                        "11	0	462216879	chr2D	0	462216879\n" +
                        "12	0	189635730	chr2D	462216879	651852609\n" +
                        "13	0	454103970	chr3A	0	454103970\n" +
                        "14	0	296739669	chr3A	454103970	750843639\n" +
                        "15	0	448155269	chr3B	0	448155269\n" +
                        "16	0	382674495	chr3B	448155269	830829764\n" +
                        "17	0	476235359	chr3D	0	476235359\n" +
                        "18	0	139317064	chr3D	476235359	615552423\n" +
                        "19	0	452555092	chr4A	0	452555092\n" +
                        "20	0	292033065	chr4A	452555092	744588157\n" +
                        "21	0	451014251	chr4B	0	451014251\n" +
                        "22	0	222603248	chr4B	451014251	673617499\n" +
                        "23	0	451004620	chr4D	0	451004620\n" +
                        "24	0	58852447	chr4D	451004620	509857067\n" +
                        "25	0	453230519	chr5A	0	453230519\n" +
                        "26	0	256543224	chr5A	453230519	709773743\n" +
                        "27	0	451372872	chr5B	0	451372872\n" +
                        "28	0	261776885	chr5B	451372872	713149757\n" +
                        "29	0	451901030	chr5D	0	451901030\n" +
                        "30	0	114179647	chr5D	451901030	566080677\n" +
                        "31	0	452440856	chr6A	0	452440856\n" +
                        "32	0	165638404	chr6A	452440856	618079260\n" +
                        "33	0	452077197	chr6B	0	452077197\n" +
                        "34	0	268911281	chr6B	452077197	720988478\n" +
                        "35	0	450509124	chr6D	0	450509124\n" +
                        "36	0	23083594	chr6D	450509124	473592718\n" +
                        "37	0	450046986	chr7A	0	450046986\n" +
                        "38	0	286659250	chr7A	450046986	736706236\n" +
                        "39	0	453822637	chr7B	0	453822637\n" +
                        "40	0	296797748	chr7B	453822637	750620385\n" +
                        "41	0	453812268	chr7D	0	453812268\n" +
                        "42	0	184873787	chr7D	453812268	638686055";

        HashMap<Integer, Integer> ChromosomeMap = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> ChromosomeLength = new HashMap<Integer, Integer>();
        String[] temps = positionFileS.split("\n");
        String[] temp = null;
        for (int i = 0; i < temps.length; i++) {
            temp = temps[i].split("\t");
            ChromosomeMap.put(Integer.parseInt(temp[0]), Integer.parseInt(temp[5]));
            ChromosomeLength.put(Integer.parseInt(temp[0]), Integer.parseInt(temp[4]));
        }
        try {
            String infile = "/Users/yxh/Documents/RareAllele/004test/test/genotype20200115.txt";
            BufferedReader br = IOUtils.getTextReader(infile);
            BufferedWriter bw;
            bw = IOUtils.getTextWriter("/Users/yxh/Documents/RareAllele/004test/test/genotype-2.txt");
            String temp1 = null;
            while ((temp1 = br.readLine()) != null) {
                if (temp1.startsWith("#")) {
                    //System.out.println(temp1);
                    bw.write(temp1);
                    bw.newLine();
                    continue;
                }
                String[] temp2 = temp1.split("\t");
                int chr = Integer.parseInt(temp2[0]);
                int chrNumber = chr * 2 - 1;
                int position = Integer.parseInt(temp2[1]);
                if (ChromosomeMap.get(chrNumber) >= position) {
                    chr = chrNumber;
                    position = position;
                    String name = "SNP_" + String.valueOf(chr) + "_" + String.valueOf(position);
                    StringBuffer sb = new StringBuffer();
                    sb.append(chr + "\t" + position + "\t" + name + "\t");
                    for (int i = 3; i < temp2.length; i++) {
                        sb.append(temp2[i]);
                        sb.append("\t");
                    }
                    //System.out.println(sb.toString());
                    bw.write(sb.toString());
                    bw.newLine();
                    continue;
                } else if (ChromosomeMap.get(chrNumber) <= position) {
                    chr = chrNumber + 1;
                    position = position - ChromosomeLength.get(chrNumber + 1);
                    String name = "SNP_" + String.valueOf(chr) + "_" + String.valueOf(position);
                    StringBuffer sb = new StringBuffer();
                    sb.append(chr + "\t" + position + "\t" + name + "\t");
                    for (int i = 3; i < temp2.length; i++) {
                        sb.append(temp2[i]);
                        sb.append("\t");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    continue;
                }
            }
            br.close();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new VCFsplit();
    }
}
