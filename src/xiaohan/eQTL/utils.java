package xiaohan.eQTL;

import pgl.infra.range.Range;
import xiaohan.rareallele.IOUtils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class utils {

    public static BufferedReader getTextGzipReader(String infileS) {
        BufferedReader br = null;
        try {
            //br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(infileS))));
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infileS), 65536)), 65536);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return br;
    }

    public static  String[] getGenelist(String infileS){
        BufferedReader brinfo = IOUtils.getTextReader(infileS);
        String info = null;
        HashSet<String> geneSet = new HashSet<>();
        String[] genelist = null;
        try{
            while((info = brinfo.readLine())!=null){
                if(info.startsWith("#"))continue;
                geneSet.add(info.split("\t")[3]);
            }
            genelist = geneSet.toArray(new String[geneSet.size()]);
            brinfo.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
        return genelist;
    }

    public static int getGeneSnp(int chr,int start,int end) throws IOException, InterruptedException {
        int number = 0;
        String vcfdir = "/data2/junxu/genotypeMaf005_87";
        String inputDir = "/data2/xiaohan/tensorQTL/tempvcf1/";
        String pos = chr + ":" + start + "-" + end;;
        StringBuilder sb = new StringBuilder();
        sb.append("tabix " + chr + ".87.B18.maf005.recode.vcf.gz " + pos + " > " + inputDir + chr + "_temp.vcf ");
        String command = sb.toString();
//                        System.out.println(command);
        File dir = new File(new File(vcfdir).getAbsolutePath());
        String[] cmdarry = {"/bin/bash", "-c", command};
        Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
        p.waitFor();
        BufferedReader brtemp = IOUtils.getTextReader(new File(inputDir, chr + "_temp_.vcf").getAbsolutePath());
        String temp = null;
        int countline = 0;
        while ((temp = brtemp.readLine()) != null) {
            countline++;
//                            System.out.println(countline);
        }
        brtemp.close();
        StringBuilder sb2 = new StringBuilder();
        sb2.append("rm " + chr + "_temp.vcf ");
        String command2 = sb2.toString();
//                        System.out.println(command2);
        File dir2 = new File(new File(inputDir).getAbsolutePath());
        String[] cmdarry2 = {"/bin/bash", "-c", command2};
        Process p2 = Runtime.getRuntime().exec(cmdarry2, null, dir2);
        p2.waitFor();
        return countline;
    }

    public static BufferedReader getTextGzipReader (String infileS, int bufferSize) {
        BufferedReader br = null;
        try {
            //br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(infileS))));
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infileS), bufferSize)), bufferSize);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return br;
    }


    public static BufferedWriter getTextGzipWriter (String outfileS) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfileS), 65536)), 65536);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bw;
    }

    public static BufferedWriter getTextGzipWriter (String outfileS, int bufferSize) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfileS), bufferSize)), bufferSize);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bw;
    }

    public static BufferedWriter getTextWriter (String outfileS) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter (new FileWriter(outfileS), 65536);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bw;
    }

    public static BufferedWriter getNIOTextWriter (String outfileS) {
        BufferedWriter bw = null;
        try {
            bw = Files.newBufferedWriter(Paths.get(outfileS), StandardCharsets.UTF_8);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bw;
    }

    public static BufferedReader getTextReader (String infileS) {
        BufferedReader br = null;
        try {
            br = new BufferedReader (new FileReader(infileS), 65536);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return br;
    }

    public static BufferedReader getNIOTextReader (String infileS) {
        BufferedReader br = null;
        try {
            br = Files.newBufferedReader(Paths.get(infileS), StandardCharsets.UTF_8);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return br;
    }

    public static DataOutputStream getBinaryWriter (String outfileS) {
        DataOutputStream dos = null;
        try {
            dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dos;
    }

    public static DataOutputStream getBinaryWriter (String outfileS, int bufferSize) {
        DataOutputStream dos = null;
        try {
            dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), bufferSize));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dos;
    }

    public static DataOutputStream getNIOBinaryWriter (String outfileS) {
        DataOutputStream dos = null;
        try {
            dos = new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(Paths.get(outfileS)), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dos;
    }

    public static DataOutputStream getNIOBinaryWriter (String outfileS, int bufferSize) {
        DataOutputStream dos = null;
        try {
            dos = new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(Paths.get(outfileS)), bufferSize));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dos;
    }

    public static DataInputStream getBinaryReader (String infileS) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dis;
    }

    public static DataInputStream getBinaryReader (String infileS, int bufferSize) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), bufferSize));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dis;
    }

    public static DataInputStream getNIOBinaryReader (String infileS) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(Files.newInputStream(Paths.get(infileS)), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dis;
    }

    public static DataInputStream getNIOBinaryReader (String infileS, int bufferSize) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(Files.newInputStream(Paths.get(infileS)), bufferSize));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dis;
    }

    public static ObjectOutputStream getObjectWriter (String outfileS) {
        ObjectOutputStream oos = null;
        try {
            oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return oos;
    }

    public static ObjectInputStream getObjectReader (String infileS) {
        ObjectInputStream ois = null;
        try {
            ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return ois;
    }

    public static File[] listFilesContains (File[] fAll, String containStr) {
        ArrayList<File> al = new ArrayList();
        for (int i = 0; i < fAll.length; i++) {
            if (fAll[i].getName().contains(containStr)) al.add(fAll[i]);
        }
        return al.toArray(new File[al.size()]);
    }

    public static File[] listFilesStartsWith (File[] fAll, String startStr) {
        ArrayList<File> al = new ArrayList();
        for (int i = 0; i < fAll.length; i++) {
            if (fAll[i].getName().startsWith(startStr)) al.add(fAll[i]);
        }
        return al.toArray(new File[al.size()]);
    }

    public static File[] listFilesEndsWith (File[] fAll, String endStr) {
        ArrayList<File> al = new ArrayList();
        for (int i = 0; i < fAll.length; i++) {
            if (fAll[i].getName().endsWith(endStr)) al.add(fAll[i]);
        }
        return al.toArray(new File[al.size()]);
    }

    /**
     * List all the files in a directory
     * @param dir
     * @return
     */
    public static File[] listRecursiveFiles (File dir) {
        TreeSet<File> fSet = getRecursiveFiles (dir);
        return fSet.toArray(new File[fSet.size()]);
    }

    private static TreeSet<File> getRecursiveFiles (File dir) {
        TreeSet<File> fileTree = new TreeSet();
        for (File entry : dir.listFiles()) {
            if (entry.isFile()) fileTree.add(entry);
            else fileTree.addAll(getRecursiveFiles(entry));
        }
        return fileTree;
    }
}
