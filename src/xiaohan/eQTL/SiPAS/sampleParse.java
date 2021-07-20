package xiaohan.eQTL.SiPAS;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;


/**
 *
 * @author xujun
 */
public class sampleParse {
    // The mode of alignment. PE or SE. The default is PE mode.
    String alignType=null;
    //The SampleInformation file (with header), the format is Taxa\tBarcode\tPlateName\tFastqPath
    String sampleInformationFileS = null;
    //The directory of output
    String outputDirS = null;

    List<String> fqFileSList = null;

    List<String>[] barcodeLists = null;

    HashMap<String, String>[] barcodeTaxaMaps = null;

    List<String>[] taxaLists = null;

    List<String>[] barcodeLengthLists = null;

    List<String> allTaxaList = new ArrayList<String>();

    int[] barcodeLengths = null;

    public sampleParse(String arg) {
        this.parseParameters(arg);
        if(alignType=="SE"){
            this.SEmode();
        }else{
            this.PEmode();
        }
    }

    private void SEmode () {
        long startTimePoint = System.nanoTime();
        fqFileSList.parallelStream().forEach(f -> {
            int fqIndex = Collections.binarySearch(this.fqFileSList, f);
            String subFqDirS = new File (this.outputDirS).getAbsolutePath();
            List<String> barcodeList = barcodeLists[fqIndex];
            String[] subFqFileS = new String[barcodeList.size()];
            HashMap<String, String> btMap = barcodeTaxaMaps[fqIndex];
            Set<String> barcodeSet = btMap.keySet();
            BufferedWriter[] bws = new BufferedWriter[subFqFileS.length];
            HashMap<String, BufferedWriter> barcodeWriterMap = new HashMap<>();
            for (int i = 0; i < subFqFileS.length; i++) {
                String taxon = btMap.get(barcodeList.get(i));
                subFqFileS[i] = new File(subFqDirS, taxon+"_R2.fq").getAbsolutePath();
                bws[i] = IOUtils.getTextWriter(subFqFileS[i]);
                barcodeWriterMap.put(barcodeList.get(i), bws[i]);
            }
            int barcodeLength = this.barcodeLengths[fqIndex];
            try {
                BufferedReader br1 = null;
                BufferedReader br2 = null;
                String f2=f.replace("R1","R2");
                if (f.endsWith(".gz")) {
                    br1 = IOUtils.getTextGzipReader(f);
                    br2 = IOUtils.getTextGzipReader(f2);
                }
                else {
                    br1 = IOUtils.getTextReader(f);
                    br2 = IOUtils.getTextGzipReader(f2);
                }
                String temp = null;
                String seq = null;
                String currentBarcode = null;
                BufferedWriter tw = null;
                int cnt = 0;
                int cnt2 = 0;
                while((temp = br1.readLine())!=null){
                    cnt2++;
                    seq = br1.readLine();
                    currentBarcode = seq.substring(0, barcodeLength);
                    int cutIndex = 0;
                    if (barcodeSet.contains(currentBarcode)) {
                        tw = barcodeWriterMap.get(currentBarcode);
                        tw.write(br2.readLine());tw.newLine();
                        tw.write(br2.readLine());tw.newLine();
                        tw.write(br2.readLine());tw.newLine();
                        tw.write(br2.readLine());tw.newLine();
                    }else {
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                        continue;
                    }
                    br1.readLine();br1.readLine();
                }
                StringBuilder sb = new StringBuilder();
                sb.append(cnt).append(" out of ").append(cnt2).append(", ").append(((float)(double)cnt/cnt2)).append(" of total reads were parsed from " + f);
                System.out.println(sb.toString());
                for (int i = 0; i < subFqFileS.length; i++) {
                    bws[i].flush();
                    bws[i].close();
                }
                br1.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("Distinguish samples according to barcode and trim the barcode.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    private void PEmode () {
        long startTimePoint = System.nanoTime();
        fqFileSList.parallelStream().forEach(f -> {
            int fqIndex = Collections.binarySearch(this.fqFileSList, f);
            String subFqDirS = new File (this.outputDirS).getAbsolutePath();
            List<String> barcodeList = barcodeLists[fqIndex];
            String[] subFqFileS = new String[barcodeList.size()];
            HashMap<String, String> btMap = barcodeTaxaMaps[fqIndex];
            Set<String> barcodeSet = btMap.keySet();
            BufferedWriter[] bws1 = new BufferedWriter[subFqFileS.length];
            BufferedWriter[] bws2 = new BufferedWriter[subFqFileS.length];
            HashMap<String, BufferedWriter> barcodeWriterMap1 = new HashMap<>();
            HashMap<String, BufferedWriter> barcodeWriterMap2 = new HashMap<>();
            for (int i = 0; i < subFqFileS.length; i++) {
                String taxon = btMap.get(barcodeList.get(i));
                subFqFileS[i] = new File(subFqDirS, taxon+".fq").getAbsolutePath();
                bws1[i] = IOUtils.getTextWriter(new File(subFqDirS, taxon+"_R1.fq").getAbsolutePath());
                bws2[i] = IOUtils.getTextWriter(new File(subFqDirS, taxon+"_R2.fq").getAbsolutePath());
                barcodeWriterMap1.put(barcodeList.get(i), bws1[i]);
                barcodeWriterMap2.put(barcodeList.get(i), bws2[i]);
            }
            int barcodeLength = this.barcodeLengths[fqIndex];
            try {
                BufferedReader br1 = null;
                BufferedReader br2 = null;
                String f2=f.replace("R1","R2");
                if (f.endsWith(".gz")) {
                    br1 = IOUtils.getTextGzipReader(f);
                    br2 = IOUtils.getTextGzipReader(f2);
                }
                else {
                    br1 = IOUtils.getTextReader(f);
                    br2 = IOUtils.getTextReader(f2);
                }
                String temp = null;
                String seq = null;
                String currentBarcode = null;
                BufferedWriter tw1 = null;
                BufferedWriter tw2 = null;
                int cnt = 0;
                int cnt2 = 0;
                while((temp = br1.readLine())!=null){
                    cnt2++;
                    seq = br1.readLine();
                    currentBarcode = seq.substring(0, barcodeLength);
                    int cutIndex = 0;
                    if (barcodeSet.contains(currentBarcode)) {
                        tw1 = barcodeWriterMap1.get(currentBarcode);
                        tw1.write(temp);tw1.newLine();
                        tw1.write(seq);tw1.newLine();
                        tw1.write(br1.readLine());tw1.newLine();
                        tw1.write(br1.readLine());tw1.newLine();
                        tw2 = barcodeWriterMap2.get(currentBarcode);
                        tw2.write(br2.readLine());tw2.newLine();
                        tw2.write(br2.readLine());tw2.newLine();
                        tw2.write(br2.readLine());tw2.newLine();
                        tw2.write(br2.readLine());tw2.newLine();

                    }
                    else {
                        br1.readLine();br1.readLine();
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                        continue;
                    }
                }
                StringBuilder sb = new StringBuilder();
                sb.append(cnt).append(" out of ").append(cnt2).append(", ").append(((float)(double)cnt/cnt2)).append(" of total reads were parsed from " + f);
                System.out.println(sb.toString());
                for (int i = 0; i < subFqFileS.length; i++) {
                    bws1[i].flush();bws2[i].flush();
                    bws1[i].close();bws2[i].close();
                }
                br1.close();br2.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("Distinguish samples according to barcode and trim the barcode.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }

    private void parseParameters(String infileS) {
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("Three' Expression Profiler (TEP)")) ifOut = true;
            if (!(temp = br.readLine()).equals("Author: Jun Xu, Fei Lu")) ifOut = true;
            if (!(temp = br.readLine()).equals("Email: liuzhongxujun@163.com; flu@genetics.ac.cn")) ifOut = true;
            if (!(temp = br.readLine()).equals("Homepage: http://plantgeneticslab.weebly.com/")) ifOut = true;
            if (ifOut) {
                System.out.println("Thanks for using xiaohan.eQTL.SiPAS.sampleParse.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.alignType=pLineList.get(0);
        this.sampleInformationFileS = pLineList.get(1);
        this.outputDirS = pLineList.get(2);
        new File(this.outputDirS, "subFastqs").mkdir();
    }
    private void processTaxaAndBarcode () {
        RowTable<String> t = new RowTable<>(this.sampleInformationFileS);
        Set<String> fqSet = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            fqSet.add(t.getCell(i, 3));
        }
        fqFileSList = new ArrayList<>(fqSet);
        Collections.sort(fqFileSList);

        barcodeLengths = new int[fqFileSList.size()]; //不同的样本可能带有不同长度的barcode
        barcodeLists = new ArrayList[fqFileSList.size()];
        taxaLists = new ArrayList[fqFileSList.size()];
        barcodeTaxaMaps = new HashMap[fqFileSList.size()];
        int[] cnts = new int[fqFileSList.size()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Collections.binarySearch(fqFileSList, t.getCell(i, 3));
            cnts[index]++;
        }
        for (int i = 0; i < cnts.length; i++) {
            barcodeLists[i] = new ArrayList<>();
            taxaLists[i] = new ArrayList<>();
            barcodeTaxaMaps[i] = new HashMap<>();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Collections.binarySearch(fqFileSList, t.getCell(i, 3));
            String taxon = t.getCell(i, 0) + "_"+ t.getCell(i, 2);//这个连接起来的taxon就是我们要找的index信息
            allTaxaList.add(taxon);
            taxaLists[index].add(taxon);
            barcodeLists[index].add(t.getCell(i, 1));
            barcodeTaxaMaps[index].put(t.getCell(i, 1), taxon);
            barcodeLengths[index] = t.getCell(i, 1).length();
        }
        Collections.sort(allTaxaList);
    }

    public static void main(String args[]) {
        new sampleParse(args[0]);
    }
}
