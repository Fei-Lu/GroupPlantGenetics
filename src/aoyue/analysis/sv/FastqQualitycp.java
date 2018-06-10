/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.sv;

import analysis.maize2k.*;
import format.alignment.ShortreadPEAlignment;
import format.dna.Fastq;
import format.table.RowTable;
import format.table.TableInterface;
import graphcis.r.DensityPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import utils.IOUtils;

/**
 *
 * @author feilu
 */
public class FastqQualitycp {
    
    public FastqQualitycp () {
        //this.covergage();
        this.covergagePlot();
        //this.sampleReads();
        //this.fastQC();
        //this.fastQCsummary();
        //this.contamination();
        //this.trimming();
        //this.alignBWA();
        //this.estimateLibrarySize();
        //this.plastid();
    }
    
    public void plastid () {
        String inputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/alignment/sam";
        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/plastid/summary.txt";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        double[][] plastid = new double[fs.length][2];
        for (int i = 0; i < fs.length; i++) {
            ShortreadPEAlignment sa = new ShortreadPEAlignment ();
            sa.readFromBWAMEM(fs[i].getAbsolutePath());
            for (int j = 0; j < sa.getAlignmentNumber(); j++) {
                if (sa.getHitF(j).equals("11")) plastid[i][0]++;
                if (sa.getHitF(j).equals("12")) plastid[i][1]++;
            }
            plastid[i][0] = plastid[i][0]/sa.getAlignmentNumber();
            plastid[i][1] = plastid[i][1]/sa.getAlignmentNumber();
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Sample\tMt\tPt");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                StringBuilder sb = new StringBuilder(fs[i].getName());
                sb.append("\t").append(plastid[i][0]).append("\t").append(plastid[i][1]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void estimateLibrarySize () {
        String inputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/alignment/sam";
        String outputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/alignment/librarysize";
        new File (outputDirS).mkdir();
        int sampleSize = 3000;
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        int[] insertSize = new int [fs.length];
        for (int i = 0; i < fs.length; i++) {
            ShortreadPEAlignment sa = new ShortreadPEAlignment ();
            sa.readFromBWAMEM(fs[i].getAbsolutePath());
            double[] fragmentSize = new double[sampleSize];
            int cnt = 0;
            for (int j = 0; j < sa.getAlignmentNumber(); j++) {
                if (sa.getMappingQualityF(j) < 30) continue;
                if (sa.getMappingQualityB(j) < 30) continue;
                if (!sa.isPerfectMatchF(j)) continue;
                if (!sa.isPerfectMatchF(j)) continue;
                if (!sa.getHitF(j).equals(sa.getHitB(j))) continue;
                fragmentSize[cnt] = sa.getPEFragmentSize(j);
                cnt++;
                if (cnt == sampleSize) break;
            }
            String outfileS = new File (outputDirS, fs[i].getName()+".pdf").getAbsolutePath();
            DensityPlot d = new DensityPlot (fragmentSize);
            d.setXLim(0, 800);
            d.saveGraph(outfileS);
        }
        
    }
    
    private void alignBWA () {
        String bwaPath = "/Users/feilu/Software/bwa-0.7.15/bwa";
        String indexFileS = "/Users/feilu/Documents/database/maize/reference/bwaLib/maizeAGPv4.fa";
        String inputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/paired";
        String outputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/alignment/sam/";
        String outputPerlS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/alignment/runAlign.pl";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        String[] names = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(names);
        //bwa mem -t  nthreads ref.fa read1.fq read2.fq > aln-pe.sam
        int numCores = Runtime.getRuntime().availableProcessors();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputPerlS);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder(bwaPath+" mem ");
                sb.append(" -t ").append(numCores).append(" ").append(indexFileS).append(" ");
                sb.append(new File(inputDirS, names[i]+"_1.paired.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(inputDirS, names[i]+"_2.paired.fq.gz").getAbsolutePath());
                sb.append(" > ").append(new File(outputDirS, names[i]+".pe.sam").getAbsolutePath());
                String cmd = sb.toString();
                String command = "system(\""+cmd+"\");";
                bw.write(command);
                bw.newLine();
                
            }
            bw.flush();
            bw.close();
            StringBuilder sb = new StringBuilder("perl ");
            sb.append(outputPerlS);
            Runtime run = Runtime.getRuntime();
            Process p = run.exec(sb.toString());
            System.out.println(sb.toString());
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            p.waitFor();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        
    }
    
    private void trimming () {
        //java -jar /Users/feilu/Software/Trimmomatic-0.36/trimmomatic-0.36.jar
        //PE -phred33 -trimlog /Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/log.txt  
        //Users/feilu/Documents/analysisL/pipelineTt/maize2k/sampleSeq/K16BJS0001_1.fq.gz /Users/feilu/Documents/analysisL/pipelineTest/maize2k/sampleSeq/K16BJS0001_2.fq.gz 
        ///Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/r1_paired.txt.gz /Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/r1_unpaired.txt.gz 
        ///Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/r2_paired.txt.gz /Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/r2_unpaired.txt.gz 
        //ILLUMINACLIP:/Users/feilu/Software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10
        String trimJarPath = "/Users/feilu/Software/Trimmomatic-0.36/trimmomatic-0.36.jar";
        String infileDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/sampleSeq";
        String outfilePairedDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/paired";
        String outfileUnpairedDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/unpaired";
        String libraryPath = "/Users/feilu/Software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa";
        String outputFileS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/trimmedSeq/summary.txt";
        File[] fs = new File (infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        String[] names = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(names);
        int numCores = Runtime.getRuntime().availableProcessors();
        new File (outfilePairedDirS).mkdir();
        new File (outfileUnpairedDirS).mkdir();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            bw.write("Samples\tSurvivingRate");
            bw.newLine();
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder("java -jar ");
                sb.append(trimJarPath).append(" PE -phred33 -threads ").append(numCores).append(" ");
                sb.append(new File(infileDirS, names[i]+"_1.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(infileDirS, names[i]+"_2.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(outfilePairedDirS, names[i]+"_1.paired.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(outfileUnpairedDirS, names[i]+"_1.unpaired.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(outfilePairedDirS, names[i]+"_2.paired.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(outfileUnpairedDirS, names[i]+"_2.unpaired.fq.gz").getAbsolutePath()).append(" ");
                sb.append("ILLUMINACLIP:").append(libraryPath).append(":2:30:10");
                String cmd = sb.toString();
                System.out.println(cmd);
                Runtime run = Runtime.getRuntime();
                Process p = run.exec(cmd);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if (!temp.startsWith("Input Read Pairs:")) continue;
                    String[] tem = temp.split("\\(");
                    tem = tem[1].split("\\)");
                    bw.write(names[i]+"\t"+tem[0]);
                    bw.newLine();
                }
                p.waitFor();
                if ((i+1)%10 == 0) System.out.println(String.valueOf(i+1) + " files processed");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
            
        
        
    }
    
    private void contamination () {
        String infileS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/sampleSeq/K16BJS0001_1.fq.gz";
        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/fastQC/contamination/K16BJS0001_1.fq.txt";
        double lowGC = 0.6;
        double highGC = 0.65;      
        Fastq fc = new Fastq (infileS);
        boolean[] ifOut = new boolean[fc.getReadNum()];
        for (int i = 0; i < fc.getReadNum(); i++) {
            double gc = fc.getRead(i).getGCContent();
            if (gc > lowGC && gc < highGC) ifOut[i] = true;
        }
        fc.writeFasta(outfileS, ifOut);
        //no contamination identified
    }
    
    private void fastQCsummary () {
        String inputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/fastQC/autoReport";
        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/fastQC/summary.txt";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".zip");
        String[][] result = null;
        String[] column = null;
        String[] fileNames = new String[fs.length];
        for (int i = 0; i < fs.length; i++) {
            try {
                ZipFile zf = new ZipFile(fs[i].getAbsolutePath());
                Enumeration<? extends ZipEntry> entries = zf.entries();
                while (entries.hasMoreElements()) {
                    ZipEntry ze = entries.nextElement();
                    if (!ze.getName().endsWith("summary.txt")) continue;
                    if ( i == 0) {
                        BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ze), "UTF-8"));
                        String temp = null;
                        ArrayList<String> columnList = new ArrayList();
                        while ((temp = br.readLine()) != null) {
                            String[] tem = temp.split("\t");
                            columnList.add(tem[1]);
                        }
                        column = columnList.toArray(new String[columnList.size()]);
                        result = new String[fs.length][column.length];
                        br.close();
                    }
                    BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ze), "UTF-8"));
                    String temp = null;
                    ArrayList<String> columnList = new ArrayList();
                    int cnt = 0;
                    while ((temp = br.readLine()) != null) {
                        String[] tem = temp.split("\t");
                        result[i][cnt] = tem[0];
                        fileNames[i] = tem[2];
                        cnt++;
                    }
                    br.close();
                    break;    
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("FileName");
            for (int i = 0; i < column.length; i++) {
                sb.append("\t").append(column[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < result.length; i++) {
                sb = new StringBuilder(fileNames[i]);
                for (int j = 0; j < column.length; j++) {
                    sb.append("\t").append(result[i][j]);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void fastQC () {
        String inputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/sampleSeq";
        String outputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/fastQC/autoReport";
        try {
            StringBuilder sb = new StringBuilder("/Users/feilu/Software/FastQC/fastqc");
            File[] fs = new File (inputDirS).listFiles();
            fs =  IOUtils.listFilesEndsWith(fs, ".gz");
            for (int i = 0; i < fs.length; i++) {
                sb.append(" ").append(fs[i].getAbsoluteFile());
            }
            sb.append(" -o ").append(outputDirS);
            String cmd = sb.toString();
            System.out.println(cmd);
            Runtime run = Runtime.getRuntime();
            Process p = run.exec(cmd);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
            System.out.println("Fastqc evalutation is finished at" + outputDirS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    private void sampleReads () {
        String infileDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/rawdata";
        String outputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/sampleSeq";
        int readNum = 100000;
        int startPoint = 1000000;
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfile1 = new File (outputDirS, name+"_1.fq.gz").getAbsolutePath();
            String outfile2 = new File (outputDirS, name+"_2.fq.gz").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfile1);
                BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfile2);
                String temp = null;
                int cnt = 0;
                while ((temp = br1.readLine()) != null) {
                    cnt++;
                    if (cnt < startPoint) {
                        br1.readLine();br1.readLine();br1.readLine();
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                    }
                    else {
                        bw1.write(temp+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                        bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                        for (int i = 0; i < readNum-1; i++) {
                            bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                            bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                        }
                        bw1.flush();bw1.close();
                        bw2.flush();bw2.close();
                        br1.close();
                        br2.close();
                        break;
                    }
                }
                System.out.println(name+ " completed");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            
        });
    }
    
    private void covergagePlot () {
        ///coverage 
        String infileS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/coverage/coverage.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/coverage/coverage.pdf";
        //SimpleTable t = new SimpleTable (infileS);
        TableInterface<String> t = new RowTable (infileS);
        double[] value = t.getColumnAsDoubleArray(3);
        DensityPlot h = new DensityPlot(value);
        h.setTitle("Coverage of 40 maize samples");
        h.setXLab("Coverage");
        h.setYLab("Density");
        h.saveGraph(outfileS);
    }
    
    private void covergage () {
        String infileDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/rawdata";
        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/coverage/coverage.txt";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        int genomeSize = 2135098301;
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        ConcurrentHashMap<String, Integer> fileSizeMap = new ConcurrentHashMap();
        AtomicInteger acnt = new AtomicInteger();
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.fq.gz").getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                String temp = br.readLine();
                temp = br.readLine();
                int len = temp.length();
                acnt.set(len);
                br.close();
                br = IOUtils.getTextGzipReader(infile1);
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    br.readLine();br.readLine();br.readLine();
                }
                fileSizeMap.put(name, cnt);
                System.out.println(name+" has "+String.valueOf(cnt)+ " reads. Read length: " +String.valueOf(len));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Sample\tReadsNum\tReadsLength\tCoverage");
            bw.newLine();
            String[] names = nameSet.toArray(new String[nameSet.size()]);
            Arrays.sort(names);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder(names[i]);
                int readNum = fileSizeMap.get(names[i]);
                int readLength = acnt.intValue();
                sb.append("\t").append(readNum).append("\t").append(readLength).append("\t").append((double)readNum*readLength*2/genomeSize);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
