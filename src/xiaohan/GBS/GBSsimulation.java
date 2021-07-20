package xiaohan.GBS;

import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.*;

public class GBSsimulation {
    String cutter1 = null;
    String cutter2 = null;
    String index1 = "-1";
    String index2 = "-1";
    String thread = null;
    String outfile = null;
    String outfileDir = null;

    public GBSsimulation(String[] args) {
//        this.parsecutters(args);
        this.test();
    }

    public void test() {
        String read = "TCCGGGC";
        String cutter1 = "CCGG";
        int index1 = read.indexOf(cutter1);
        System.out.println(index1);
    }

    public void parsecutters(String[] args) {
        String infile = "/data1/home/xiaohan/nam/GBSlengthtest/length.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        HashSet<String> nameSet = new HashSet<>();
        HashMap<String,String> cutterMap = new HashMap<>();
        HashMap<String,String> cutterIndexMap = new HashMap<>();
        try{
            String temp = null;
            while((temp = br.readLine())!=null){
                nameSet.add(temp.split("\t")[0]);
                cutterMap.put(temp.split("\t")[0],temp.split("\t")[1]);
                cutterIndexMap.put(temp.split("\t")[0],temp.split("\t")[2]);
            }
            br.close();
        }catch (Exception e){
            e.printStackTrace();
        }
        thread = args[0];
        String[] names = nameSet.toArray(new String[0]);
        for (int i = 0; i < names.length-1; i++) {
            for (int j = i+1; j < names.length; j++) {
                String name1 = names[i];
                String name2 = names[j];
                cutter1 = cutterMap.get(name1);
                cutter2 = cutterMap.get(name2);
                index1 = cutterIndexMap.get(name1);
                index2 = cutterIndexMap.get(name2);
                outfile = name1+"_"+name2;
                outfileDir = new File("/data1/home/xiaohan/nam/GBSlengthtest",outfile).getAbsolutePath();
                File out = new File(outfileDir);
                if(!out.exists()) {
                    out.mkdir();
                }
                this.recordposition(cutter1,cutter2,index1,index2);
                this.lengthsummary(outfileDir);
            }
        }

    }

    public void recordposition(String cutter1,String cutter2,String index1,String index2) {
        HashSet<Integer> nameSet = new HashSet<>();
        List<Future<Command>> resultList = new ArrayList<>();
        ExecutorService pool = Executors.newFixedThreadPool(Integer.parseInt(thread));
        File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());

        for (int i = 1; i < 43; i++) {
            nameSet.add(i);
            String infile = new File("/data1/home/xiaohan/hapscanner/ref", "chr" + i + ".fa").getAbsolutePath();
            String outfile = new File(outfileDir, "chr" + i + ".txt").getAbsolutePath();
            String command = null;
            Command com = new Command(infile, outfile, cutter1, cutter2, index1, index2);
            Future<Command> chrom = pool.submit(com);
        }
        try {
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    static class Command implements Callable<Command> {
        String inputfile = null;
        String outputfile = null;
        String cutter1 = null;
        String cutter2 = null;
        int index1 = -1;
        int index2 = -1;
        int pos = -1;
        int pos1 = -1;
        int pos2 = -1;
        int countline = -1;

        public Command(String infileS, String outfileS, String cutter1, String cutter2, String index1, String index2) {
            this.inputfile = infileS;
            this.outputfile = outfileS;
            this.cutter1 = cutter1;
            this.cutter2 = cutter2;
            this.index1 = Integer.parseInt(index1);
            this.index2 = Integer.parseInt(index2);
        }

        @Override
        public Command call() throws Exception {
            BufferedReader br = IOUtils.getTextReader(inputfile);
            BufferedWriter bw = IOUtils.getTextWriter(outputfile);
            try {
                String temp = null;
                String temp1 = null;
                String temps = null;
                temp1 = br.readLine();
                while ((temp = temp1) != null) {
                    temp1 = br.readLine();
                    if (temp.startsWith(">")) continue;
                    countline++;
                    temps = temp+temp1;
                    pos1 = temps.indexOf(cutter1);
                    pos2 = temps.indexOf(cutter2);
                    if (pos1 != -1) {
                        pos1 = countline * 60 + pos1 + index1;
                    }
                    if (pos2 != -1) {
                        pos2 = countline * 60 + pos2 + index2;
                    }
                    if (pos1 < pos2) {
                        if (pos1 != -1) {
                            bw.write(pos1 + "\n");
                        }
                        if (pos2 != -1) {
                            bw.write(pos2 + "\n");
                        }
                    } else {
                        if (pos2 != -1) {
                            bw.write(pos2 + "\n");
                        }
                        if (pos1 != -1) {
                            bw.write(pos1 + "\n");
                        }
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            return this;
        }
    }

    public void lengthsummary(String outfileDir) {
        for (int i = 1; i < 43; i++) {
            String infile = new File(outfileDir, "chr" + i + ".txt").getAbsolutePath();
            String outfile = new File(outfileDir, "chr" + i + "_length.txt").getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(infile);
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            String temp = null;
            int pos1 = -1;
            int pos2 = -1;
            int length = -1;
            try {
                pos1 = Integer.parseInt(br.readLine());
                while ((temp = br.readLine()) != null) {
                    pos2 = Integer.parseInt(temp);
                    length = pos2 - pos1;
                    if (length > 0) {
                        bw.write(length + "\n");
                    }
                    pos1 = pos2;
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public static void main(String[] args) {
        new GBSsimulation(args);
    }
}
