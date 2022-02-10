package daxing;

import daxing.common.utiles.IOTool;
import daxing.v2.loter.LoterMisc;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) {
        String inputDir=args[0];
        String genotypeDir=args[1];
        String taxaInfo=args[2];
        LoterMisc.GroupInTaxaInfo groupInTaxaInfo= LoterMisc.GroupInTaxaInfo.valueOf(args[3]);
        String outDir =args[4];
        LoterMisc.calculateIndividualSitePercent(inputDir, genotypeDir, taxaInfo, groupInTaxaInfo, outDir);
    }

    public static void imputation(String inputDir, String outDir){
        List<File> files= IOTool.getFileListInDirEndsWith(inputDir, ".gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz","_imputation.vcf.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                List<String> temp;
                String line;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                StringBuilder sb = new StringBuilder();
                while ((line=br.readLine())!=null){
                    sb.setLength(0);
                    temp = PStringUtils.fastSplit(line);
                    sb.append(String.join("\t", temp.subList(0, 9))).append("\t");
                    for (int i = 9; i < temp.size(); i++) {
                        if (temp.get(i).startsWith("0/1") || temp.get(i).equals("./.")){
                            sb.append(".|.").append("\t");
                        }else if (temp.get(i).startsWith("0/0")){
                            sb.append("0|0").append("\t");
                        }else {
                            sb.append("1|1").append("\t");
                        }
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void extractFuture(String v2Dir, String ldPrunedDir, String outDir, int threadsNum){
        long start =System.nanoTime();
        List<File> ldPrunedFiles =IOTool.getFileListInDirEndsWith(ldPrunedDir, ".prune.in");
        List<File> vcfFiles = IOTool.getFileListInDirEndsWith(v2Dir, ".vcf.gz");
        String[] outNames = ldPrunedFiles.stream().map(File::getName).map(s -> s.replaceAll(".prune.in",".prune.in.vcf.gz")).toArray(String[]::new);
        List<Callable<Integer>> callableList = new ArrayList<>();
        for (int i = 0; i < vcfFiles.size(); i++) {
            int finalI = i;
            callableList.add(()-> extractFuture(vcfFiles.get(finalI), ldPrunedFiles.get(finalI), new File(outDir, outNames[finalI])));
        }
        ExecutorService executorService = Executors.newFixedThreadPool(threadsNum);
        TIntArrayList exitCodes = new TIntArrayList();
        try {
            List<Future<Integer>> futureList=executorService.invokeAll(callableList);
            executorService.shutdown();
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            for (Future<Integer> future : futureList){
                exitCodes.add(future.get());
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        for (int i = 0; i < exitCodes.size(); i++) {
            if (exitCodes.get(i)!=0){
                System.out.println(vcfFiles.get(i).getName()+" failed");
            }
        }
        System.out.println(Benchmark.getTimeSpanSeconds(start)+ " s");
    }

    public static int extractFuture(File v2GenotypeFile, File ldPrunedFile, File outFile){
        try (BufferedReader br = IOTool.getReader(v2GenotypeFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            List<String> markerIDList = Files.readAllLines(Paths.get(ldPrunedFile.getAbsolutePath()));
            Map<String, Integer> markerIDMap = new HashMap<>(3000000);
            for (String str: markerIDList){
                markerIDMap.put(str, 0);
            }
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            bw.write(line);
            bw.newLine();
            while ((line=br.readLine())!=null){
                temp =PStringUtils.fastSplit(line.substring(0, 30));
                if (markerIDMap.containsKey(temp.get(2))){
                    bw.write(line);
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException ioException) {
            ioException.printStackTrace();
        }
        return 0;
    }
}