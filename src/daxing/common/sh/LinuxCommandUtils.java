package daxing.common.sh;

import daxing.common.DateTime;
import daxing.common.IOTool;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * This class is used to run linux scripts in java
 */
public class LinuxCommandUtils {

    /**
     *
     * @param command one simple linux command with or without options, can't be a compound command: et al. ls | wc
     * @param workingDirectory current working dir
     * @param logFile log file
     * @return 0 indicates normal termination
     */
    public static Integer runOneSH(String command,
                                   String workingDirectory, File logFile){
        List<String> commandList= PStringUtils.fastSplit(command, " ");
        ProcessBuilder processBuilder = new ProcessBuilder(commandList);
        processBuilder.directory(new File(workingDirectory));
        processBuilder.redirectErrorStream(true);
        processBuilder.redirectOutput(logFile);
        Process process;
        int exitCode=-1;
        try {
            process = processBuilder.start();
            exitCode = process.waitFor();
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
        return exitCode;
    }

    /**
     *
     * @param title log title
     * @param shFile sh script, one command per line
     * @param workingDirectory current working dir
     * @param logDir log dir
     * @param threadsNum thread number
     */
    public static void runSH(String title, String shFile, String workingDirectory, String logDir,
                             int threadsNum){
        System.out.println(DateTime.getDateTimeOfNow());
        long start=System.nanoTime();
        List<String> commandList= IOTool.readAllLines(shFile);
        List<File> logFiles= IOTool.getLogFile(title, logDir, commandList.size());
        List<Callable<Integer>> callableList = new ArrayList<>();
        for (int i = 0; i <callableList.size(); i++) {
            int finalI = i;
            callableList.add(()-> runOneSH(commandList.get(finalI), workingDirectory, logFiles.get(finalI)));
        }
        ExecutorService executorService = Executors.newFixedThreadPool(threadsNum);
        List<Integer> exitCodes = new ArrayList<>();
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
        List<String> failCommandList= new ArrayList<>();
        for (int i = 0; i < exitCodes.size(); i++) {
            if (exitCodes.get(i)!=0){
                failCommandList.add(commandList.get(i));
            }
        }
        if (failCommandList.size()==0){
            System.out.println(shFile+" all commands had completed in "+ Benchmark.getTimeSpanSeconds(start)+ " hours");
        }else {
            IOTool.writeAllLines(new File(logDir, "failedRunCommands.sh"), failCommandList);
            System.out.println(shFile+" had "+ failCommandList.size()+ "commands run failed");
            System.out.println("Commands run failed had been written to "+logDir);
            System.out.println(shFile+" total spend "+Benchmark.getTimeSpanHours(start)+ " hours");
        }
    }
}
