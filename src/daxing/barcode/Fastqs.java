package daxing.barcode;

import format.table.RowTable;
import utils.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import static java.util.stream.Collectors.toSet;

/**
 *
 * @author xudaxing
 */
public class Fastqs {
    private List<Fastq> fastqList;
    private List<String> taxonNameList;
    private boolean isEveryFastqHasOnlyOneFlowcellLaneIndex=false;

    Fastqs(String sampleFastqMapFileS){
        this.readAllFastq(sampleFastqMapFileS);
    }

    Fastqs(File FlowcellLaneIndexFastqDirS){
        this.readFlowcellLaneIndexFastq(FlowcellLaneIndexFastqDirS);
    }

    private void readAllFastq(String sampleFastqMapFileS){
        fastqList=new ArrayList<>();
        taxonNameList=new ArrayList<>();
        RowTable<String> t=new RowTable<>(sampleFastqMapFileS);
        String read1=null;
        String read2=null;
        String sampleID=null;
        Fastq fq=null;
        for(int i=0; i<t.getRowNumber();i++){
            read1=t.getCell(i, 1);
            read2=t.getCell(i, 2);
            sampleID=t.getCell(i, 0);
            fq=new Fastq(read1,read2,sampleID);
            fastqList.add(fq);
            taxonNameList.add(sampleID);
        }
        isEveryFastqHasOnlyOneFlowcellLaneIndex=fastqList.stream().map(Fastq::getAllDifferentFlowcellLaneIndex).allMatch(l -> l.size()==1);
        read1=null;read2=null;sampleID=null;fq=null;
//        List<String> allDifferentFlowcellLaneIndex=fastqList.stream()
//            .map(Fastq::getFlowcellLaneIndexList)
//            .flatMap(List::stream).distinct()
//            .collect(Collectors.toList());
//        Collections.sort(allDifferentFlowcellLaneIndex);
//        for(int i=0;i<allDifferentFlowcellLaneIndex.size();i++){
//
//        }
    }

    private void readFlowcellLaneIndexFastq(File AllFlowcellLaneIndexDirS){
        fastqList=new ArrayList<>();
        taxonNameList=new ArrayList<>();
        File[] file = AllFlowcellLaneIndexDirS.listFiles();
        file=IOUtils.listFilesEndsWith(file, "sfq");
        Arrays.sort(file);
        String fastq1=null;
        String fastq2=null;
        String sampleID=null;
        Fastq fq=null;
        for(int i=0;i<file.length;i=i+2){
            fastq1=file[i].getAbsolutePath();
            fastq2=file[i+1].getAbsolutePath();
            sampleID=file[i].getName().split("_")[3];
            fq=new Fastq(fastq1, fastq2, sampleID);
            this.fastqList.add(fq);
            this.taxonNameList.add(sampleID);
        }
        isEveryFastqHasOnlyOneFlowcellLaneIndex=true;
        fastq1=null;fastq2=null;sampleID=null;fq=null;
    }

    //将所有fastq文件分别按照其各自的flowcellLaneIndex进行输出
    public void writeFlowcellLaneIndexFastq(String outfileDirs){
        fastqList.stream().forEach(fq->{
            fq.writeFlowcellLaneIndexFastq(outfileDirs);
        });
        int count=new File(outfileDirs).listFiles().length;
        System.out.println(outfileDirs+" has "+count/2+" simple fastq files");
    }

    public Set<String> getDifferentFlowcellLaneIndex(){
        Set<String> s=fastqList.stream()
                .map(Fastq::getAllDifferentFlowcellLaneIndex)
                .flatMap(Set::stream)
                .collect(toSet());
        System.out.println("Current sampleFastqMap have "+s.size()+" differnet flowcellLaneIndex");
        return s;
    }

    public List<Fastq> getFastqList(){
        return this.fastqList;
    }

    public boolean getIsEveryFastqHasOnlyOneFlowcellLaneIndex(){
        return this.isEveryFastqHasOnlyOneFlowcellLaneIndex;
    }

}

