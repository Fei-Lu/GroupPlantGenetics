package daxing;

import format.table.RowTable;
import utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static java.util.stream.Collectors.groupingBy;

/**
 *
 * @author xudaxing
 */
public class MergeFlowcellLaneIndexFastq {
    Fastqs fastqs =null;
    Map<String, List<Fastq>> flowcellLaneIndexFastqMap =null;

    MergeFlowcellLaneIndexFastq(Fastqs fastqs){
        this.fastqs=fastqs;
        if(fastqs.isEveryFastqHasOnlyOneFlowcellLaneIndex){
            this.readFlowcellLaneIndexFastqMap();
        }
        else{
            System.out.println("There are some fastq match at least two flowcellLaneIndex");
            System.exit(1);
        }
    }

    private void readFlowcellLaneIndexFastqMap(){
        flowcellLaneIndexFastqMap=new HashMap<>();
        List<Fastq> fastqList=this.fastqs.fastqList;
        flowcellLaneIndexFastqMap = fastqList.stream().collect(groupingBy(Fastq::getFlowcellLaneIndex));
    }

    //向所有fastq添加barcode,并返回barcodeFile
    public void addBarcodeAndgetLibraries(RowTable inputWellBarcode, String wellBarcodeOutputFileS, String librariesOutputDirS, String LibrariesFastqMap) {
        try(BufferedWriter bw1 = IOUtils.getTextWriter(wellBarcodeOutputFileS);
            BufferedWriter bw2 =IOUtils.getTextWriter(LibrariesFastqMap)
        ){
            List<Fastq> l =null;
            String[] tempForFlowcellLaneIndex=null;
            String lineForWellBarcodeOutputFileS="SampleID"+"\t"+"Flowcell-ID"+"\t"+"Lane"+"\t"+"Library-index"+"\t"+"Well-ID"+"\t"+"E1-Barcode"+"\t"+"E2-Barcode"+"\t"+"SampleID";
            String lineForLibrariesFastqMap="Flowcell-ID"+"\t"+"Lane"+"\t"+"Library-index"+"\t"+"R1Path"+"\t"+"R2Path";
            bw1.write(lineForWellBarcodeOutputFileS); bw1.newLine();
            bw2.write(lineForWellBarcodeOutputFileS); bw2.newLine();

            for(String str:flowcellLaneIndexFastqMap.keySet()){
                l=flowcellLaneIndexFastqMap.get(str);
                this.addBarcode(str,l,inputWellBarcode,bw1);
            }
            for (String str:flowcellLaneIndexFastqMap.keySet()){
                //merge同一个flowcellLaneIndex对应的所有fastq
                l=flowcellLaneIndexFastqMap.get(str);
                Fastq fq=l.get(0);
                for(int i=1;i<l.size();i++){
                    fq.mergeFastq(l.get(i));
                }
                String r1FastqOutput=new File(librariesOutputDirS, fq.getFlowcellLaneIndex()+"_"+"R1.fastq").getAbsolutePath();
                String r2FastqOutput=new File(librariesOutputDirS, fq.getFlowcellLaneIndex()+"_"+"R2.fastq").getAbsolutePath();
                fq.writeFastq(r1FastqOutput, r2FastqOutput);
                tempForFlowcellLaneIndex=str.split("_");
                bw2.write(tempForFlowcellLaneIndex[0]+"\t"+tempForFlowcellLaneIndex[1]+"\t"+tempForFlowcellLaneIndex[2]+"\t"+r1FastqOutput+"\t"+r2FastqOutput);
                bw2.newLine();
            }
            bw1.flush();bw2.flush();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    private void addBarcode(String flowcellLaneIndex, List<Fastq> fqList, RowTable inputWellBarcode, BufferedWriter wellBarcodeOutputStream) throws IOException {
        String wellID;
        String barcodeForR1;
        String barcodeForR2;
        String line;
        String[] temp=flowcellLaneIndex.split("_");
        if(fqList.size()<=96){
            for(int i=0;i<fqList.size();i++){
                wellID=inputWellBarcode.getCellAsString(i, 0);
                barcodeForR1=inputWellBarcode.getCellAsString(i, 1);
                barcodeForR2=inputWellBarcode.getCellAsString(i, 2);
                fqList.get(i).addBarcode(wellID, barcodeForR1, barcodeForR2);
                line=fqList.get(i).getTaxonName()+"\t"+temp[0]+"\t"+temp[1]+"\t"+temp[2]+"\t"+wellID+"\t"+barcodeForR1+"\t"+barcodeForR2+"\t"+fqList.get(i).getTaxonName();
                wellBarcodeOutputStream.write(line);wellBarcodeOutputStream.newLine();
            }
            System.out.println(flowcellLaneIndex+" has "+fqList.size()+" samples");
        }
        else{
            System.out.println("There are libraries have at least 97 samples");
            System.exit(1);
        }
    }

}

