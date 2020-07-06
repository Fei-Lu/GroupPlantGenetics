package yafei;
//import utils.IOUtils;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.IOUtils;
/**
 * @author Yafei Guo
 * @create 2020-07-06 2:27 PM
 */
public class para {
    public static List getBamPath(String infileS) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(new File(infileS)));
        String line;
        List path = new ArrayList<>();
        while((line = reader.readLine())!=null){
            path.add(String.valueOf(line));
        }
        return path;
    }
    public static List getBamName(String inFileS) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(new File(inFileS)));
        String line;
        List name = new ArrayList<>();
        while((line = reader.readLine())!=null){
            List list1 = Arrays.asList(line.split("/"));
//            String s = (String) list1.get(list1.size() - 1);
//            List list2 = Arrays.asList(s.split("\\."));
//            String bam = String.valueOf(list2.get(0));
            String bam = String.valueOf(Arrays.asList(((String) list1.get(list1.size() - 1)).split("\\.")).get(0));
            name.add(bam);
        }
        return name;
    }
    public static void main(String[] args) {
        int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38};
        getNewTaxafile("/Users/guoyafei/Desktop/test1/AA_bam.txt",chrA,"/Users/guoyafei/Desktop/test1/parameters_hapScanner_001.txt","/Users/guoyafei/Desktop/test1/out");
//        int chrB[] = {3,4,9,10,15,16,21,22,27,28,33,34,39,40};
//        getNewTaxafile("/data2/xuebo/Projects/Speciation/More_accessions/hapScan/SS_bam.txt",chrB,"/data2/xuebo/Projects/Speciation/More_accessions/hapScan/parameters_hapScanner_001.txt","/data2/xuebo/Projects/Speciation/More_accessions/hapScan");
//        int chrD[] = {5,6,11,12,17,18,23,24,29,30,35,36,41,42};
//        getNewTaxafile("/data2/xuebo/Projects/Speciation/More_accessions/hapScan/DD_bam.txt",chrD,"/data2/xuebo/Projects/Speciation/More_accessions/hapScan/parameters_hapScanner_001.txt","/data2/xuebo/Projects/Speciation/More_accessions/hapScan");
//        int chrAB[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40};
//        getNewTaxafile("/data2/xuebo/Projects/Speciation/More_accessions/hapScan/AABB_bam.txt",chrAB,"/data2/xuebo/Projects/Speciation/More_accessions/hapScan/parameters_hapScanner_001.txt","/data2/xuebo/Projects/Speciation/More_accessions/hapScan");
//        int chrABD[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40,5,6,11,12,17,18,23,24,29,30,35,36,41,42};
//        getNewTaxafile("/data2/xuebo/Projects/Speciation/More_accessions/hapScan/AABBDD_bam.txt",chrABD,"/data2/xuebo/Projects/Speciation/More_accessions/hapScan/parameters_hapScanner_001.txt","/data2/xuebo/Projects/Speciation/More_accessions/hapScan");
    }
    public static void getNewTaxafile(String bamPath, int chr[], String infileS,String outfileS){
        try {
            BufferedWriter bw = null;
            String temp = "";
            int chrA[] = chr;
//            int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38};
            //int chrA[] = {3,4,9,10,15,16,21,22,27,28,33,34,39,40};
            //int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40};
            //int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40,5,6,11,12,17,18,23,24,29,30,35,36,41,42};
            //int chrA[] = {5,6,11,12,17,18,23,24,29,30,35,36,41,42};
            //int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,5,6,11,12,17,18,23,24,29,30,35,36,41,42};
            List path = getBamPath(bamPath);
            List name = getBamName(bamPath);
            for (int j = 0; j < path.size() ; j++) {
                for(int i = 0; i < chrA.length;i++ ){
                    BufferedReader br = IOUtils.getTextReader(infileS);
                    int line = 1;
                    //bw = IOUtils.getTextWriter(outfileS+"/parameters_hapScannerAABBDD_chr"+chrA[i]+".txt");
                    //bw = IOUtils.getTextWriter(outfileS+"/parameters_hapScannerBarley_chr"+chrA[i]+".txt");
                    bw = IOUtils.getTextWriter(outfileS+"/parameters_hapScanner_"+name.get(j)+"_"+chrA[i]+".txt");
                    while((temp = br.readLine())!=null){
                        if(line==9){
                            bw.write(String.valueOf(path.get(j)));
                            //bw.write("/data2/xuebo/Projects/Speciation/E5/hapSacn/TaxaRefBam/TaxaRefBamAABBDD.txt");
                            bw.newLine();
                        }else if(line==11){
                            //bw.write("/data2/xuebo/Projects/Speciation/E5/hapSacn/posAllele/chr"+chrA[i]+".allele.txt");
                            bw.write("/data2/xuebo/Projects/Speciation/More_accessions/hapScan/chr"+chrA[i]+"_posAllele.txt");
                            bw.newLine();
                        }
                        else if(line==13){
                            //bw.write("/data2/xuebo/Projects/Speciation/E5/hapSacn/hapPos/chr"+chrA[i]+".pos.txt");
                            bw.write("//data2/xuebo/Projects/Speciation/More_accessions/hapScan/chr"+chrA[i]+"_pos.txt");
                            bw.newLine();
                        }
                        else if(line==15){
                            bw.write(Integer.toString(chrA[i]));
                            bw.newLine();
                        }else if(line==17){
                            bw.write("/data1/home/xuebo/anaconda3/bin/samtools");
                            bw.newLine();
                        }else if(line==21){
                            bw.write("/data2/xuebo/Projects/Speciation/More_accessions/hapScan/out/outchr"+name.get(j)+"_"+chrA[i]);
                            bw.newLine();
                        }else{
                            bw.write(temp);
                            bw.newLine();
                        }
                        line++;
                    }
                    bw.flush();
                    bw.close();
                    br.close();
                }
            }
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
}
