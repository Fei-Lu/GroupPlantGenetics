/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import com.itextpdf.text.pdf.AcroFields.Item;
import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import pgl.infra.dna.Fastq;
import pgl.infra.dna.Read;
import pgl.infra.table.RowTable;
import static java.lang.Integer.sum;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.io. * ;
import static java.lang.String.format;
import java.util.Collection;
import java.util.Collections;
import static java.util.Collections.list;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import static xujun.analysis.rnaseq.GBS.getKeyList;

/**
 *
 * @author xujun
 */
public class DistinguishSample {

    public String ID;
    Read[] reads=null;
    int phredScale = Integer.MIN_VALUE;
    public DistinguishSample()  {
 //         this.fastqToFasta1();
//          this.testxlsx();
//          this.sampleFastq1 ();
//         this.test();
 //          this.genecount();
//            this.cut12base();
//            this.anyWrong();
            this.changewheatgtf();
    }
     /**
     * Build a byte converter to convert AscII byte following the BaseCoder rules
     * A(00000000), C(00000001), G(00000010), T(0000000011), others(00000100)
     * @return 
     */
    public static void changewheatgtf(){//按照老师的将一条染色体分成两条
        String sampleInfor="/Users/xujun/Desktop/apartchro.txt";
        RowTable rt=new RowTable(sampleInfor);
        HashMap chrnum=new HashMap();
        HashMap numnum = new HashMap () ;
        HashMap chroStart = new HashMap ();
        List<String> h1= new ArrayList();
        for(int i=0;i<rt.getRowNumber();i++){
            chrnum.put( rt.getCell(i, 5),rt.getCell(i, 3));
            numnum.put(rt.getCell(i, 5),rt.getCell(i, 0));    
            chroStart.put(rt.getCell(i, 0),rt.getCell(i, 4));
        }
//        String inputFile="/Users/xujun/Downloads/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3";
        String inputFile="/Users/xujun/Downloads/wheat.gtf";
        String temp=null;String [] tem = null;
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter("/Users/xujun/Downloads/rightchangewheat.gtf");
            while((temp=br.readLine())!=null){  
                List<String> tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if(!tem[0].equals("chrUn")){
                    h1=getKeyList(chrnum,tem[0]);
                    Collections.sort(h1);
                    if(Integer.valueOf(tem[4].toString()) < Integer.valueOf(h1.get(0))){
                        bw.write(String.valueOf(numnum.get( String.valueOf(h1.get(0)))));
                        String start = String.valueOf(chroStart.get(String.valueOf(numnum.get( String.valueOf(h1.get(0))))));
                        String out=temp.replace(tem[0],"");
                        String out1 = out.replace(tem[3],String.valueOf(Integer.valueOf(tem[3])-Integer.valueOf(start)));
                        bw.write(out1.replace(tem[4],String.valueOf(Integer.valueOf(tem[4])-Integer.valueOf(start))));
                        bw.newLine();
                    }else{
                        bw.write(String.valueOf(numnum.get( String.valueOf(h1.get(1)))));
                        String start = String.valueOf(chroStart.get(String.valueOf(numnum.get( String.valueOf(h1.get(1))))));
                        String out=temp.replace(tem[0],"");
                        String out1 = out.replace(tem[3],String.valueOf(Integer.valueOf(tem[3])-Integer.valueOf(start)));
                        bw.write(out1.replace(tem[4],String.valueOf(Integer.valueOf(tem[4])-Integer.valueOf(start))));
//                        bw.write(temp.replace(tem[0],""));
                        bw.newLine();
                    }
                }else{
                    bw.write('0');
                    bw.write(temp.replace(tem[0],""));
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();bw.close();
//            System.out.println(count);
        }
        catch (Exception ex) {
//               System.out.println(count);
               ex.printStackTrace();
               
        }
    }
    public static List<String> getKeyList(HashMap<String,String> map,String value){//the method to get key from value
         List<String> keyList = new ArrayList();
         for(String getKey: map.keySet()){
             if(map.get(getKey).equals(value)){
                 keyList.add(getKey);
             }
         }
         return keyList;
     }
    public static void anyWrong(){
        String inputFile="/Users/xujun/Desktop/PLATE-seq_1.clean.fq";
        String temp=null;int count=0;
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter("/Users/xujun/Desktop/mayWrong.txt");
            while((temp=br.readLine())!=null){  
                count++;
                if(count>=5754881){
                    bw.write(temp);bw.newLine();
                    bw.write(br.readLine());bw.newLine();
                    bw.write(br.readLine());bw.newLine();
                    bw.write(br.readLine());bw.newLine();
                    continue;
                }
                br.readLine();br.readLine();br.readLine();
            }
            System.out.println(count);
        }
        catch (Exception ex) {
               System.out.println(count);
               ex.printStackTrace();
               
        }
    }
    public static HashByteByteMap getAscIIByteMap () {
        int size = 96;
        byte[] key = new byte[size];
        byte[] value = new byte[size];
        for (int i = 0; i < key.length; i++) {
            key[i] = (byte)i;
            value[i] = (byte)4;//所有的value[i]都是4
        }
        value[65] = 0;//A 如果转换成char类型的话
        value[67] = 1;//C
        value[71] = 2;//G
        value[84] = 3;//T
//        HashByteByteMap ascIIByteMap = HashByteByteMaps.newImmutableMap(key, value);
        HashByteByteMap ascIIByteMap = HashByteByteMaps.newMutableMap(key,value);
        for (int i = 0; i < 128; i++) {
            ascIIByteMap.put((byte)i, (byte)4);//Associates the specified value with the specified key in this map.将k和v联系起来
        }
        ascIIByteMap.put((byte)65, (byte)0);//前面就是key 后面就是value  后面这个value就是指的ATCG所代表的byte 只不过是以十进制的形式展示出来 如果是以二进制的话就和上面一样了 只要将key转换成char类型 就是对应的ATCG
        ascIIByteMap.put((byte)67, (byte)1);
        ascIIByteMap.put((byte)71, (byte)2);
        ascIIByteMap.put((byte)84, (byte)3);
        return ascIIByteMap;
    }
 
    public static void getbyte(){
        HashMap<Character,String> hm=new HashMap<Character,String>();
        hm.put('A',"00");
        hm.put('T',"01");
        hm.put('C',"10");
        hm.put('G',"11");
        String[] a ={"ATCGACTC","ATTCCGTA"};
        String []barcode=new String[96];
        for(int t=0;t<96;t++){
            for(int b=0;b<8;b++){
               if(a[t].charAt(b)=='A'){
                    barcode[t]=hm.get('A').concat(barcode[t]);                   
               }else{ if(a[t].charAt(b)=='T'){
                        barcode[t]=hm.get('T').concat(barcode[t]);
                      }else{ if(a[t].charAt(b)=='C'){
                                barcode[t]=hm.get('C').concat(barcode[t]);
                             }else{
                                  barcode[t]=hm.get('G').concat(barcode[t]);
                             }                                                                              
                      }
               }
            }
            System.out.print(t+":"+barcode);
        }
    }
    public static void trans(int num){
        char[] chs={'0','1'};
        char[] arr=new char[16];
        int pos=arr.length;
        for(int i=0;i<2;i++){
            int temp=num&1;
            arr[--pos]=chs[temp];        
        }
        for(int x=pos;x<arr.length;x++){
            System.out.println(arr[x]);
        }        
    }
    public static int[] changtoint(String[] inputFileS){
        int[]barcode=new int[96];
        for(int t=0;t<96;t++){
            for(int b=8,sum=0;b>0;b++){
               if(!(inputFileS[t].charAt(b)=='A')){
                      sum=sum+0;                 
               }else{ if(!(inputFileS[t].charAt(b)=='T')){
                        sum=sum+2^(16-2*b);
                      }else{ if(!(inputFileS[t].charAt(b)=='C')){
                                sum=sum+2^(17-2*b);
                             }else{
                                  sum=sum+2^(16-2*b)+2^(17-2*b);
                             }                                                                              
                      }
               }
               barcode[t]=sum;
            }    
        }
        Arrays.sort(barcode);
        return barcode;   
//      Arrays.binarySearch(bacord, );
    }
    public static void sortReads(ArrayList<String> inputFiles){
//      HashMap<String,Integer> hm=new HashMap<String,Integer>();
        ArrayList<Integer> dl=new ArrayList<Integer>();  
        for(int i=1;i<=96;i++){
            dl.add(new Integer(i));    
        }
        ArrayList<String> ds=new ArrayList<String>();  
        ds.add("a");  
        ds.add("b");  
        ds.add("c");  
 //       comp.rebuildRealation(dl, ds);  
        System.out.println("重构关系之后");  
        for(int i=0;i<dl.size();i++)  
        {  
            System.out.println(dl.get(i)+"\t"+ds.get(i));  
        }  
        
        
    }
     public static void getsampleMap(String[] inputfiles) {
        int size = 96;
        Integer[] key = new Integer[size];
        String[] value = new String[size];
 //       for (int i = 0; i < key.length; i++) {
 //           key[i] = i;
 //          value[i]=inputfiles[i];
 //      }
        HashMap<Integer,String> hm=new HashMap<Integer,String>();
        for (int i = 0; i < 96; i++) {
            hm.put(i, inputfiles[i]);//Associates the specified value with the specified key in this map.将k和v联系起来
        }
    }
     
//     public void testxlsx() throws FileNotFoundException, IOException{
//         String barcodeFileS="E:\\experimental data\\RNA_seq\\3'RNAseq barcode.xlsx";
//         InputStream is = new FileInputStream(barcodeFileS); 
//         Workbook wb = Workbook.getWorkbook(is); //得到工作薄 have some question.but i can not find.
//         Sheet ft = wb.getSheet(0);
//         int rownumber =ft.getRows(); 
//         int columnumber=ft.getColumns(); 
//         HashMap<String,Integer> barcodeIndexMap=new HashMap<>();
//         for(int i=0;i<=ft.getRows();i++){
//             barcodeIndexMap.put(ft.getCell(i,1).getContents(), i);
//         }
//         String s=null;
//         int index=barcodeIndexMap.get(s.subSequence(0,8));
//         BufferedWriter[] bwa=new BufferedWriter[rownumber];//BufferedWriter[]这个数组是用于输出的 根据96个barcode我们会输出96个文件 所以建立了一个barcode输出数组      
//     }
     public void test () throws IOException {//区分barcode的最终版本 运行时间3min8s
        String barcodeFileS = "E:\\experimental data\\RNA_seq\\3'RNAseq barcode.txt";
        String inputDirS ="E:\\experimental data\\RNA_seq\\twice\\clean_data\\";
        String outputDirS = "E:\\experimental data\\RNA_seq\\twice\\small";        
        RowTable<String> rt = new RowTable<>(barcodeFileS);
        int rowNumber = rt.getRowNumber();
        int columnNumber = rt.getColumnNumber();
        HashMap<String, Integer> barcodeIndexMap = new HashMap<>();
        List<String> nameList = new ArrayList<>();
        for (int i = 0; i < rt.getRowNumber(); i++) {
            barcodeIndexMap.put(rt.getCell(i, 1), i);
            nameList.add(rt.getCell(i, 0));
        }
        File[] fs = new File(inputDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        
        nameSet.stream().forEach(p -> {
            String infile1 = new File (inputDirS, p+"_1.clean.fq").getAbsolutePath();
            String outfile=new File(outputDirS,"count").getAbsolutePath();
            String seq = null;
            try {
                BufferedReader br1 = IOUtils.getTextReader(infile1);                               
                BufferedWriter[] bws = new BufferedWriter[rowNumber];//这里只是创建了一个bufferedreader类型的数组 并没有对里面的各个new
                BufferedWriter bw = IOUtils.getTextWriter(outfile);                               
                for (int i = 0; i < bws.length; i++) {
                    String outfileS = new File (outputDirS, nameList.get(i)).getAbsolutePath();
                    bws[i] = IOUtils.getTextWriter(outfileS);
                }
                String temp=null;
                //String seq=null;
                int cnt=0;
                List w=new ArrayList();//记录T的个数 并从中随机抽取100个数 从2，3，4，5四个里面分别产生100个随机数
                while ((temp = br1.readLine()) != null){
                    cnt++;                   
                    seq=br1.readLine(); 
                    Integer index = barcodeIndexMap.get(seq.substring(0, 8));
                    if (index == null) {
                        br1.readLine();br1.readLine();
                        continue;
                    }
                    int i=0;//用charat的方法 运行31s 试了两次都是30s左右
                    for(int a=9;a<seq.length();a++){
                        if(seq.charAt(a)=='T'){
                            i++;
                        }else{
                            if(i<8){
                                i++;
                            }else{
                                seq=seq.substring(i+8);
                                w.add(i);
                                break;
                            }
                        }
                    }
                
                    
                    
//                   System.out.println(seq);
/*                    seq=seq.substring(8);//每一次都取4个碱基变成byte类型 速度很慢 两次都是4min
                    for(int i=0;i+4<seq.length();i=i+4){
                        byte[] oriByte = seq.substring(i,i+4).getBytes();
                        HashByteByteMap ascMap = format.dna.BaseEncoder.getAscIIByteMap();
                        byte[] tranByte = new byte[oriByte.length];
                        for (int a = 0; a < oriByte.length; a++) {
                            tranByte[a] = ascMap.get(oriByte[a]);
                        }
                        int v=0;
                        for (int b = 0; b < tranByte.length; b++) {            
                            v = (v << 2) + tranByte[b];
                        }
                        if(v!=255){
                            seq=seq.substring(i);
                            break;
                        }
                    }*/
/*                    seq=seq.substring(8);//除barcode以外的其他碱基全都变成byte 每次取8个bite 速度要快很多 1min以内
                    byte[] oriByte = seq.substring(0).getBytes();
                    HashByteByteMap ascMap = format.dna.BaseEncoder.getAscIIByteMap();
                    byte[] tranByte = new byte[oriByte.length];
                    for (int a = 0; a < oriByte.length; a++) {
                        tranByte[a] = ascMap.get(oriByte[a]);
                    }
                    for(int i=0;i+4<seq.length();i+=4){                        
                        int v=0;
                        for (int b = i; b < i+4; b++) {            
                            v = (v << 2) + tranByte[b];
                        }
                        if(v!=255){
                                seq=seq.substring(i);
                                break;                           
                        }
                    }*/
                    
                    
                    bws[index].write(temp + "\n");
                    bws[index].write(seq + "\n");
                    bws[index].write(br1.readLine() + "\n");
                    bws[index].write(br1.readLine()+"\n");
                    
                }
                Collections.shuffle(w); 
                int randomSeriesLength = 100;
                List<Integer> randomSeries = w.subList(0, randomSeriesLength);
                System.out.println(randomSeries);
                
                br1.close();
                for (int i = 0; i < bws.length; i++) {
                    bws[i].flush();
                    bws[i].close();
                }
                    
                       
        }
                
        catch (Exception ex) {
               System.out.println(seq+"\t1234");
               ex.printStackTrace();
               
        }
        });      
     }
     public void fastqToFasta1() throws IOException {
         String inputFileS="E:\\experimental data\\RNA_seq\\twice\\clean_data.fq\\s1116-1_HGVHVCCXY_L2_1.clean.fq";
         String outputFileS="E:\\experimental data\\RNA_seq\\twice\\clean_data.fq\\L2_1.clean_test.txt";
//         BufferedReader br = new BufferedReader(new FileReader(inputFileS));
         BufferedReader br = IOUtils.getTextReader(inputFileS);
         String temp=null;
         BufferedWriter bw = null;
         int count = 0;
         bw = new BufferedWriter(new FileWriter(outputFileS), 65536);
         while ((temp = br.readLine()) != null) {
             bw.write(br.readLine());br.readLine();br.readLine();
             bw.newLine();
             count++;
         }
         bw.flush();
         bw.close();
         System.out.println(count);
     }
     public void sampleFastq () {
        String infileDirS = "E:\\experimental data\\RNA_seq\\twice\\clean_data.fq";
        String outputDirS ="E:\\experimental data\\RNA_seq\\twice\\clean_data.fq)";
        String outputFastaDirS = "E:\\experimental data\\RNA_seq\\twice\\clean_data.fq)";
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfile1 = new File (outputDirS, name+"_1.fq.gz").getAbsolutePath();
            String outfile2 = new File (outputDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfileFasta = new File (outputFastaDirS, name+"_1.fa").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfile1);
                BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfile2);
                BufferedWriter bwf = IOUtils.getTextGzipWriter(outfileFasta);
                String temp = null;
                while ((temp = br1.readLine()) != null) {
                    bw1.write(temp);
                    temp = br1.readLine();bw1.write(temp+"\n");
                    br1.readLine();br1.readLine();
                    bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");br2.readLine();br2.readLine();
                    
                   
                    bw1.flush();bw1.close();
                    bw2.flush();bw2.close();
                    bwf.flush();bwf.close();
                    br1.close();
                    br2.close();
                    break;      
                }
                System.out.println(String.valueOf(1) + " reads are sampled from"+ name);
            }
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
     
     private void sampleFastq1 () {
        String infileDirS = "E:\\experimental data\\RNA_seq\\twice\\clean_data";
        String outputDirS ="E:\\experimental data\\RNA_seq\\twice\\clean_data\\test";
        String outputFastaDirS = "E:\\experimental data\\RNA_seq\\twice\\clean_data\\test";
        int readNum = 1000;
        int startPoint = 1000;
        int fastaNum = 1000;
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.clean.fq").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.clean.fq").getAbsolutePath();
            String outfile1 = new File (outputDirS, name+"_1.fq").getAbsolutePath();
            String outfile2 = new File (outputDirS, name+"_2.fq").getAbsolutePath();
            String outfileFasta = new File (outputFastaDirS, name+"_1.fa").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextReader(infile1);
                BufferedReader br2 = IOUtils.getTextReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextWriter(outfile1);
                BufferedWriter bw2 = IOUtils.getTextWriter(outfile2);
                BufferedWriter bwf = IOUtils.getTextWriter(outfileFasta);
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
                            bw1.write(br1.readLine()+"\n");
                            temp = br1.readLine();bw1.write(temp+"\n");
                            bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                            bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                            if (i > fastaNum) continue;
                            bwf.write(">"+String.valueOf(i));
                            bwf.newLine();
                            bwf.write(temp);
                            bwf.newLine();
                        }
                        bw1.flush();bw1.close();
                        bw2.flush();bw2.close();
                        bwf.flush();bwf.close();
                        br1.close();
                        br2.close();
                        break;
                    }
                }
                System.out.println(String.valueOf(readNum) + " reads are sampled from"+ name);
            }
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
    public static void genecount(){
        String outputdir="/Users/xujun/Desktop/genecount/count";
        String inputdir ="/Users/xujun/Desktop/genecount/inputfile";       
        File[] fs = new File(inputdir).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        nameSet.stream().forEach((String name) -> {
            String infile = new File (inputdir, name).getAbsolutePath();
            String outfile = new File (outputdir, name).getAbsolutePath();
            RowTable<String> rt = new RowTable<>(infile);
            int rowNumber = rt.getRowNumber();
            int columnNumber = rt.getColumnNumber();
            int cnt=0;
            int sum=0;
            try {
                BufferedReader br = IOUtils.getTextReader(infile);
                BufferedWriter bw = IOUtils.getTextWriter(outfile);
                br.readLine();br.readLine();br.readLine();br.readLine();                
                while (br.readLine() != null) {
                    for (int i = 3; i < rt.getRowNumber(); i++) {
                        if(!rt.getCell(i, 1).equals("0")){
                            bw.write(/*rt.getCell(i, 0)+" "+*/rt.getCell(i, 1));
                            sum=sum+rt.getCellAsInteger(i, 1);
                            bw.newLine();
                            cnt=cnt+1;
                        }else{
                            br.readLine();
                        }
                    }
                    bw.flush();bw.close();
                    br.close();
                    System.out.println(cnt);
                    System.out.println(sum);
                    break;      
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
    public static void cut12base(){
        String inputfile="/Users/xujun/Desktop/SRR5909940_1.fastq";
        String outputfile="/Users/xujun/Desktop/SRR5909940cut.fastq";
        String seq=null;
        String name=null;
        try {
                BufferedReader br = IOUtils.getTextReader(inputfile);
                BufferedWriter bw = IOUtils.getTextWriter(outputfile);              
                while ((name=br.readLine()) != null) {
                    bw.write(name+"\n");
                    seq=br.readLine().substring(12);
                    bw.write(seq+"\n");
                    bw.write(br.readLine()+"\n");bw.write(br.readLine().substring(12)+"\n");
                         
                }
                br.close();
                bw.flush();bw.close();
                
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        
    }
    public static void main(String[] args) throws IOException, FileNotFoundException{
        new DistinguishSample();
    }
}
