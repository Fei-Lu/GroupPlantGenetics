/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import format.genomeAnnotation.GeneFeature;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import jxl.read.biff.BiffException;

/**
 *
 * @author xujun
 */
public class RNAseq {
    public RNAseq() throws IOException, FileNotFoundException, BiffException {
//            this.getgenename();
//           this.genenumber();
          this.notSDI();
    }      
    public void getgenename(){
         String kfffile="/Users/xujun/Desktop/Zea_mays.AGPv4.38.kgf";
         GeneFeature fs=new GeneFeature(kfffile);
         String inputfile="/Users/xujun/Desktop/news.fq"; 
         String outputfile="/Users/xujun/Desktop/notSIDname.text";
         BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
         BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);
         String a=null;
         int index=0;
         String name=null;
         String yon=null;
         try{
            while((a=br1.readLine())!=null){
                index=fs.getGeneIndex(Integer.parseInt(a.split("\t")[0]),Integer.parseInt(a.split("\t")[1]) );
                if (index < 0) continue;
                name=fs.getGeneName(index);
                bw.write(name+"\t"+a.split("\t")[2]+"\t"+a.split("\t")[3]);               
                bw.newLine();
                }
            br1.close();
            bw.flush();bw.close();
         }
         catch(Exception ex){
             System.out.println(a);
             ex.printStackTrace();
         }        
     }
    public void genenumber(){//对map到的每个基因上的reads进行计数 第一列基因名 第二列计数 第三列平均长度 再加一个第四列测序质量值的均值
         String inputfile="/Users/xujun/Desktop/notSIDname.text";
         String outputfile="/Users/xujun/Desktop/notSDIgenecount.text";
         String news=null;
         List w=new ArrayList();
         int count []=new int[40000];
         int averagelength[]=new int[40000];
         int basenumber[]=new int[40000];
         int averagephred[]=new int[40000];
         String phred=null;         
         int location=0;
         int sum =0;
         try{
             BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
             BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);            
             while((news= br1.readLine()) != null){
                int phred1=0;
                if(!w.contains(news.split("\t")[0])){
                    w.add(news.split("\t")[0]);                    
                }
                location=w.indexOf(news.split("\t")[0]);
                count[location]=count[location]+1;
                averagelength[location]=averagelength[location]+Integer.parseInt(news.split("\t")[1]);
                phred=news.split("\t")[2];
                for(int i=0;i<phred.length();i++){
                    phred1=phred1+(int)phred.charAt(i)-33;
                }
                basenumber[location]=basenumber[location]+phred.length();
                averagephred[location]=averagephred[location]+phred1;
             }
             for(int i=0;i<w.size();i++){
                bw.write((String)w.get(i)+"\t"+count[i]+"\t"+averagelength[i]/count[i]+"\t"+averagephred[i]/basenumber[i]);
                bw.newLine();
             }
             for(int i=1;i<count.length;i++){
                 sum=sum+count[i];
             }
             System.out.println(sum);
             br1.close();
             bw.flush();bw.close();
             
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         }
     }
    public void notSDI(){
         String inputfile1="/Users/xujun/Desktop/news(MD).text";
         String inputfile2="/Users/xujun/Desktop/normal.fq";
         String outputfile="/Users/xujun/Desktop/fuck.text";
         String a=null;
         String seq=null;
         String q=null;
         String seq1=null;                   
         String q1=null;
         String startpos=null;
         String startpos1=null;
         String chro=null;
         String cigar=null;
         String MD=null;
         String [][] base=new String[13][];
         String SNPinsertion [][]=new String[13][];
         try{
            BufferedReader br1 = utils.IOUtils.getTextReader(inputfile1);
            BufferedReader br2 = utils.IOUtils.getTextReader(inputfile2);
            BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile); 
            while ((a= br1.readLine()) != null){ 
                startpos=a.split("\t")[1];
                chro=a.split("\t")[0];
                cigar=a.split("\t")[2];
                MD=a.split("\t")[3];
                br2.readLine();
                seq=br2.readLine();
                br2.readLine();
                q=br2.readLine();
                int cont=0;
                int cont2=0;
                int cont3=0;
                int cont4=0;
                boolean deletion=false;
                for(int i=0;i<cigar.length();i++){                           
                            if (!Character.isDigit(cigar.charAt(i))){ 
                                if(cigar.charAt(i)=='M'){
                                    startpos1=String.valueOf(Integer.parseInt(startpos)+cont3);
                                        seq1=seq.substring(cont2,cont2+Integer.parseInt(cigar.substring(cont,i)));  
                                        q1=q.substring(cont2,cont2+Integer.parseInt(cigar.substring(cont,i)));
                                        for(int j=0;j<MD.length();j++){
                                           
                                        }
                                        bw.write(chro+"\t");bw.write(startpos1+"\t");bw.write(seq1.length()+"\t");bw.write(q1);
                                        bw.newLine();
                                        cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
                                        cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                        cont=i+1;
                                }else{
                                    if(cigar.charAt(i)=='D'|cigar.charAt(i)=='N'){
                                        cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                        cont=i+1;
                                    }else{ 
                                        if(cigar.charAt(i)=='I'){
//                                            base[Integer.parseInt(chro)][cont4]=seq.substring(cont2+cont,cont2+i);
//                                            SNPinsertion[Integer.parseInt(chro)][cont4]=String.valueOf(Integer.parseInt(startpos)+cont3);
//                                            cont4++;
                                              cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
                                              cont=i+1;
                                        }else{
                                            cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                        }
                                        cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
                                        cont=i+1;
                                    }
                                }                                
                            }else{
                                continue;
                            }   
                }                              
            }
            br1.close();br2.close();
            bw.flush();bw.close();
         }
         catch(Exception ex){
             System.out.println(cigar);
             ex.printStackTrace();
         }
     }
    public static void main(String[] args) throws IOException, FileNotFoundException, BiffException{
        new DemoSample();
    }
    
}
