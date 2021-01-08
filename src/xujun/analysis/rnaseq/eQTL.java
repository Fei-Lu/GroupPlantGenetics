/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import com.koloboke.collect.map.hash.HashDoubleIntMap;
import com.koloboke.collect.map.hash.HashDoubleIntMaps;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import pgl.infra.range.Range;
import pgl.infra.table.ColumnTable;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import rcaller.RCaller;
import rcaller.RCode;
import xuebo.analysis.annotation.FStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/**
 *
 * @author xujun
 */
public class eQTL {
    int[] barcodeLengths = null;
    List<String> indexList = null;
    List<String>[] barcodeLists = null;

    HashMap<String, String>[] barcodeMethodMaps = null;
    //    HashMap<Integer,Integer>
    HashMap barcodeStrain = new HashMap();
    HashMap indexStrain = new HashMap();
    List<String>[] methodLists = null;

    List<String> allMethodList = new ArrayList<String>();
    String inputDirS=null;
    String outputDirS=null;
    int geneNumber=92;

    public eQTL(){//String parameterFileS
        this.filterSample();
//        this.filterSampleDS();
        this.change();
        this.subSetVCF();//根据文件当中的SNP数据，按照数目或者是比例进行抽取
//        this.countToBed();
          this.geneRegion();
        this.sort();
        this.changeVCF();
        this.fastQTL();
        this.nominalsMode();//sun fastQTL nominals mode
        this.parseBed();
        this.subSetPermutations();
        this.parseNorminals();
        this.parseHomo();
        this.parseHomoNorminals();
        this.distri();
        this.intergenicPattern();
        this.homoExpre();
        this.perseHomoCis();
        this.subFastqRatio();
        this.subFastqNum();
//        this.dataCommxand();//从文件中定量抽取一定数目的reads
        this.nominalVCF();
        this.getRawSingle();
        this.changVCF();
        this.getThreeSet();
        this.expreCis();
        this.getPairFile();//从表达的数据当中按照亚基因组的表达情况分开（标准化或不标准化）
        this.getExpreCis();
        this.getExpreCisOne();
        this.getExpreCisCat();
        this.cateExprePatten();
        this.getGeneChro();//为getGeneUpVCF服务，为了提高速度而讲这些基因分到不同的文件之下
        this.getGeneUpVCF();//得到只含有基因上游1M的snp的VCF
        this.getTajimaD();//实现多个线程同时计算tajima'D的至
        this.getTajimaDMerge();
        this.categoryCis();
        this.getGeneUpTajima();//通过对单位点进行tajima'D的计算，整合基因上游一定区域内的值
        this.getCisEffectSize();//使用GTEx的方法来就算effect size
//        this.removeBandT();
        this.batch();
    }
    public void batch (){
        String outputDirS="/data1/home/junxu/eQTL/7check/";
        String inputDirS="/data1/home/junxu/eQTL/7phenotype/sams/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.bam");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        try{
            File dir = new File(new File ("/data1/home/junxu/eQTL/7check/bam/").getAbsolutePath());
            int count=0;List<File> fL1 = new ArrayList(Arrays.asList());
            for(int i=0;i<fList.size();i++){
                if((i+1)/20==count && i!=fList.size()-1){//实现了同时跑10个HTSeq
                    fL1.add(fList.get(i));
                }else{
                    fL1.add(fList.get(i));
                    fList.parallelStream().forEach(f ->{
//                        BufferedWriter bw = IOUtils.getTextWriter(outputDirS+"/"+f.getName().split("\\.")[0]+".sh");
                        try{
                            for(int j=1;j<11;j++){
                                StringBuilder sb = new StringBuilder();
                                sb.append("java -jar /data1/home/junxu/software/picard.jar SortSam I="+f.getAbsolutePath()+" O="+f.getAbsolutePath().replace("Aligned.out.bam", "Aligned.sorted.bam"));
                                sb.append(" SORT_ORDER=coordinate");
                                String command = sb.toString();
                                try {
                                    String []cmdarry ={"/bin/bash","-c",command};
                                    Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                                    p.waitFor();
                                }
                                catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
//                            bw.flush();bw.close();
                        }
                        catch(Exception ex){
                            ex.getStackTrace();
                        }

                    });
                }
            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }

    public void getCisEffectSize(){
        GeneFeature gf = new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
        DecimalFormat decFor = new DecimalFormat("0.000");
        try{
            BufferedReader br = IOUtils.getTextReader("/data1/home/junxu/tajima/BD/1.txt");
            BufferedWriter bw = IOUtils.getTextWriter("/data1/home/junxu/tajima/BD/1.effect.txt");
            String temp=null; String geneName=null;
            HashSet pos = new HashSet();
            geneName=br.readLine().split(" ")[0];
            while((temp=br.readLine())!=null){
                if(temp.split(" ")[0]==geneName){
                    pos.add(temp.split(" ")[1].split("_")[2]);continue;
                }
                int index = gf.getGeneIndex(geneName);
                BufferedReader brExpre = IOUtils.getTextGzipReader("/data1/home/junxu/tajima/BD/7_nor_countResult.txt.gz");
                String tempExpre=null;String expre=null;
                brExpre.readLine();
                while((tempExpre=brExpre.readLine())!=null){
                    if(tempExpre.split("\t")[0].equals(geneName)){
                        expre=tempExpre;break;
                    }
                }
                brExpre.close();
                BufferedReader brVCF = IOUtils.getTextGzipReader("/data1/home/junxu/tajima/BD/geneUpVCF/"+index+".vcf.gz");
                String tempVCF=null;double u=0;double d=0;
                while((tempVCF=brVCF.readLine())!=null){
                    if(temp.startsWith("#"))continue;
                    if(pos.contains(tempVCF.split("\t")[1])){
                        for(int i=9;i<tempVCF.split("\t").length;i++){
                            if(tempVCF.split("\t")[i].contains("\\./\\."))continue;
                            if(tempVCF.split("\t")[i].split("/")[0].equals("1")){
                                u+=Double.valueOf(expre.split("\t")[i-8]);
                            }else{
                                d+=Double.valueOf(expre.split("\t")[i-8]);
                            }
                        }
                        if(u==0 || d==0){
                            bw.write(temp+" 0");//这里把纯合的位点定义为0
                        }else{
                            bw.write(temp+" "+decFor.format((u/d)*1000/1000)+"\n");
                        }
                    }
                }
                brVCF.close();
                geneName=temp.split(" ")[0];pos.clear();
            }
            br.close();
            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }

    }
    public void getGeneUpTajima(){
        HashIntIntMap[] posGeneMaps = new HashIntIntMap [45];
        for(int i=0;i<posGeneMaps.length;i++){
            posGeneMaps[i]=HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        }
        GeneFeature gf = new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
//        GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/eQTL/N344/wheat_v1.1_Lulab.gff3");
        String temp1=null;
        try{
            BufferedReader br= IOUtils.getTextReader("/data1/home/junxu/tajima/BD/1.txt");
//            BufferedReader br= IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/BD/1.txt");
            while((temp1=br.readLine())!=null){
                int chro=Integer.valueOf(temp1.split(" ")[1].split("_")[1]);
                int pos=Integer.valueOf(temp1.split(" ")[1].split("_")[2]);
                posGeneMaps[chro].put(pos,gf.getGeneIndex(temp1.split(" ")[0]));
            }
            br.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
        String inputDirS = "/data1/home/junxu/tajima/BD/nega/tajima/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".D");
        List<File> fList = Arrays.asList(fs);
        int []count = new int[1000]; double [] tajima = new double[1000];
        fList.stream().forEach(f->{
            try{
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                int index=Integer.valueOf(f.getName().split("\\.")[0]);
                String temp=null;int chro=0;int pos=0;
                if(gf.getGeneStrand(index)==1){
                    int start = gf.getGeneStart(index);
                    br.readLine();
                    while((temp=br.readLine())!=null){
                        if(temp.split("\t")[3].contains("nan"))continue;
                        chro=Integer.valueOf(temp.split("\t")[0]);
                        pos=Integer.valueOf(temp.split("\t")[1]);
                        if(posGeneMaps[chro].get(pos)==-1) continue;
                        int chunk=(start-pos)/1000;
                        count[chunk]++;
                        double t=Double.valueOf(temp.split("\t")[3]);
                        tajima[chunk]+=t;
                    }
                    br.close();
                }else{
                    int end = gf.getGeneEnd(index);
                    br.readLine();
                    while((temp=br.readLine())!=null){
                        if(temp.split("\t")[3].contains("nan"))continue;
                        chro=Integer.valueOf(temp.split("\t")[0]);
                        pos=Integer.valueOf(temp.split("\t")[1]);
                        if(posGeneMaps[chro].get(pos)==-1) continue;
                        int chunk=(pos-end)/1000;
                        count[chunk]++;
                        double t=Double.valueOf(temp.split("\t")[3]);
                        tajima[chunk]+=t;
                    }
                    br.close();
                }
            }
            catch(Exception ex){
                ex.getStackTrace();
            }
        });
        try{
            BufferedWriter bw= IOUtils.getTextWriter("/data1/home/junxu/tajima/BD/nega.cis.Tajima.D");
            BufferedWriter bw1= IOUtils.getTextWriter("/data1/home/junxu/tajima/BD/nega.cis.count.txt");
            DecimalFormat decFor = new DecimalFormat("0.000");
            for(int i=0;i<count.length;i++){
                if(tajima[i]==0){
                    bw.write(0+"\n");
                }else{
                    bw.write(decFor.format((tajima[i]/count[i])*1000/1000)+"\n");
                }
                bw1.write(count[i]+"\n");
            }
            bw.flush();bw.close();
            bw1.flush();bw1.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void categoryCis(){
        HashMap<String,String> hp = new HashMap();
        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/catCis/index.txt");
        for(int i=0;i<rt.getRowNumber();i++){
            hp.put(rt.getCell(i, 0).toString(), rt.getCell(i, 1).toString());
        }
        String inputDirS="/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/";
        HashSet A = new HashSet(); HashSet B= new HashSet(); HashSet D = new HashSet();
        String temp=null;
        try{
            BufferedReader brA= IOUtils.getTextReader(inputDirS+"0_gene.txt");
            BufferedReader brB= IOUtils.getTextReader(inputDirS+"1_gene.txt");
            BufferedReader brD= IOUtils.getTextReader(inputDirS+"2_gene.txt");
            while((temp=brA.readLine())!=null){
                A.add(temp);
            }
            while((temp=brB.readLine())!=null){
                B.add(temp);
            }
            while((temp=brD.readLine())!=null){
                D.add(temp);
            }
            brA.close();brB.close();brD.close();

            BufferedReader br= IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/Ddomi.txt");
            BufferedWriter bw= IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/catCis/Ddomi.txt");
            while((temp=br.readLine())!=null){
                StringBuilder sb = new StringBuilder();
                if(A.contains(temp.split("\t")[0])){
                    sb.append("A");
                }
                if(B.contains(temp.split("\t")[1])){
                    sb.append("B");
                }
                if(D.contains(temp.split("\t")[2])){
                    sb.append("D");
                }
                bw.write(hp.get(sb.toString())+"\n");
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }

    }
    public void getTajimaDMerge(){
        String inputDirS = "/data1/home/junxu/tajima/nega/tajima/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".D");
        List<File> fList = Arrays.asList(fs);
        double [] value = new double[1000];
        fList.stream().forEach(f->{
            try{
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp=null;int row=0;
                br.readLine();
                while((temp=br.readLine())!=null){
                    value[row]+=Double.valueOf(temp.split("\t")[3]);
                    row++;
                }
                br.close();
            }
            catch(Exception ex){
                ex.getStackTrace();
            }
        });
        try{
            BufferedWriter bw= IOUtils.getTextWriter("/data1/home/junxu/tajima/nega/tajima/merge.Tajima.D");
            DecimalFormat decFor = new DecimalFormat("0.000");
            for(int i=0;i<value.length;i++){
                bw.write(decFor.format((value[i]/fList.size())*1000/1000));bw.newLine();
            }
//            for(int i=value.length-1;i>=0;i--){
//               bw.write(decFor.format((value[i]/fList.size())*1000/1000));bw.newLine();
//            }
            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void getTajimaD(){
        String inputDirS = "/data1/home/junxu/tajima/BD/nega/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fList = Arrays.asList(fs);
        int count=0;List<File> fL1 = new ArrayList(Arrays.asList());
        for(int i=0;i<fList.size();i++){
            if((i+1)/30==count && i!=fList.size()-1){//实现了同时跑10个HTSeq
                fL1.add(fList.get(i));
            }else{
                fL1.add(fList.get(i));
                fL1.parallelStream().forEach(f -> {
                    StringBuilder sb = new StringBuilder();
                    sb.append("vcftools --gzvcf ").append(f);
                    sb.append(" --out ").append(f.getName().replace(".vcf.gz", ""));
                    sb.append(" --TajimaD 1");
//                    sb.append("sed -i 's/nan/0/g' ").append(f);
                    String command = sb.toString();
                    System.out.println(command);
                    try {
                        File dir = new File(new File ("/data1/home/junxu/tajima/BD/nega/tajima/").getAbsolutePath());
                        String []cmdarry ={"/bin/bash","-c",command};
                        Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                        p.waitFor();
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                });
                count++;
                fL1.clear();
            }
        }
    }
    public void getGeneChro(){
        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/homoGene/BD/1_gene_times.txt");
        GeneFeature gf =new GeneFeature("/Users/xujun/Desktop/eQTL/N344/wheat_v1.1_Lulab.gff3");
        String geneName=gf.getGeneName(19374);
        try{
            BufferedWriter[] bw = new BufferedWriter[45];
            for(int i =0;i<bw.length;i++){
                bw[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/BD/parseB/"+i+".txt");
                bw[i].write("num\tgene");bw[i].newLine();
            }
            for(int i=0;i<rt.getRowNumber();i++){
                int index= gf.getGeneIndex(rt.getCell(i, 1).toString());
                int chro=gf.getGeneChromosome(index);
                bw[chro].write(rt.getCellAsString(i, 0)+"\t"+rt.getCellAsString(i, 1));bw[chro].newLine();
            }
            for(int i=0;i<bw.length;i++){
                bw[i].flush();bw[i].close();
            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void getGeneUpVCF(){
        GeneFeature gf = new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
//        GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/eQTL/N344/wheat_v1.1_Lulab.gff3");
        String inputFiles="/data1/home/junxu/tajima/BD/parseB/";
        File[] fs = new File(inputFiles).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        try{
            fList.stream().forEach(f -> {
                String temp1=null;String temp2=null;
                HashIntIntMap[] posGeneMaps = new HashIntIntMap [45];
                for(int i=0;i<posGeneMaps.length;i++){
                    posGeneMaps[i]=HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
                }
                BufferedReader br1=IOUtils.getTextReader(f.getAbsolutePath());
                BufferedReader br2=IOUtils.getTextGzipReader("/data1/home/junxu/eQTL/7FastQTL/5nominals/92genotype/"+f.getName().replace("txt", "snp.maf001.mis01.recode.vcf.gz"));
                BufferedWriter bw[]= new BufferedWriter[gf.getGeneNumber()];
                try{
                    br1.readLine(); HashSet indexSet = new HashSet();
                    while((temp1=br1.readLine())!=null){
                        int index=gf.getGeneIndex(temp1.split("\t")[1]);
                        indexSet.add(index);
                        int chro=gf.getGeneChromosome(index);
                        if(gf.getGeneStrand(index)==1){
                            bw[index]=IOUtils.getTextGzipWriter("/data1/home/junxu/tajima/BD/posi/"+index+".vcf.gz");
                            int start=gf.getGeneStart(index)-1000000;
                            int end=gf.getGeneStart(index);
                            for(int i =start-1;i< end+1;i++){
                                posGeneMaps[chro].put(i, index);
                            }
                        }else{
                            bw[index]=IOUtils.getTextGzipWriter("/data1/home/junxu/tajima/BD/nega/"+index+".vcf.gz");
                            int start=gf.getGeneEnd(index);
                            int end=gf.getGeneEnd(index)+1000000;
                            for(int i =start-1;i< end+1;i++){
                                posGeneMaps[chro].put(i, index);
                            }
                        }
                    }
                    while((temp2=br2.readLine())!=null){
                        if(temp2.startsWith("#")){
                            for(int i =0;i<bw.length;i++){
                                if(indexSet.contains(i)){
                                    bw[i].write(temp2);bw[i].newLine();
                                }
                            }
                            continue;
                        }
                        int chro=Integer.valueOf(f.getName().replace(".txt", ""));
                        int pos=Integer.valueOf(temp2.split("\t")[1]);
                        int index=posGeneMaps[chro].get(pos);
                        if(index!=-1){
                            bw[index].write(temp2);bw[index].newLine();
                        }
                    }
                    for(int i=0;i<bw.length;i++){
                        if(indexSet.contains(i)){
                            bw[i].flush();bw[i].close();
                        }
                    }
                }
                catch(Exception ex){
                    ex.getStackTrace();
                }

            });
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void cateExprePatten(){
        try{
            String subG="ABD";
            BufferedReader br = IOUtils.getTextGzipReader("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/category.txt");
            BufferedReader br1 = IOUtils.getTextGzipReader("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/expre.txt");
            BufferedWriter bw [] = new BufferedWriter[7];
            bw[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/balance.txt");
            bw[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/Asup.txt");
            bw[2]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/Bsup.txt");
            bw[3]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/Dsup.txt");
            bw[4]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/Adomi.txt");
            bw[5]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/Bdomi.txt");
            bw[6]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/"+subG+"/Ddomi.txt");
            String temp=null;String temp1=null;
            br.readLine();br1.readLine();
            while((temp=br.readLine())!=null){
                temp1=br1.readLine();
                int type = Integer.valueOf(temp.split("\t")[4]);
                bw[type-1].write(temp1);bw[type-1].newLine();
            }
            br.close();br1.close();
            for(int i=0;i<bw.length;i++){
                bw[i].flush();bw[i].close();
            }
        }
        catch( Exception ex){
            ex.getStackTrace();
        }
    }
    public void getExpreCisCat(){
        HashSet A= new HashSet(); HashSet B = new HashSet(); HashSet D= new HashSet();
        String temp =null;
        try{
            BufferedReader br1 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/0_gene.txt");
            BufferedReader br2 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/1_gene.txt");
            BufferedReader br = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/2_gene.txt");
            while((temp=br1.readLine())!=null){
                A.add(temp);
            }
            br1.close();
            while((temp=br2.readLine())!=null){
                B.add(temp);
            }
            br2.close();
            while((temp=br.readLine())!=null){
                D.add(temp);
            }
            br.close();
            BufferedReader br3 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/category.txt");
            BufferedWriter[] bw = new BufferedWriter[8];
            bw[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/cis/zero.txt");
            bw[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/cis/A.txt");
            bw[2]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/cis/B.txt");
            bw[3]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/cis/D.txt");
            bw[4]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/cis/AB.txt");
            bw[5]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/cis/AD.txt");
            bw[6]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/cis/BD.txt");
            bw[7]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGeneNomPop/ABD/cis/ABD.txt");
            br3.readLine();
            while((temp=br3.readLine())!=null){
                String geneName=temp.split("\t")[0];
                String geneName1=temp.split("\t")[1];
                String geneName2=temp.split("\t")[2];
                if(A.contains(geneName)){
                    if(B.contains(geneName1)){
                        if(D.contains(geneName2)){
                            bw[7].write(temp);bw[7].newLine();
                        }else{
                            bw[4].write(temp);bw[4].newLine();
                        }
                    }else{
                        if(D.contains(geneName2)){
                            bw[5].write(temp);bw[5].newLine();
                        }else{
                            bw[1].write(temp);bw[1].newLine();
                        }
                    }
                }else{
                    if(B.contains(geneName1)){
                        if(D.contains(geneName2)){
                            bw[6].write(temp);bw[6].newLine();
                        }else{
                            bw[2].write(temp);bw[2].newLine();
                        }
                    }else{
                        if(D.contains(geneName2)){
                            bw[3].write(temp);bw[3].newLine();
                        }else{
                            bw[0].write(temp);bw[0].newLine();
                        }
                    }
                }
            }
            br3.close();
            for(int i=0;i<bw.length;i++){
                bw[i].flush();bw[i].close();
            }

        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void getExpreCisOne(){
        HashSet A= new HashSet(); HashSet B = new HashSet(); HashSet D= new HashSet();
        String temp =null;
        try{
            BufferedReader br1 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/one/0_gene.txt");
            BufferedReader br2 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/one/1_gene.txt");
            BufferedReader br = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/one/2_gene.txt");
            while((temp=br1.readLine())!=null){
                A.add(temp);
            }
            br1.close();
            while((temp=br2.readLine())!=null){
                B.add(temp);
            }
            br2.close();
            while((temp=br.readLine())!=null){
                D.add(temp);
            }
            br.close();
            BufferedReader br3 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/expre/one/0.txt");
            BufferedReader br4 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/expre/one/1.txt");
            BufferedReader br5 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/expre/one/2.txt");
            BufferedWriter[] bw = new BufferedWriter[3];BufferedWriter[] bwH = new BufferedWriter[3];
            for(int i=0;i<bw.length;i++){
                bw[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/one/cis/No/"+i+".txt");
                bwH[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/one/cis/Have/"+i+".txt");
            }
            String temp1=null;String temp2=null;
            while((temp=br3.readLine())!=null){
                String geneName=temp.split("\t")[0];
                temp1=br4.readLine();
                String geneName1=temp1.split("\t")[0];
                temp2=br5.readLine();
                String geneName2=temp2.split("\t")[0];
                if(A.contains(geneName)){
                    bwH[0].write(temp);bwH[0].newLine();
                }else{
                    bw[0].write(temp);bw[0].newLine();
                }
                if(B.contains(geneName1)){
                    bwH[1].write(temp1);bwH[1].newLine();
                }else{
                    bw[1].write(temp1);bw[1].newLine();
                }
                if(D.contains(geneName2)){
                    bwH[2].write(temp2);bwH[2].newLine();
                }else{
                    bw[2].write(temp2);bw[2].newLine();
                }
            }
            br3.close();br4.close();br5.close();
            for(int i=0;i<bw.length;i++){
                bw[i].flush();bw[i].close();
                bwH[i].flush();bwH[i].close();
            }

        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void getExpreCis(){
//        HashSet ABA= new HashSet(); HashSet ABB = new HashSet();
        HashSet ABDA= new HashSet(); HashSet ABDB = new HashSet(); HashSet ABDD= new HashSet();
        String temp =null;
        try{
            BufferedReader br1 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/0_gene.txt");
            BufferedReader br2 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/1_gene.txt");
            BufferedReader br = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/2_gene.txt");
            while((temp=br1.readLine())!=null){
                ABDA.add(temp);
            }
            br1.close();
            while((temp=br2.readLine())!=null){
                ABDB.add(temp);
            }
            br2.close();
            while((temp=br.readLine())!=null){
                ABDD.add(temp);
            }
            br.close();
            BufferedReader br3 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/expre/ABD/0.txt");
            BufferedReader br4 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/expre/ABD/1.txt");
            BufferedReader br5 = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/expre/ABD/2.txt");
//            BufferedWriter[] bw = new BufferedWriter[2];BufferedWriter[] bwA = new BufferedWriter[2];
//            BufferedWriter[] bwB = new BufferedWriter[2];BufferedWriter[] bwAB = new BufferedWriter[2];
//            for(int i=0;i<bwA.length;i++){
//                bw[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/BD/cis/zero/"+i+".txt");
//                bwA[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/BD/cis/B/"+i+".txt");
//                bwB[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/BD/cis/D/"+i+".txt");
//                bwAB[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/BD/cis/BD/"+i+".txt");
//            }
            BufferedWriter[] bw = new BufferedWriter[3];BufferedWriter[] bwA = new BufferedWriter[3];
            BufferedWriter[] bwB = new BufferedWriter[3];BufferedWriter[] bwD = new BufferedWriter[3];
            BufferedWriter[] bwAB = new BufferedWriter[3];BufferedWriter[] bwAD = new BufferedWriter[3];
            BufferedWriter[] bwBD = new BufferedWriter[3];BufferedWriter[] bwABD = new BufferedWriter[3];
            for(int i=0;i<bwA.length;i++){
                bw[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/cis/zero/"+i+".txt");
                bwA[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/cis/A/"+i+".txt");
                bwB[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/cis/B/"+i+".txt");
                bwD[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/cis/D/"+i+".txt");
                bwAB[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/cis/AB/"+i+".txt");
                bwAD[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/cis/AD/"+i+".txt");
                bwBD[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/cis/BD/"+i+".txt");
                bwABD[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/cis/ABD/"+i+".txt");
            }

            String temp1=null;String temp2=null;
            while((temp=br3.readLine())!=null){
                String geneName=temp.split("\t")[0];
                temp1=br4.readLine();
                String geneName1=temp1.split("\t")[0];
                temp2=br5.readLine();
                String geneName2=temp2.split("\t")[0];
//                if(ABA.contains(geneName)){
//                    if(ABB.contains(geneName1)){
//                        bwAB[0].write(temp);bwAB[0].newLine();
//                        bwAB[1].write(temp1);bwAB[1].newLine();
//                    }else{
//                        bwA[0].write(temp);bwA[0].newLine();
//                        bwA[1].write(temp1);bwA[1].newLine();
//                    }
//                }else{
//                    if(ABB.contains(geneName1)){
//                        bwB[0].write(temp);bwB[0].newLine();
//                        bwB[1].write(temp1);bwB[1].newLine();
//                    }else{
//                        bw[0].write(temp);bw[0].newLine();
//                        bw[1].write(temp1);bw[1].newLine();
//                    }
//                }
                if(ABDA.contains(geneName)){
                    if(ABDB.contains(geneName1)){
                        if(ABDD.contains(geneName2)){
                            bwABD[0].write(temp);bwABD[0].newLine();
                            bwABD[1].write(temp1);bwABD[1].newLine();
                            bwABD[2].write(temp2);bwABD[2].newLine();
                        }else{
                            bwAB[0].write(temp);bwAB[0].newLine();
                            bwAB[1].write(temp1);bwAB[1].newLine();
                            bwAB[2].write(temp2);bwAB[2].newLine();
                        }
                    }else{
                        if(ABDD.contains(geneName2)){
                            bwAD[0].write(temp);bwAD[0].newLine();
                            bwAD[1].write(temp1);bwAD[1].newLine();
                            bwAD[2].write(temp2);bwAD[2].newLine();
                        }else{
                            bwA[0].write(temp);bwA[0].newLine();
                            bwA[1].write(temp1);bwA[1].newLine();
                            bwA[2].write(temp2);bwA[2].newLine();
                        }

                    }
                }else{
                    if(ABDB.contains(geneName1)){
                        if(ABDD.contains(geneName2)){
                            bwBD[0].write(temp);bwBD[0].newLine();
                            bwBD[1].write(temp1);bwBD[1].newLine();
                            bwBD[2].write(temp2);bwBD[2].newLine();
                        }else{
                            bwB[0].write(temp);bwB[0].newLine();
                            bwB[1].write(temp1);bwB[1].newLine();
                            bwB[2].write(temp2);bwB[2].newLine();
                        }
                    }else{
                        if(ABDD.contains(geneName2)){
                            bwD[0].write(temp);bwD[0].newLine();
                            bwD[1].write(temp1);bwD[1].newLine();
                            bwD[2].write(temp2);bwD[2].newLine();
                        }else{
                            bw[0].write(temp);bw[0].newLine();
                            bw[1].write(temp1);bw[1].newLine();
                            bw[2].write(temp2);bw[2].newLine();
                        }
                    }
                }
            }
            br3.close();br4.close();
            for(int i=0;i<bw.length;i++){
                bw[i].flush();bw[i].close();bwA[i].flush();bwA[i].close();
                bwB[i].flush();bwB[i].close();bwD[i].flush();bwD[i].close();
                bwAB[i].flush();bwAB[i].close();bwAD[i].flush();bwAD[i].close();
                bwBD[i].flush();bwBD[i].close();bwABD[i].flush();bwABD[i].close();
            }

        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void getPairFile(){
        String ouptutDirS=null;String geneName=null;
        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/homoGene/7_homoExpre.txt");
        HashSet [] one = new HashSet[3];HashSet[] AB=new HashSet[2];HashSet[] AD=new HashSet[2];HashSet[] BD=new HashSet[2];
        List [] oneList= new ArrayList[3];List [] ABList= new ArrayList[2];List [] ADList= new ArrayList[2];List [] BDList= new ArrayList[2];
        for(int i=0;i<AB.length;i++){
            AB[i]=new HashSet();AD[i]=new HashSet();BD[i]=new HashSet();
            ABList[i]=new ArrayList();ADList[i]=new ArrayList();BDList[i]=new ArrayList();
        }
        HashSet[] ABD=new HashSet[3];List [] ABDList= new ArrayList[3];
        for(int i=0;i<ABD.length;i++){
            ABD[i]=new HashSet();one[i]= new HashSet();
            ABDList[i]=new ArrayList(); oneList[i]=new ArrayList();
        }
        for( int i=0;i<rt.getRowNumber();i++){
            if(rt.getCell(i, 3).equals("1")){
                if(rt.getCell(i, 4).equals("A")){
                    one[0].add(rt.getCell(i, 0));oneList[0].add(rt.getCell(i, 0));
                }else{
                    if(rt.getCell(i, 4).equals("B")){
                        one[1].add(rt.getCell(i, 1));oneList[1].add(rt.getCell(i, 1));
                    }else{
                        one[2].add(rt.getCell(i, 2));oneList[2].add(rt.getCell(i, 2));
                    }
                }
            }
            if(rt.getCell(i, 3).equals("2")){
                if(rt.getCell(i, 4).equals("AB")){
                    AB[0].add(rt.getCell(i, 0));ABList[0].add(rt.getCell(i, 0));
                    AB[1].add(rt.getCell(i, 1));ABList[1].add(rt.getCell(i, 1));
                }else{
                    if(rt.getCell(i, 4).equals("AD")){
                        AD[0].add(rt.getCell(i, 0));ADList[0].add(rt.getCell(i, 0));
                        AD[1].add(rt.getCell(i, 2));ADList[1].add(rt.getCell(i, 2));
                    }else{
                        BD[0].add(rt.getCell(i, 1));BDList[0].add(rt.getCell(i, 1));
                        BD[1].add(rt.getCell(i, 2));BDList[1].add(rt.getCell(i, 2));
                    }
                }
                continue;
            }
            if(rt.getCell(i, 3).equals("3")){
                ABD[0].add(rt.getCell(i, 0));ABDList[0].add(rt.getCell(i, 0));
                ABD[1].add(rt.getCell(i, 1));ABDList[1].add(rt.getCell(i, 1));
                ABD[2].add(rt.getCell(i, 2));ABDList[2].add(rt.getCell(i, 2));
            }
        }
        try{
            BufferedReader br =IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/7_nor_countResult.txt");
            BufferedWriter[] bwAB = new BufferedWriter[2];
            bwAB[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/AB/0.txt");
            bwAB[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/AB/1.txt");
            BufferedWriter[] bwAD = new BufferedWriter[2];
            bwAD[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/AD/0.txt");
            bwAD[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/AD/2.txt");
            BufferedWriter[] bwBD = new BufferedWriter[2];
            bwBD[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/BD/1.txt");
            bwBD[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/BD/2.txt");
            BufferedWriter[] bwABD =new BufferedWriter[3];
            BufferedWriter[] bwOne =new BufferedWriter[3];
            for(int i=0;i<bwABD.length;i++){
                bwABD[i] = IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/ABD/"+i+".txt");
                bwOne[i] = IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/expre/one/"+i+".txt");
            }
            String temp=null;
            List<String> [][] infoOne = new ArrayList[ABD[0].size()][3];
            for(int i=0;i<infoOne.length;i++){
                for(int j=0;j<infoOne[i].length;j++){
                    infoOne[i][j]= new ArrayList();
                }
            }
            List<String> [][] infoAB = new ArrayList[AB[0].size()][2];
            List<String> [][] infoAD = new ArrayList[AD[0].size()][2];
            List<String> [][] infoBD = new ArrayList[BD[0].size()][2];
            for(int i=0;i<infoAB.length;i++){
                for(int j=0;j<infoAB[i].length;j++){
                    infoAB[i][j]= new ArrayList();
                    infoAD[i][j]= new ArrayList();
                    infoBD[i][j]= new ArrayList();
                }
            }
            for(int i=0;i<infoAD.length;i++){
                for(int j=0;j<infoAD[i].length;j++){
                    infoAD[i][j]= new ArrayList();
                }
            }
            for(int i=0;i<infoBD.length;i++){
                for(int j=0;j<infoBD[i].length;j++){
                    infoBD[i][j]= new ArrayList();
                }
            }
            List<String> [][] infoABD = new ArrayList[ABD[0].size()][3];
            for(int i=0;i<infoABD.length;i++){
                for(int j=0;j<infoABD[i].length;j++){
                    infoABD[i][j]= new ArrayList();
                }
            }
            br.readLine();
            while((temp=br.readLine())!=null){
                geneName=temp.split("\t")[0];
                for(int j=0;j<one.length;j++){
                    if(one[j].contains(geneName)){
                        int index=oneList[j].indexOf(geneName);
                        infoOne[index][j].add(temp);
                        continue;
                    }
                }
                for(int j=0;j<AB.length;j++){
                    if(AB[j].contains(geneName)){
                        int index=ABList[j].indexOf(geneName);
                        infoAB[index][j].add(temp);
                        continue;
                    }
                }
                for(int j=0;j<AD.length;j++){
                    if(AD[j].contains(geneName)){
                        int index = ADList[j].indexOf(geneName);
                        infoAD[index][j].add(temp);
                        continue;
                    }
                }
                for(int j=0;j<BD.length;j++){
                    if(BD[j].contains(geneName)){
                        int index = BDList[j].indexOf(geneName);
                        infoBD[index][j].add(temp);
                        continue;
                    }
                }
                for(int j=0;j<ABD.length;j++){
                    if(ABD[j].contains(geneName)){
                        int index = ABDList[j].indexOf(geneName);
                        infoABD[index][j].add(temp);
                        continue;
                    }
                }
            }
            br.close();
            for(int i =0;i<infoOne.length;i++){
                bwOne[0].write(infoOne[i][0].toString().replace("[", "").replace("]", ""));bwOne[0].newLine();
                bwOne[1].write(infoOne[i][1].toString().replace("[", "").replace("]",""));bwOne[1].newLine();
                bwOne[2].write(infoOne[i][2].toString().replace("[", "").replace("]",""));bwOne[2].newLine();
            }
            for(int i =0;i<infoAB.length;i++){
                bwAB[0].write(infoAB[i][0].toString().replace("[", "").replace("]",""));bwAB[0].newLine();
                bwAB[1].write(infoAB[i][1].toString().replace("[","").replace("]",""));bwAB[1].newLine();
            }
            for(int i =0;i<infoAD.length;i++){
                bwAD[0].write(infoAD[i][0].toString().replace("[","").replace("]",""));bwAD[0].newLine();
                bwAD[1].write(infoAD[i][1].toString().replace("[","").replace("]",""));bwAD[1].newLine();
            }
            for(int i =0;i<infoBD.length;i++){
                bwBD[0].write(infoBD[i][0].toString().replace("[", "").replace("]",""));bwBD[0].newLine();
                bwBD[1].write(infoBD[i][1].toString().replace("[", "").replace("]", ""));bwBD[1].newLine();
            }
            for(int i =0;i<infoABD.length;i++){
                bwABD[0].write(infoABD[i][0].toString().replace("[","").replace("]", ""));bwABD[0].newLine();
                bwABD[1].write(infoABD[i][1].toString().replace("[", "").replace("]", ""));bwABD[1].newLine();
                bwABD[2].write(infoABD[i][2].toString().replace("[","").replace("]", ""));bwABD[2].newLine();
            }
//            List<String> tList = new ArrayList(); String [] tem = null;
//            List<String> tList1 = new ArrayList(); String [] tem1 = null;
//            List<String> tList2 = new ArrayList(); String [] tem2 = null;
//            for(int i=0;i<infoAB.length;i++){
//                StringBuilder sb = new StringBuilder();
//                tList= PStringUtils.fastSplit(infoAB[i][0].toString());
//                tem = tList.toArray(new String[tList.size()]);
//                sb.append(tem[0].replace("[", "")+"\t");
//                StringBuilder sb1 = new StringBuilder();
//                tList1= PStringUtils.fastSplit(infoAB[i][1].toString());
//                tem1 = tList1.toArray(new String[tList1.size()]);
//                sb1.append(tem1[0].replace("[", "")+"\t");
//                for(int j=1;j<tem.length;j++){
//                    if(tem[j].replace("]", "").equals("0")){
//                        if(tem1[j].replace("]", "").equals("0")){
//                            sb.append(0+"\t");
//                            sb1.append(0+"\t");
//                        }else{
//                            sb.append(0+"\t");
//                            sb1.append(1+"\t");
//                        }
//                    }else{
//                        if(tem1[j].replace("]", "").equals("0")){
//                            sb.append(1+"\t");
//                            sb1.append(0+"\t");
//                        }else{
//                            double sum= Double.valueOf(tem[j].replace("]",""))+Double.valueOf(tem1[j].replace("]", ""));
//                            DecimalFormat decFor = new DecimalFormat("0.00");
//                            sb.append(decFor.format((Double.valueOf(tem[j].replace("]", ""))/sum)*100/100)+"\t");
//                            sb1.append(decFor.format((Double.valueOf(tem1[j].replace("]", ""))/sum)*100/100)+"\t");
//                        }
//                    }
//                }
//                bwAB[0].write(sb.toString().replaceAll("\\s+$", ""));bwAB[0].newLine();
//                bwAB[1].write(sb1.toString().replaceAll("\\s+$", ""));bwAB[1].newLine();
//            }
//            for(int i=0;i<infoAD.length;i++){
//                StringBuilder sb = new StringBuilder();
//                tList= PStringUtils.fastSplit(infoAD[i][0].toString());
//                tem = tList.toArray(new String[tList.size()]);
//                sb.append(tem[0].replace("[", "")+"\t");
//                StringBuilder sb1 = new StringBuilder();
//                tList1= PStringUtils.fastSplit(infoAD[i][1].toString());
//                tem1 = tList1.toArray(new String[tList1.size()]);
//                sb1.append(tem1[0].replace("[", "")+"\t");
//                for(int j=1;j<tem.length;j++){
//                    if(tem[j].replace("]", "").equals("0") ){
//                        if(tem1[j].replace("]", "").equals("0")){
//                            sb.append(0+"\t");
//                            sb1.append(0+"\t");
//                        }else{
//                            sb.append(0+"\t");
//                            sb1.append(1+"\t");
//                        }
//                    }else{
//                        if(tem1[j].replace("]", "").equals("0")){
//                            sb.append(1+"\t");
//                            sb1.append(0+"\t");
//                        }else{
//                            double sum= Double.valueOf(tem[j].replace("]", ""))+Double.valueOf(tem1[j].replace("]", ""));
//                            DecimalFormat decFor = new DecimalFormat("0.00");
//                            sb.append(decFor.format((Double.valueOf(tem[j].replace("]", ""))/sum)*100/100)+"\t");
//                            sb1.append(decFor.format((Double.valueOf(tem1[j].replace("]", ""))/sum)*100/100)+"\t");
//                        }
//                    }
//                }
//                bwAD[0].write(sb.toString().replaceAll("\\s+$", ""));bwAD[0].newLine();
//                bwAD[1].write(sb1.toString().replaceAll("\\s+$", ""));bwAD[1].newLine();
//            }
//            for(int i=0;i<infoBD.length;i++){
//                StringBuilder sb = new StringBuilder();
//                tList= PStringUtils.fastSplit(infoBD[i][0].toString());
//                tem = tList.toArray(new String[tList.size()]);
//                sb.append(tem[0].replace("[", "")+"\t");
//                StringBuilder sb1 = new StringBuilder();
//                tList1= PStringUtils.fastSplit(infoBD[i][1].toString());
//                tem1 = tList1.toArray(new String[tList1.size()]);
//                sb1.append(tem1[0].replace("[", "")+"\t");
//                for(int j=1;j<tem.length;j++){
//                    if(tem[j].replace("]", "").equals("0")){
//                        if(tem1[j].replace("]", "").equals("0")){
//                            sb.append(0+"\t");
//                            sb1.append(0+"\t");
//                        }else{
//                            sb.append(0+"\t");
//                            sb1.append(1+"\t");
//                        }
//                    }else{
//                        if(tem1[j].replace("]", "").equals("0")){
//                            sb.append(1+"\t");
//                            sb1.append(0+"\t");
//                        }else{
//                            double sum= Double.valueOf(tem[j].replace("]", ""))+Double.valueOf(tem1[j].replace("]",""));
//                            DecimalFormat decFor = new DecimalFormat("0.00");
//                            sb.append(decFor.format((Double.valueOf(tem[j].replace("]",""))/sum)*100/100)+"\t");
//                            sb1.append(decFor.format((Double.valueOf(tem1[j].replace("]",""))/sum)*100/100)+"\t");
//                        }
//                    }
//                }
//                bwBD[0].write(sb.toString().replaceAll("\\s+$", ""));bwBD[0].newLine();
//                bwBD[1].write(sb1.toString().replaceAll("\\s+$", ""));bwBD[1].newLine();
//            }
//            for(int i=0;i<infoABD.length;i++){
//                StringBuilder sb = new StringBuilder();
//                tList= PStringUtils.fastSplit(infoABD[i][0].toString());
//                tem = tList.toArray(new String[tList.size()]);
//                sb.append(tem[0].replace("[", "")+"\t");
//                StringBuilder sb1 = new StringBuilder();
//                tList1= PStringUtils.fastSplit(infoABD[i][1].toString());
//                tem1 = tList1.toArray(new String[tList1.size()]);
//                sb1.append(tem1[0].replace("[", "")+"\t");
//                StringBuilder sb2 = new StringBuilder();
//                tList2= PStringUtils.fastSplit(infoABD[i][2].toString());
//                tem2 = tList2.toArray(new String[tList2.size()]);
//                sb2.append(tem2[0].replace("[", "")+"\t");
//                for(int j=1;j<tem.length;j++){
//                    if(tem[j].replace("]","").equals("0")){
//                        if(tem1[j].replace("]", "").equals("0") ){
//                            if(tem2[j].replace("]", "").equals("0") ){
//                                sb.append(0+"\t");
//                                sb1.append(0+"\t");
//                                sb2.append(0+"\t");
//                            }else{
//                                sb.append(0+"\t");
//                                sb1.append(0+"\t");
//                                sb2.append(1+"\t");
//                            }
//                        }else{
//                            if(tem2[j].replace("]", "").equals("0") ){
//                                sb.append(0+"\t");
//                                sb1.append(1+"\t");
//                                sb2.append(0+"\t");
//                            }else{
//                                double sum= Double.valueOf(tem1[j].replace("]", ""))+Double.valueOf(tem2[j].replace("]",""));
//                                DecimalFormat decFor = new DecimalFormat("0.00");
//                                sb.append(0+"\t");
//                                sb1.append(decFor.format((Double.valueOf(tem1[j].replace("]",""))/sum)*100/100)+"\t");
//                                sb2.append(decFor.format((Double.valueOf(tem2[j].replace("]",""))/sum)*100/100)+"\t");
//                            }
//                        }
//                    }else{
//                        if(tem1[j].replace("]", "").equals("0") ){
//                            if(tem2[j].replace("]", "").equals("0") ){
//                                sb.append(1+"\t");
//                                sb1.append(0+"\t");
//                                sb2.append(0+"\t");
//                            }else{
//                                double sum= Double.valueOf(tem[j].replace("]", ""))+Double.valueOf(tem2[j].replace("]",""));
//                                DecimalFormat decFor = new DecimalFormat("0.00");
//                                sb1.append(0+"\t");
//                                sb.append(decFor.format((Double.valueOf(tem[j].replace("]",""))/sum)*100/100)+"\t");
//                                sb2.append(decFor.format((Double.valueOf(tem2[j].replace("]",""))/sum)*100/100)+"\t");
//                            }
//                        }else{
//                            if(tem2[j].replace("]", "").equals("0") ){
//                                double sum= Double.valueOf(tem[j].replace("]", ""))+Double.valueOf(tem1[j].replace("]",""));
//                                DecimalFormat decFor = new DecimalFormat("0.00");
//                                sb2.append(0+"\t");
//                                sb.append(decFor.format((Double.valueOf(tem[j].replace("]",""))/sum)*100/100)+"\t");
//                                sb1.append(decFor.format((Double.valueOf(tem1[j].replace("]",""))/sum)*100/100)+"\t");
//                            }else{
//                                double sum= Double.valueOf(tem[j].replace("]", ""))+Double.valueOf(tem1[j].replace("]", ""))+Double.valueOf(tem2[j].replace("]", ""));
//                                DecimalFormat decFor = new DecimalFormat("0.00");
//                                sb.append(decFor.format((Double.valueOf(tem[j].replace("]", ""))/sum)*100/100)+"\t");
//                                sb1.append(decFor.format((Double.valueOf(tem1[j].replace("]", ""))/sum)*100/100)+"\t");
//                                sb2.append(decFor.format((Double.valueOf(tem2[j].replace("]", ""))/sum)*100/100)+"\t");
//                            }
//                        }
//                    }
//                }
//                bwABD[0].write(sb.toString().replaceAll("\\s+$", ""));bwABD[0].newLine();
//                bwABD[1].write(sb1.toString().replaceAll("\\s+$", ""));bwABD[1].newLine();
//                bwABD[2].write(sb2.toString().replaceAll("\\s+$", ""));bwABD[2].newLine();
//            }
            for(int i=0;i<bwOne.length;i++){
                bwOne[i].flush();bwOne[i].close();
            }
            for(int i=0;i<bwAB.length;i++){
                bwAB[i].flush();bwAB[i].close();
            }
            for(int i=0;i<bwAD.length;i++){
                bwAD[i].flush();bwAD[i].close();
            }
            for(int i=0;i<bwBD.length;i++){
                bwBD[i].flush();bwBD[i].close();
            }
            for(int i=0;i<bwABD.length;i++){
                bwABD[i].flush();bwABD[i].close();
            }
        }
        catch(Exception ex){
            System.out.println(geneName);
            ex.getStackTrace();
        }
    }
    public void expreCis(){
        try{
            BufferedReader brA = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/0_gene.txt");
            BufferedReader brB = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/1_gene.txt");
//             BufferedReader brD = IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/2_gene.txt");
            HashSet A=new HashSet();HashSet B =new HashSet(); HashSet D = new HashSet();
            String temp=null;
            while((temp=brA.readLine())!=null){
                A.add(temp.replaceAll("\t", ""));
            }
            brA.close();
            while((temp=brB.readLine())!=null){
                B.add(temp.replaceAll("\t", ""));
            }
            brB.close();
//             while((temp=brD.readLine())!=null){
//                 D.add(temp.replaceAll("\t", ""));
//             }
//             brD.close();
            BufferedWriter bw = IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/expre_cis.txt");
            RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/homoGene/7_homoExpre.txt");
            for(int i=0;i< rt.getRowNumber();i++){
                if(rt.getCell(i, 4).equals("AB")){//changeable
                    int count=0; String sub="";
                    StringBuilder sb = new StringBuilder();
                    if(A.contains(rt.getCell(i, 0))){
                        count++;sub+="A";
                    }
                    sb.append(rt.getCell(i, 0)+"\t");
                    if(B.contains(rt.getCell(i, 1))){
                        count++; sub+="B";
                    }
                    sb.append(rt.getCell(i, 1)+"\t");
//                     if(D.contains(rt.getCell(i, 2))){
//                         count++; sub+="D";
//                     }
//                     sb.append(rt.getCell(i, 2)+"\t");
                    bw.write(sb.toString()+"0"+"\t"+count+"\t"+sub);
                    bw.newLine();
                }
            }
            bw.flush();bw.close();

        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void getThreeSet(){
//        RowTable rt = new RowTable("/data1/home/junxu/eQTL/7FastQTL/5nominals/getSub7.txt");
        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/getSub7.txt");
        HashMap chroChro= new HashMap();
        for(int i=0;i<rt.getRowNumber();i++){
            chroChro.put(rt.getCell(i, 0), rt.getCell(i, 1));
        }
//        String inputFiles="/data1/home/junxu/eQTL/7FastQTL/5nominals/plink/test.ld";
//        String outputFiles1="/data1/home/junxu/eQTL/7FastQTL/5nominals/within.ld.txt";
//        String outputFiles2="/data1/home/junxu/eQTL/7FastQTL/5nominals/across.ld.txt";
//        String outputFiles3="/data1/home/junxu/eQTL/7FastQTL/5nominals/homo.ld.txt";
        String inputFiles="/Users/xujun/Desktop/eQTL/N344/test.ld";
        String outputFiles1="/Users/xujun/Desktop/eQTL/N344/within.ld.txt";
        String outputFiles2="/Users/xujun/Desktop/eQTL/N344/across.ld.txt";
        String outputFiles3="/Users/xujun/Desktop/eQTL/N344/homo.ld.txt";
        try{
            String temp=null;
//            BufferedReader brHomo=IOUtils.getTextReader("/data1/home/junxu/eQTL/7FastQTL/5nominals/homoGene.txt");
            BufferedReader brHomo=IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/homoGene.txt");
            HashMap AB=new HashMap();HashMap AD= new HashMap();HashMap BD= new HashMap();
            while((temp=brHomo.readLine())!=null){
                AB.put(temp.split("\t")[0], temp.split("\t")[1]);
                AD.put(temp.split("\t")[0], temp.split("\t")[2]);
                BD.put(temp.split("\t")[1], temp.split("\t")[2]);
            }
            brHomo.close();
//            BufferedReader brNom=IOUtils.getTextGzipReader("/data1/home/junxu/eQTL/7FastQTL/5nominals/homoGene.nominals.all.txt.gz");
            BufferedReader brNom=IOUtils.getTextGzipReader("/Users/xujun/Desktop/eQTL/N344/homoGene.nominals.all.txt.gz");
            HashMap[] SNPGene = new HashMap[45];
            for(int i=0;i<SNPGene.length;i++){
                SNPGene[i]=new HashMap();
            }
            while((temp=brNom.readLine())!=null){
                int chro=Integer.valueOf(temp.split(" ")[1].split("_")[1]);
                SNPGene[chro].put(temp.split(" ")[1].split("_")[2], temp.split(" ")[0]);
            }
            brNom.close();
            BufferedReader br=IOUtils.getTextReader(inputFiles);
            BufferedWriter bw1 = IOUtils.getTextWriter(outputFiles1);
            BufferedWriter bw2 = IOUtils.getTextWriter(outputFiles2);
            BufferedWriter bw3 = IOUtils.getTextWriter(outputFiles3);
            br.readLine();
            while((temp=br.readLine())!=null){
                String chro1=chroChro.get(temp.split(" ")[1]).toString();
                String chro2=chroChro.get(temp.split(" ")[4]).toString();
                if(chro1.substring(1).equals(chro2.substring(1))){
                    bw1.write(temp);bw1.newLine();
                }else{
                    if(!chro1.substring(0, 1).equals(chro2.substring(0,1))){
                        bw2.write(temp);bw2.newLine();
                    }else{
                        int c1=Integer.valueOf(temp.split(" ")[1]);
                        String gene1=SNPGene[c1].get(temp.split(" ")[2]).toString();
                        int c2=Integer.valueOf(temp.split(" ")[4]);
                        String gene2=SNPGene[c2].get(temp.split(" ")[5]).toString();
                        if(chro1.substring(1).equals("A") && chro2.substring(1).equals("B")){
                            if(AB.get(gene1).equals(gene2)){
                                bw3.write(temp);bw3.newLine();
                                continue;
                            }
                        }
                        if(chro1.substring(1).equals("A") && chro2.substring(1).equals("D")){
                            if(AD.get(gene1).equals(gene2)){
                                bw3.write(temp);bw3.newLine();
                                continue;
                            }
                        }
                        if(chro1.substring(1).equals("B") && chro2.substring(1).equals("D")){
                            if(BD.get(gene1).equals(gene2)){
                                bw3.write(temp);bw3.newLine();
                                continue;
                            }
                        }
                    }
                }
            }
            br.close();
            bw1.flush();bw1.close();bw2.flush();bw2.close();bw3.flush();bw3.close();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    public void changVCF(){
        String inputFiles="/data1/home/junxu/eQTL/7FastQTL/5nominals/92HomoGeneVCF/";
        String outputFiles="/data1/home/junxu/eQTL/7FastQTL/5nominals/92HomoGeneVCF/";
        File[] fs = new File(inputFiles).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "recode.vcf");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        fList.stream().forEach(f -> {
            try{
                BufferedReader br=IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputFiles,f.getName().replace("recode.vcf","92.change.vcf")).getAbsolutePath());
//             BufferedReader br=IOUtils.getTextGzipReader(inputFiles);
//             BufferedWriter bw = IOUtils.getTextGzipWriter(outputFiles);
                String chro=f.getName().split("\\.")[0];
                String temp=null; String [] tem = null;List<String> tList=new ArrayList();
                StringBuilder header = new StringBuilder();
                HashSet colSet=new HashSet();
                while((temp=br.readLine())!=null){
                    if(temp.startsWith("##")){
                        bw.write(temp);bw.newLine();
                        continue;
                    }
                    if(temp.startsWith("#")){
                        StringBuilder out =new StringBuilder();
                        tList= PStringUtils.fastSplit(temp);
                        tem = tList.toArray(new String[tList.size()]);
                        for(int i=0;i<tem.length;i++){
                            if(i<9){
                                out.append(tem[i]+"\t");
                            }else{
                                out.append(chro+tem[i]+"\t");
                            }
                        }
                        bw.write(out.toString().replaceAll("\\s+$", ""));
                        bw.newLine();
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        });
    }
    public void getRawSingle(){
        String inputDirS="/data2/junxu/SiPASData/three-lanes-firsttime-copydata0129/data_P101SC18112845-01-F004-B4-21/1.rawdata/S20190118S1-7-20190118S1-7_FKDL190665547-1a/";
        String outputDirS ="/data1/home/junxu/eQTL/phenotype/";
//         RowTable rt = new RowTable("/data2/junxu/SiPASData/index.txt");
//         HashMap indexName = new HashMap();
//         List<String> nameList = new ArrayList();
//         for (int i=0;i<rt.getRowNumber();i++){
//             indexName.put(rt.getCellAsString(i, 1), rt.getCellAsString(i, 0));
//             nameList.add(rt.getCellAsString(i, 0));
//         }
//        File[] subDirS = new File(inputDirS).listFiles();
//        HashSet<String> nameSet = new HashSet();
//        for(int i =0;i<subDirS.length;i++){
//            if(subDirS[i].isDirectory()){
//                File[] fs = subDirS[i].listFiles();
//                fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
//                for(int j=0;j<fs.length;j++){
//                    if (fs[j].isHidden()) continue;
//                    nameSet.add(fs[j].getName().replace("_"+fs[j].getName().split("_")[fs[j].getName().split("_").length-1],""));
//                }
//            }
//        }
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().replace("_"+fs[i].getName().split("_")[fs[i].getName().split("_").length-1],""));
        }
//        List<File> fList = Arrays.asList(fs);
//        BufferedWriter [] bw = new BufferedWriter[nameList.size()];
//        BufferedWriter [] bw1 = new BufferedWriter[nameList.size()];
//        for(int i=0;i<nameList.size();i++){
//            bw[i]=IOUtils.getTextGzipWriter(outputDirS+"/"+nameList.get(i)+"_R1.fq.gz");
//            bw1[i]=IOUtils.getTextGzipWriter(outputDirS+"/"+nameList.get(i)+"_R2.fq.gz");
//        }
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS,p+"_1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"_2.fq.gz").getAbsolutePath();
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw=IOUtils.getTextGzipWriter(outputDirS+"/"+p+"_R1.fq.gz");
                BufferedWriter bw1=IOUtils.getTextGzipWriter(outputDirS+"/"+p+"_R2.fq.gz");
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    String index=temp.split(":")[9];
                    if(index.equals("TCCCGA")){
                        bw.write(temp);bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                        bw1.write(br1.readLine());bw1.newLine();
                        bw1.write(br1.readLine());bw1.newLine();
                        bw1.write(br1.readLine());bw1.newLine();
                        bw1.write(br1.readLine());bw1.newLine();
                    }else{
                        br.readLine();br.readLine();br.readLine();
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    }
                }
                br.close();br1.close();
                bw.flush();bw.close();
                bw1.flush();bw1.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    public void getRaw(){
        String inputDirS="/data2/junxu/SiPASData/three-lanes-firsttime-copydata0129/data_P101SC18112845-01-F004-B4-21/1.rawdata/";
        String outputDirS ="/data2/junxu/SiPASData/20190319/Rawdata/";
        RowTable rt = new RowTable("/data2/junxu/SiPASData/index.txt");
        HashMap indexName = new HashMap();
        List<String> nameList = new ArrayList();
        for (int i=0;i<rt.getRowNumber();i++){
            indexName.put(rt.getCellAsString(i, 1), rt.getCellAsString(i, 0));
            nameList.add(rt.getCellAsString(i, 0));
        }
        File[] subDirS = new File(inputDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for(int i =0;i<subDirS.length;i++){
            if(subDirS[i].isDirectory()){
                File[] fs = subDirS[i].listFiles();
                fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
                for(int j=0;j<fs.length;j++){
                    if (fs[j].isHidden()) continue;
                    nameSet.add(fs[j].getName().replace("_"+fs[j].getName().split("_")[fs[j].getName().split("_").length-1],""));
                }
            }
        }
//        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
//        HashSet<String> nameSet = new HashSet();
//        for (int i = 0; i < fs.length; i++) {
//            if (fs[i].isHidden()) continue;
//            nameSet.add(fs[i].getName().split("_")[0]);
//        }
//        List<File> fList = Arrays.asList(fs);
        BufferedWriter [] bw = new BufferedWriter[nameList.size()];
        BufferedWriter [] bw1 = new BufferedWriter[nameList.size()];
        for(int i=0;i<nameList.size();i++){
            bw[i]=IOUtils.getTextGzipWriter(outputDirS+"/"+nameList.get(i)+"_R1.fq.gz");
            bw1[i]=IOUtils.getTextGzipWriter(outputDirS+"/"+nameList.get(i)+"_R2.fq.gz");
        }
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS, p+"/"+p+"_1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"/"+p+"_2.fq.gz").getAbsolutePath();
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    String index=temp.split(":")[9];
                    if(indexName.get(index)!=null){
                        int num =nameList.indexOf(indexName.get(index));
                        bw[num].write(temp);bw[num].newLine();
                        bw[num].write(br.readLine());bw[num].newLine();
                        bw[num].write(br.readLine());bw[num].newLine();
                        bw[num].write(br.readLine());bw[num].newLine();
                        bw1[num].write(br1.readLine());bw1[num].newLine();
                        bw1[num].write(br1.readLine());bw1[num].newLine();
                        bw1[num].write(br1.readLine());bw1[num].newLine();
                        bw1[num].write(br1.readLine());bw1[num].newLine();
                    }else{
                        br.readLine();br.readLine();br.readLine();
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    }
                }
                br.close();br1.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try{
            for (int i=0;i<nameList.size();i++){
                bw[i].flush();bw[i].close();
                bw1[i].flush();bw1[i].close();
            }
        }
        catch(Exception ex){
            ex.printStackTrace();
        }
    }
    public void nominalVCF(){
        HashSet[] hs = new HashSet[45];
        for(int i=0;i<hs.length;i++){
            hs[i]=new HashSet();
        }
        String temp=null;
//        String inputFileDirS="/data2/junxu/SNP/";
        String inputFileDirS="/data1/home/junxu/eQTL/7FastQTL/5nominals/";
        File[] fs = new File(inputFileDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "snp.maf001.mis01.recode.vcf.gz");
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        System.out.println();
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        try{
            BufferedReader br =IOUtils.getTextGzipReader("/data1/home/junxu/eQTL/7FastQTL/5nominals/001homoGene.nominals.all.txt.gz");
            while((temp=br.readLine())!=null){
                int chr=Integer.valueOf(temp.split(" ")[1].split("_")[1]);
                hs[chr].add( temp.split(" ")[1].split("_")[2]);
            }
            br.close();
//            BufferedWriter bw =IOUtils.getTextWriter("/data1/home/junxu/eQTL/7FastQTL/5nominals/homoGene.vcf");
            fList.stream().forEach(f -> {
                try{
                    BufferedReader vcf =IOUtils.getTextReader(f.getAbsolutePath());
                    BufferedWriter bw =IOUtils.getTextWriter("/data1/home/junxu/eQTL/7FastQTL/5nominals/001"+f.getName());
                    String t=null;
                    while((t=vcf.readLine())!=null){
                        char s = t.charAt(0);
                        if ((int)s < 48 || (int)s > 57) {
                            bw.write(t);bw.newLine();
                            continue;
                        }
                        int chr=Integer.valueOf(t.split("\t")[0]);
                        if(hs[chr].contains(t.split("\t")[1])){
                            bw.write(t);bw.newLine();
                        }
                    }
                    vcf.close();
                    bw.flush();bw.close();
                }
                catch(Exception ex){
                    ex.getStackTrace();
                }
            });
//            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void dataCommand() {
        String subFqDirS = new File ("/data1/home/junxu/wheat/SiPAS1.1/test/subFastqs/").getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, "R1.fq.gz");
        for(int i=0;i<fs.length;i++){
//            if (fs[i].length() < 450000000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            if (fs[i].length() < 159000000) continue;
            fList.add(fs[i]);
        }
        try{
            BufferedWriter bw =IOUtils.getTextWriter("/data1/home/junxu/wheat/SiPAS1.1/test/2M/file.txt");
            for(int i=0;i<fList.size();i++){
                bw.write(fList.get(i).getName().replace("_R1.fq.gz",""));
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.printStackTrace();
        }
    }
    public void subFastqRatio(){
        String inputDirS="/data2/junxu/SiPASResult/SiPASReverse/subFastqs/";
        String outputDirS="/data2/junxu/SiPASResult/SiPASReverse/2G/subFastqs/";
//        String inputDirS="/Users/xujun/Desktop/eQTL/N344/";
//        String outputDirS="/data2/junxu/SiPASResult/SiPASReverse/2G/subFastqs/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq.gz");
//        fs = IOUtils.listFilesStartsWith(fs, "homoGene.nominals");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        RowTable rtAm = new RowTable("/data2/junxu/SiPASResult/SiPASReverse/subFastqs/countam.txt");
        HashMap hmAm = new HashMap ();
        for(int i=0;i<rtAm.getRowNumber();i++){
            hmAm.put("am"+(i+1), rtAm.getCell(i, 0));
        }
        RowTable rtPm = new RowTable("/data2/junxu/SiPASResult/SiPASReverse/subFastqs/countpm.txt");
        HashMap hmPm = new HashMap ();
        for(int i=0;i<rtPm.getRowNumber();i++){
            hmPm.put("pm"+(i+1), rtAm.getCell(i, 0));
        }
        try{
            int count=0;List<File> fL1 = new ArrayList(Arrays.asList());
            for(int i=0;i<fList.size();i++){
                if((i+1)/2==count && i!=fList.size()-1){//实现了同时跑10个HTSeq
                    fL1.add(fList.get(i));
                }else{
                    fL1.add(fList.get(i));
                    fL1.parallelStream().forEach(f -> {
                        try{
                            int num=0;
                            BufferedWriter bwR1 =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName()).getAbsolutePath());
                            if(hmAm.get(f.getName().split("_")[0].split("-")[1])!=null){
                                num=Integer.valueOf(hmAm.get(f.getName().split("_")[0].split("-")[1]).toString());
                            }else{
                                num=Integer.valueOf(hmPm.get(f.getName().split("_")[0].split("-")[1]).toString());
                            }
                            int fqNumber=num/4;
                            List<Integer> total =new ArrayList();
                            for(int a=1;a<fqNumber+1;a++){
                                total.add(a);
                            }
                            Collections.shuffle(total);
                            int sub = (int)Math.round(fqNumber*0.4);
                            System.out.println(f.getName()+"\t"+fqNumber+"\t"+sub);
                            List subList = total.subList(0, sub+1);
                            String tempR1=null; int rowR1=1;
                            BufferedReader brR1 =IOUtils.getTextGzipReader(f.getAbsolutePath());
                            while((tempR1=brR1.readLine())!=null){
                                if(subList.contains(rowR1)){
                                    bwR1.write(tempR1);bwR1.newLine();
                                    bwR1.write(brR1.readLine());bwR1.newLine();
                                    bwR1.write(brR1.readLine());bwR1.newLine();
                                    bwR1.write(brR1.readLine());bwR1.newLine();
                                }
                                rowR1++;
                            }
                            brR1.close();
                            bwR1.flush();bwR1.close();

                            BufferedWriter bwR2 =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName().replace("_R1", "_R2")).getAbsolutePath());
                            String tempR2=null; int rowR2=1;
                            BufferedReader brR2 =IOUtils.getTextGzipReader(f.getAbsolutePath().replace("_R1.fq.gz", "_R2.fq.gz"));
                            while((tempR2=brR2.readLine())!=null){
                                if(subList.contains(rowR2)){
                                    bwR2.write(tempR2);bwR1.newLine();
                                    bwR2.write(brR2.readLine());bwR2.newLine();
                                    bwR2.write(brR2.readLine());bwR2.newLine();
                                    bwR2.write(brR2.readLine());bwR2.newLine();
                                }
                                rowR2++;
                            }
                            brR2.close();
                            bwR2.flush();bwR2.close();
                        }
                        catch(Exception ex){
                            ex.getStackTrace();
                        }
                    });
                    count++;
                    fL1.clear();
                }
            }

        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void subFastqNum(){
//        String inputDirS="/data2/junxu/SiPASResult/SiPASAm2Pm1/difVolumn/";
//        String inputDirS="/data2/junxu/SiPASResult/SiPASUR/difVolumn/";
//        String inputDirS="/data2/junxu/SiPASResult/SiPASReverse/2G/difVolumn/";
        String inputDirS="/data1/home/junxu/wheat/SiPAS1.1/test/";
        String outputDirS="/data1/home/junxu/wheat/SiPAS1.1/testR/repTest/subFastqs";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq.gz");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
//        RowTable rt = new RowTable("/data2/junxu/SiPASResult/SiPASReverse/2G/difVolumn/count1.txt");
//        RowTable rt = new RowTable("/data2/junxu/SiPASResult/SiPASAm2Pm1/difVolumn/count.txt");
        RowTable rt = new RowTable("/data1/home/junxu/wheat/SiPAS1.1/test/count.txt");
//        RowTable rt = new RowTable("/data2/junxu/SiPASResult/SiPASUR/difVolumn/count.txt");
        HashMap hm = new HashMap ();
        for(int i=0;i<rt.getRowNumber();i++){
            hm.put(rt.getCell(i, 0), rt.getCell(i, 1));
        }
        try{
            int count=0;List<File> fL1 = new ArrayList(Arrays.asList());
            List<Integer>total =new ArrayList();
            for(int i=0;i<fList.size();i++){
                if((i+1)/1==count && i!=fList.size()-1){//实现了同时跑10个HTSeq
                    fL1.add(fList.get(i));
                }else{
                    fL1.add(fList.get(i));
                    fL1.parallelStream().forEach(f -> {
                        try{
                            int num=0;
                            if(hm.get(f.getName().split("_")[0])!=null){
                                num=Integer.valueOf(hm.get(f.getName().split("_")[0]).toString());
                            }
                            int fqNumber=num/4;
                            for(int a=1;a<fqNumber+1;a++){
                                total.add(a);
                            }
                            for(int rep=1;rep<4;rep++){
                                Collections.shuffle(total);
                                int sub = 2000000;
                                List<Integer> subList = total.subList(0, sub);
                                Set<Integer> set = new HashSet<Integer>(subList);
                                BufferedWriter bwR1 =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName().split("_")[0]+rep+"_R1.fq.gz").getAbsolutePath());
                                String tempR1=null; int rowR1=1;
                                BufferedReader brR1 =IOUtils.getTextGzipReader(f.getAbsolutePath());
                                while((tempR1=brR1.readLine())!=null){
                                    if(set.contains(rowR1)){
                                        bwR1.write(tempR1);bwR1.newLine();
                                        bwR1.write(brR1.readLine());bwR1.newLine();
                                        bwR1.write(brR1.readLine());bwR1.newLine();
                                        bwR1.write(brR1.readLine());bwR1.newLine();
                                    }else{
                                        brR1.readLine();brR1.readLine();brR1.readLine();
                                    }
                                    rowR1++;
                                }
                                brR1.close();
                                bwR1.flush();bwR1.close();

                                BufferedWriter bwR2 =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName().split("_")[0]+rep+"_R2.fq.gz").getAbsolutePath());
                                String tempR2=null; int rowR2=1;
                                BufferedReader brR2 =IOUtils.getTextGzipReader(f.getAbsolutePath().replace("_R1.fq.gz", "_R2.fq.gz"));
                                while((tempR2=brR2.readLine())!=null){
                                    if(set.contains(rowR2)){
                                        bwR2.write(tempR2);bwR2.newLine();
                                        bwR2.write(brR2.readLine());bwR2.newLine();
                                        bwR2.write(brR2.readLine());bwR2.newLine();
                                        bwR2.write(brR2.readLine());bwR2.newLine();
                                    }else{
                                        brR2.readLine();brR2.readLine();brR2.readLine();
                                    }
                                    rowR2++;
                                }
                                brR2.close();
                                bwR2.flush();bwR2.close();
                            }
                        }
                        catch(Exception ex){
                            ex.getStackTrace();
                        }
                    });
                    count++;
                    fL1.clear();
                }
            }

        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void perseHomoCis(){
        String ouptutDirS=null;
        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/expre_cis.txt");
        HashSet[] zero=new HashSet[3];HashSet[] one=new HashSet[3];
        HashSet[] AB=new HashSet[2];HashSet[] AD=new HashSet[2];HashSet[] BD=new HashSet[2];
        for(int i=0;i<AB.length;i++){
            AB[i]=new HashSet();AD[i]=new HashSet();BD[i]=new HashSet();
        }
        HashSet[] ABD=new HashSet[3];
        for(int i=0;i<one.length;i++){
            zero[i]=new HashSet();
            one[i]=new HashSet();
            ABD[i]=new HashSet();
        }
        for( int i=0;i<rt.getRowNumber();i++){
            if(rt.getCell(i, 3).equals("0")){
                zero[0].add(rt.getCell(i, 0));
                zero[1].add(rt.getCell(i, 1));
                zero[2].add(rt.getCell(i, 2));
                continue;
            }
            if(rt.getCell(i, 3).equals("1")){
                if(rt.getCell(i, 4).equals("A")){
                    one[0].add(rt.getCell(i, 0));
                }else{
                    if(rt.getCell(i, 4).equals("B")){
                        one[1].add(rt.getCell(i, 1));
                    }else{
                        one[2].add(rt.getCell(i, 2));
                    }
                }
                continue;
            }
            if(rt.getCell(i, 3).equals("2")){
                if(rt.getCell(i, 4).equals("AB")){
                    AB[0].add(rt.getCell(i, 0));
                    AB[1].add(rt.getCell(i, 1));
                }else{
                    if(rt.getCell(i, 4).equals("AD")){
                        AD[0].add(rt.getCell(i, 0));
                        AD[1].add(rt.getCell(i, 2));
                    }else{
                        BD[0].add(rt.getCell(i, 1));
                        BD[1].add(rt.getCell(i, 2));
                    }
                }
                continue;
            }
            if(rt.getCell(i, 3).equals("3")){
                ABD[0].add(rt.getCell(i, 0));
                ABD[1].add(rt.getCell(i, 1));
                ABD[2].add(rt.getCell(i, 2));
            }
        }
        try{
//            BufferedReader br =IOUtils.getTextGzipReader("/Users/xujun/Desktop/eQTL/N344/nominals.all.txt.gz");
            BufferedReader br =IOUtils.getTextReader("/Users/xujun/Desktop/eQTL/N344/7_nor_countResult.txt");
            BufferedWriter[] bwzero = new BufferedWriter[3];
            BufferedWriter[] bwone = new BufferedWriter[3];
            BufferedWriter[] bwAB = new BufferedWriter[2];
            bwAB[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/expre_cis/AB/0.txt");
            bwAB[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/expre_cis/AB/1.txt");
//            BufferedWriter[] bwAD = new BufferedWriter[2];
//            bwAD[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/expre_cis/AD/0.txt");
//            bwAD[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/expre_cis/AD/2.txt");
//            BufferedWriter[] bwBD = new BufferedWriter[2];
//            bwBD[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/expre_cis/BD/1.txt");
//            bwBD[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/expre_cis/BD/2.txt");
//            BufferedWriter[] bwABD =new BufferedWriter[3];
            for(int i=0;i<bwone.length;i++){
                bwzero[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/expre_cis/zero/"+i+".txt");
                bwone[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/expre_cis/one/"+i+".txt");
//                bwABD[i] = IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/expre_cis/ABD/"+i+".txt");
            }
            String temp=null;String geneName=null;
            List<String> tList = new ArrayList(); String [] tem = null;
            while((temp=br.readLine())!=null){
                geneName=temp.split("\t")[0];
                for(int j=0;j<zero.length;j++){
                    if(zero[j].contains(geneName)){
                        bwzero[j].write(temp);bwzero[j].newLine();
                        continue;
                    }
                }
                for(int j=0;j<one.length;j++){
                    if(one[j].contains(geneName)){
                        bwone[j].write(temp);bwone[j].newLine();
                        continue;
                    }
                }
                for(int j=0;j<AB.length;j++){
                    if(AB[j].contains(geneName)){
                        bwAB[j].write(temp);bwAB[j].newLine();
                        continue;
                    }
                }
//                    for(int j=0;j<AD.length;j++){
//                        if(AD[j].contains(geneName)){
//                            for(int i=0;i<infoList.size();i++){
//                                bwAD[j].write(infoList.get(i).toString());bwAD[j].newLine();
//                            }
//                            infoList.clear();
//                            geneName=temp.split(" ")[0];
//                            infoList.add(temp);
//                            continue;
//                        }
//                    }
//                    for(int j=0;j<BD.length;j++){
//                        if(BD[j].contains(geneName)){
//                            for(int i=0;i<infoList.size();i++){
//                                bwBD[j].write(infoList.get(i).toString());bwBD[j].newLine();
//                            }
//                            infoList.clear();
//                            geneName=temp.split(" ")[0];
//                            infoList.add(temp);
//                            continue;
//                        }
//                    }
//                    for(int j=0;j<ABD.length;j++){
//                        if(ABD[j].contains(geneName)){
//                            for(String s : infoList){
//                                bwABD[j].write(s);bwABD[j].newLine();
//                            }
//                            infoList.clear();
//                            geneName=temp.split(" ")[0];
//                            infoList.add(temp);
//                            continue;
//                        }
//                    }
//                    infoList.clear();
//                    geneName=temp.split(" ")[0];
//                    infoList.add(temp);

            }
            br.close();
            for(int i=0;i<one.length;i++){
                bwone[i].flush();bwone[i].close();
            }
            for(int i=0;i<bwAB.length;i++){
                bwAB[i].flush();bwAB[i].close();
            }
//            for(int i=0;i<bwAD.length;i++){
//                bwAD[i].flush();bwAD[i].close();
//            }
//            for(int i=0;i<bwBD.length;i++){
//                bwBD[i].flush();bwBD[i].close();
//            }
//            for(int i=0;i<bwABD.length;i++){
//                bwABD[i].flush();bwABD[i].close();
//            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void parseHomo(){
        String ouptutDirS=null;
        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/homoGene/7_homoExpre.txt");
        HashSet[] one=new HashSet[3];
        HashSet[] AB=new HashSet[2];HashSet[] AD=new HashSet[2];HashSet[] BD=new HashSet[2];
        for(int i=0;i<AB.length;i++){
            AB[i]=new HashSet();AD[i]=new HashSet();BD[i]=new HashSet();
        }
        HashSet[] ABD=new HashSet[3];
        for(int i=0;i<one.length;i++){
            one[i]=new HashSet();
            ABD[i]=new HashSet();
        }
        for( int i=0;i<rt.getRowNumber();i++){
            if(rt.getCell(i, 3).equals("1")){
                if(rt.getCell(i, 4).equals("A")){
                    one[0].add(rt.getCell(i, 0));
                }else{
                    if(rt.getCell(i, 4).equals("B")){
                        one[1].add(rt.getCell(i, 1));
                    }else{
                        one[2].add(rt.getCell(i, 2));
                    }
                }
                continue;
            }
            if(rt.getCell(i, 3).equals("2")){
                if(rt.getCell(i, 4).equals("AB")){
                    AB[0].add(rt.getCell(i, 0));
                    AB[1].add(rt.getCell(i, 1));
                }else{
                    if(rt.getCell(i, 4).equals("AD")){
                        AD[0].add(rt.getCell(i, 0));
                        AD[1].add(rt.getCell(i, 2));
                    }else{
                        BD[0].add(rt.getCell(i, 1));
                        BD[1].add(rt.getCell(i, 2));
                    }
                }
                continue;
            }
            if(rt.getCell(i, 3).equals("3")){
                ABD[0].add(rt.getCell(i, 0));
                ABD[1].add(rt.getCell(i, 1));
                ABD[2].add(rt.getCell(i, 2));
            }
        }
        try{
            BufferedReader br =IOUtils.getTextGzipReader("/Users/xujun/Desktop/eQTL/N344/nominals.all.aFC.txt.gz");
            BufferedWriter[] bwone = new BufferedWriter[3];
            BufferedWriter[] bwAB = new BufferedWriter[2];
            bwAB[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/0.txt");
            bwAB[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AB/1.txt");
            BufferedWriter[] bwAD = new BufferedWriter[2];
            bwAD[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AD/0.txt");
            bwAD[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/AD/2.txt");
            BufferedWriter[] bwBD = new BufferedWriter[2];
            bwBD[0]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/BD/1.txt");
            bwBD[1]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/BD/2.txt");
            BufferedWriter[] bwABD =new BufferedWriter[3];
            for(int i=0;i<bwABD.length;i++){
                bwone[i]=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/one/"+i+".txt");
                bwABD[i] = IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/"+i+".txt");
            }
            String temp=null;String geneName=null;
            ArrayList<String> infoList = new ArrayList();
            temp=br.readLine();
            geneName=temp.split(" ")[0];
            while((temp=br.readLine())!=null){
                if(geneName.equals(temp.split(" ")[0])){
                    infoList.add(temp);continue;
                }else{
                    for(int j=0;j<one.length;j++){
                        if(one[j].contains(geneName)){
                            for(String s:infoList) {
                                bwone[j].write(s);bwone[j].newLine();
                            }
                            infoList.clear();
                            geneName=temp.split(" ")[0];
                            infoList.add(temp); //这个地方
                            continue;
                        }
                    }
                    for(int j=0;j<AB.length;j++){
                        if(AB[j].contains(geneName)){
                            for(int i=0;i<infoList.size();i++){
                                bwAB[j].write(infoList.get(i).toString());bwAB[j].newLine();
                            }
                            infoList.clear();
                            geneName=temp.split(" ")[0];
                            infoList.add(temp);
                            continue;
                        }
                    }
                    for(int j=0;j<AD.length;j++){
                        if(AD[j].contains(geneName)){
                            for(int i=0;i<infoList.size();i++){
                                bwAD[j].write(infoList.get(i).toString());bwAD[j].newLine();
                            }
                            infoList.clear();
                            geneName=temp.split(" ")[0];
                            infoList.add(temp);
                            continue;
                        }
                    }
                    for(int j=0;j<BD.length;j++){
                        if(BD[j].contains(geneName)){
                            for(int i=0;i<infoList.size();i++){
                                bwBD[j].write(infoList.get(i).toString());bwBD[j].newLine();
                            }
                            infoList.clear();
                            geneName=temp.split(" ")[0];
                            infoList.add(temp);
                            continue;
                        }
                    }
                    for(int j=0;j<ABD.length;j++){
                        if(ABD[j].contains(geneName)){
                            for(String s : infoList){
                                bwABD[j].write(s);bwABD[j].newLine();
                            }
                            infoList.clear();
                            geneName=temp.split(" ")[0];
                            infoList.add(temp);
                            continue;
                        }
                    }
                    infoList.clear();
                    geneName=temp.split(" ")[0];
                    infoList.add(temp);
                }
            }
            br.close();
            for(int i=0;i<one.length;i++){
                bwone[i].flush();bwone[i].close();
            }
            for(int i=0;i<bwAB.length;i++){
                bwAB[i].flush();bwAB[i].close();
            }
            for(int i=0;i<bwAD.length;i++){
                bwAD[i].flush();bwAD[i].close();
            }
            for(int i=0;i<bwBD.length;i++){
                bwBD[i].flush();bwBD[i].close();
            }
            for(int i=0;i<bwABD.length;i++){
                bwABD[i].flush();bwABD[i].close();
            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void homoExpre(){
        String inputFileS="/Users/xujun/Desktop/eQTL/N344/7_nor_countResult.txt";
        String outputFileS="/Users/xujun/Desktop/eQTL/N344/homeExpre.txt";
        HashSet gene = new HashSet();String temp=null;
        HashMap[] geneExpre = new HashMap[8];
        for(int i=0;i<geneExpre.length;i++){
            geneExpre[i]=new HashMap();
        }
        try{
            BufferedReader br =IOUtils.getTextReader(inputFileS);
            BufferedWriter bw =IOUtils.getTextWriter(outputFileS);
            br.readLine();int count=0;double total=0;
            while((temp=br.readLine())!=null){
                for(int i =1;i<temp.split("\t").length;i++){
                    if(Double.valueOf(temp.split("\t")[i])>0){
                        count++;
                        total+=Double.valueOf(temp.split("\t")[i]);
                    }
                }
                if(count>20){
                    gene.add(temp.split("\t")[0]);
                }
                DecimalFormat decFor = new DecimalFormat("0.000");
                if(temp.split("\t")[0].charAt(7)=='U'){
                    int chr=0;
                    geneExpre[chr].put(temp.split("\t")[0], decFor.format(total/count));
                }else{
                    int chr =Integer.valueOf(temp.split("\t")[0].charAt(7));
                    geneExpre[chr].put(temp.split("\t")[0], decFor.format(total/count));
                }
                count=0;total=0;
            }
            br.close();
            RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/homologyGene.txt");
            for(int i =0;i<rt.getRowNumber();i++){
                StringBuilder sb= new StringBuilder();
                for(int j =0;j<rt.getColumnNumber();j++){
                    if(gene.contains(rt.getCell(i, j))){
                        sb.append(rt.getCell(i, j)+"\t");
                    }else{
                        sb.append(0+"\t");
                    }
                }
                bw.write(sb.toString().replaceAll("\\s+$", ""));
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }

    }
    public void homoNominal(){
        String inputFileS="/Users/xujun/Desktop/eQTL/N344/homoGene/7_homoExpre.txt";
        String outputFileS="/Users/xujun/Desktop/eQTL/N344/homoGene/nominals.txt";
        HashSet gene = new HashSet();String temp=null;
        try{
            BufferedReader br =IOUtils.getTextReader(inputFileS);
            BufferedWriter bw =IOUtils.getTextWriter(outputFileS);
            br.readLine();int count=0;
            while((temp=br.readLine())!=null){
                if(temp.split("\t")[3].equals("0"))continue;
                if(temp.split("\t")[4].contains("A"))
                    for(int i =1;i<temp.split("\t").length;i++){
                        if(Double.valueOf(temp.split("\t")[i])>0){
                            count++;
                        }
                    }
                if(count>20){
                    gene.add(temp.split("\t")[0]);
                }
                count=0;
            }
            br.close();
            RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/homologyGene.txt");
            for(int i =0;i<rt.getRowNumber();i++){
                StringBuilder sb= new StringBuilder();
                for(int j =0;j<rt.getColumnNumber();j++){
                    if(gene.contains(rt.getCell(i, j))){
                        sb.append(rt.getCell(i, j)+"\t");
                    }else{
                        sb.append(0+"\t");
                    }
                }
                bw.write(sb.toString().replaceAll("\\s+$", ""));
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }

    }
    public void distri(){
//        String inputFileS="/data1/home/junxu/eQTL/FastQTL2/vst/5nominalsThred/result/nominals.txt.gz";
//        String outputFileS="/data1/home/junxu/eQTL/FastQTL2/vst/nominals.distributation.txt";
//        GeneFeature gf = new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
        String inputFileS="/Users/xujun/Desktop/eQTL/N344/nominals.all.aFC.txt.gz";
        String outputFileS="/Users/xujun/Desktop/eQTL/N344/nominals.distributation.txt";
        GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/eQTL/N344/wheat_v1.1_Lulab.gff3");
        String temp=null;
        try{
            BufferedReader br =IOUtils.getTextGzipReader(inputFileS);
            BufferedWriter bw =IOUtils.getTextWriter(outputFileS);
            int interC=0;int threeC=0;int fiveC=0;int CDSC=0;int introC=0;
            int count=0;
            while((temp=br.readLine())!=null){
//                 if(!temp.split(" ")[0].equals("TraesCS5B02G310500"))continue;
                int i=gf.getGeneIndex(temp.split(" ")[0]);
                if(!gf.isWithinThisGene(i, Integer.valueOf(temp.split(" ")[1].split("_")[1]), Integer.valueOf(temp.split(" ")[1].split("_")[2]))){
                    interC++;
                    continue;
                }
                int j=gf.getLongestTranscriptIndex(i);
                int pos=Integer.valueOf(temp.split(" ")[1].split("_")[2]);
                if(gf.get5UTRList(i, j).size()!=0){
                    List<Range> fU=gf.get5UTRList(i, j);
                    for(int a=0;a<fU.size();a++){
                        int start=fU.get(a).getRangeStart();
                        int end = fU.get(a).getRangeEnd();
                        if(pos>=start && pos<=end){
                            fiveC++;continue;
                        }
                    }
                }
                if(gf.getCDSList(i, j).size()!=0){
                    List<Range> CDS=gf.getCDSList(i, j);
                    for(int a=0;a<CDS.size();a++){
                        int start=CDS.get(a).getRangeStart();
                        int end = CDS.get(a).getRangeEnd();
                        if(pos>=start && pos<=end){
                            CDSC++;continue;
                        }
                    }
                }
                if(gf.get3UTRList(i, j).size()!=0){
                    List<Range> tU=gf.get3UTRList(i, j);
                    for(int a=0;a<tU.size();a++){
                        int start=tU.get(a).getRangeStart();
                        int end = tU.get(a).getRangeEnd();
                        if(pos>=start && pos<=end){
                            threeC++;continue;
                        }
                    }
                }
                introC++;
            }
            bw.write("intergenetic"+"\t"+interC);bw.newLine();
            bw.write("fiveUTR"+"\t"+fiveC);bw.newLine();
            bw.write("CDS"+"\t"+CDSC);bw.newLine();
            bw.write("Intro"+"\t"+introC);bw.newLine();
            bw.write("threeUTR"+"\t"+threeC);bw.newLine();
            bw.flush();bw.close();
            br.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
            System.out.println(temp);
        }
    }
    public void intergenicPattern(){
//        String inputFileS="/data1/home/junxu/eQTL/FastQTL2/vst/5nominalsThred/result/nominals.txt.gz";
//        String outputFileS="/data1/home/junxu/eQTL/FastQTL2/vst/nominals.distributation.txt";
//        GeneFeature gf = new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
        int subGenome=-1;String pattern=null;String mode=null;
        pattern="ABD";subGenome=0;
        mode=pattern+"/"+subGenome;
        String inputFileS="/Users/xujun/Desktop/eQTL/N344/homoGene/"+mode+".txt";
        String outputFileUp="/Users/xujun/Desktop/eQTL/N344/homoGene/"+mode+".intergenic.up.distributation.txt";
        String outputFileDown="/Users/xujun/Desktop/eQTL/N344/homoGene/"+mode+".intergenic.down.distributation.txt";
        String outputFileUpEf="/Users/xujun/Desktop/eQTL/N344/homoGene/"+mode+".intergenic.up.ef.txt";
        String outputFileDownEf="/Users/xujun/Desktop/eQTL/N344/homoGene/"+mode+".intergenic.down.ef.txt";
        GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/eQTL/N344/wheat_v1.1_Lulab.gff3");
//        int a=gf.getGeneStart(gf.getGeneIndex("TraesCS1A02G004400"));
//        int b=gf.getGeneEnd(gf.getGeneIndex("TraesCS1A02G004400"));
        int[] countUp = new int [1000]; double[] efCUp = new double[1000];int up=0;
        int[] countDown = new int [1000]; double[] efCDown = new double[1000];int down=0;
        String temp=null;int pos=0;int start=0; int end =0;
        String geneName=null;double ef=0;
        try{
            BufferedReader br =IOUtils.getTextReader(inputFileS);
            BufferedWriter bwUp =IOUtils.getTextWriter(outputFileUp);
            BufferedWriter bwDown =IOUtils.getTextWriter(outputFileDown);
            BufferedWriter bwUpEf =IOUtils.getTextWriter(outputFileUpEf);
            BufferedWriter bwDownEf =IOUtils.getTextWriter(outputFileDownEf);
            while((temp=br.readLine())!=null){
                geneName=temp.split(" ")[0];
                if(temp.split("\t")[1].equals("nan")){
                    ef=0;
                }else{
                    ef=Double.valueOf(temp.split("\t")[1]);
                }
                int i=gf.getGeneIndex(geneName);
                if(!gf.isWithinThisGene(i, Integer.valueOf(temp.split(" ")[1].split("_")[1]), Integer.valueOf(temp.split(" ")[1].split("_")[2]))){
                    pos=Integer.valueOf(temp.split(" ")[1].split("_")[2]);
                    if(gf.getGeneStrand(i)==1){//1表示的是正链
                        start =gf.getGeneStart(i);
                        if(start>=pos){
                            int chunk=(start-pos)/1000;
                            countUp[chunk]++;up++;
                            efCUp[chunk]+=ef;
                        }else{
                            end=gf.getGeneEnd(i);
                            int chunk=(pos-end)/1000;
                            countDown[chunk]++;down++;
                            efCDown[chunk]+=ef;
                        }
                    }else{
                        start=gf.getGeneEnd(i);
                        if(start<=pos){
                            int chunk=(pos-start)/1000;
                            countUp[chunk]++;up++;
                            efCUp[chunk]+=ef;
                        }else{
                            end=gf.getGeneStart(i);
                            int chunk=(end-pos)/1000;
                            countDown[chunk]++;down++;
                            efCDown[chunk]+=ef;
                        }
                    }

                }
            }
            br.close();
            DecimalFormat decFor = new DecimalFormat("0.000000");
            for(int i=0;i<countUp.length;i++){
                bwUp.write(countUp[i]+"\n");
                if(efCUp[i]==0){
                    bwUpEf.write(0+"\n");
                }else{
                    bwUpEf.write(decFor.format((efCUp[i]/countUp[i])*1000000/1000000)+"\n");
                }
            }
            bwUp.flush();bwUp.close();bwUpEf.flush();bwUpEf.close();
            System.out.println(up);
            for(int i=0;i<countDown.length;i++){
                bwDown.write(countDown[i]+"\n");
                if(efCDown[i]==0){
                    bwDownEf.write(0+"\n");
                }else{
                    bwDownEf.write(decFor.format((efCDown[i]/countDown[i])*1000000/1000000)+"\n");
                }
            }
            bwDown.flush();bwDown.close();bwDownEf.flush();bwDownEf.close();
            System.out.println(down);
        }
        catch(Exception ex){
            System.out.println(geneName);
            System.out.println(pos);
            System.out.println(start);
            ex.getStackTrace();
        }
    }
    public void parseNorminals(){
//         String inputFileS="/data1/home/junxu/eQTL/7FastQTL/5nominals/homoGene.nominals.all.txt.gz";
//         String outputFileS="/data1/home/junxu/eQTL/7FastQTL/5nominals/homoGene.nominals.gene.news.txt.gz";
//         String outputFileS1="/data1/home/junxu/eQTL/7FastQTL/5nominals/homoGene.nominals.SNP.news.txt.gz";
//           RowTable rt = new RowTable("/data1/home/junxu/eQTL/7FastQTL/5nominals/geneSet.txt");
        String inputFileS="/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/2.txt";
        String outputFileS="/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/2.nominals.gene.news.txt";
        String outputFileS1="/Users/xujun/Desktop/eQTL/N344/homoGene/ABD/2.nominals.SNP.news.txt";
        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/getSub.txt");
        HashMap chroSub=new HashMap();
        for(int i=0;i<rt.getRowNumber();i++){
            chroSub.put(rt.getCell(i, 0), rt.getCell(i,0)+"_"+rt.getCell(i,1));
        }
        try{
//             BufferedReader br =IOUtils.getTextGzipReader(inputFileS);
//             BufferedWriter bw =IOUtils.getTextGzipWriter(outputFileS);
//             BufferedWriter bw1 =IOUtils.getTextGzipWriter(outputFileS1);
            BufferedReader br =IOUtils.getTextReader(inputFileS);
            BufferedWriter bw =IOUtils.getTextWriter(outputFileS);
            BufferedWriter bw1 =IOUtils.getTextWriter(outputFileS1);
            HashSet gene = new HashSet(); HashMap hmGene = new HashMap();
            HashSet SNP = new HashSet(); HashMap hmSNP = new HashMap();
            String temp=null;
            int countGene=0;int countSNP=0;
            while((temp=br.readLine())!=null){
                if(!gene.contains(temp.split(" ")[0])){
                    gene.add(temp.split(" ")[0]);
                    countGene=1;
                }else{
                    countGene++;
                }
                hmGene.put(temp.split(" ")[0], countGene);
                if(!SNP.contains(temp.split(" ")[1])){
                    SNP.add(temp.split(" ")[1]);
                    countSNP=1;
                }else{
                    countSNP++;
                }
                hmSNP.put(temp.split(" ")[1], countSNP);
            }
            br.close();
            Iterator<String> itr = gene.iterator();
            while(itr.hasNext()){
                String a= itr.next();
                bw.write(a+"\t");
                bw.write(hmGene.get(a).toString());
                bw.newLine();
            }
            bw.flush();bw.close();
            Iterator<String> itr1 = SNP.iterator();
            while(itr1.hasNext()){
                String a= itr1.next();
                bw1.write(a.split("_")[1]+"\t"+a.split("_")[2]+"\t");
                bw1.write(hmSNP.get(a).toString()+"\t");
                bw1.write(chroSub.get(a.split("_")[1]).toString().split("_")[1]);
                bw1.newLine();
            }
            bw1.flush();bw1.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void parseHomoNorminals(){
//         String inputFileS="/data1/home/junxu/eQTL/7FastQTL/5nominals/nominals.all.txt.gz";
//         String inputFileS1="/data1/home/junxu/eQTL/7FastQTL/5nominals/honoGene.txt";
//         String outputFileS="/data1/home/junxu/eQTL/7FastQTL/5nominals/D/D_nominals.gene.news.txt.gz";
        String inputFileS="/Users/xujun/Desktop/eQTL/N344/nominals.all.txt.gz";
        String inputFileS1="/Users/xujun/Desktop/eQTL/N344/homoGene.txt";
        String outputFileS="/Users/xujun/Desktop/eQTL/N344/homoGene.nominals.all.txt.gz";
//        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/N344/getSub.txt");
        try{
            BufferedReader br =IOUtils.getTextGzipReader(inputFileS);
            BufferedReader br1 =IOUtils.getTextReader(inputFileS1);
            BufferedWriter bw =IOUtils.getTextGzipWriter(outputFileS);
            String temp=null;
            HashSet homoGene=new HashSet();
            while((temp=br1.readLine())!=null){
                for(int i=0;i<temp.split("\t").length;i++){
                    homoGene.add(temp.split("\t")[i]);
                }
            }
            br1.close();
            while((temp=br.readLine())!=null){
                if(homoGene.contains(temp.split(" ")[0])){
                    bw.write(temp);bw.newLine();
                }
            }
            br.close();
            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void subSetPermutations(){
        String inputDirS="/data1/home/junxu/eQTL/FastQTL2/genotype/";
        String outputDirS="/data1/home/junxu/eQTL/FastQTL2/vst/0.05subGenotype/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "DS.right.vcf.gz");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        RowTable rt = new RowTable("/data1/home/junxu/eQTL/FastQTL2/snpNum.txt");
        HashMap hm = new HashMap ();
        for(int i=0;i<rt.getRowNumber();i++){
            hm.put(i, rt.getCell(i, 0));
//            System.out.println(hm.get(i+1));
        }
        try{
            fList.stream().forEach(f -> {
                try{
                    BufferedWriter bw =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName()).getAbsolutePath());
                    String b=hm.get(Integer.valueOf(f.getName().split("\\.")[0])).toString();
                    int snpNumber=Integer.valueOf(b)-46;
                    System.out.println(snpNumber);
                    if(snpNumber<0){
                        snpNumber=0;
                    }
                    List<Integer> total =new ArrayList();
                    for(int i=4;i<snpNumber+5;i++){
                        total.add(i);
                    }
                    Collections.shuffle(total);
                    int sub = (int)Math.round(snpNumber*0.05);
                    System.out.println(f.getName()+"\t"+snpNumber+"\t"+sub);
                    List subList = total.subList(0, sub+1);
                    String temp=null; int row=0;
                    BufferedReader br =IOUtils.getTextGzipReader(f.getAbsolutePath());
                    while((temp=br.readLine())!=null){
                        char s = temp.charAt(0);
                        if ((int)s < 48 || (int)s > 57){
                            bw.write(temp);bw.newLine();
                            continue;
                        }
                        if(subList.contains(row)){
                            bw.write(temp);bw.newLine();
                        }
                        row++;
                    }
                    br.close();
                    bw.flush();bw.close();
                }
                catch(Exception ex){
                    ex.getStackTrace();
                }
            });
        }
        catch(Exception ex){
            ex.getStackTrace();
        }

    }
    public void geneLoca(){
        GeneFeature gf = new GeneFeature("/data1/home/xujun/wheat_v1.1_Lulab.gff3");
        String inputFile="/data1/home/junxu/eQTL/FastQTL2/vst/P7_0.2_vst.sorted.bed.gz";
        String outputDirS="/data1/home/junxu/eQTL/FastQTL2/vst/P7ParseBed/";
        try{
            BufferedReader br = IOUtils.getTextGzipReader(inputFile);
            String temp=null;
            String [] tem = null;List<String> tList=new ArrayList();
            try{
                temp=br.readLine();
                BufferedWriter[] bw = new BufferedWriter[45];
                for(int i=0;i<45;i++){
                    bw[i] = IOUtils.getTextWriter(new File(outputDirS,i+".sorted.bed").getAbsolutePath());
                    bw[i].write(temp);bw[i].newLine();
                }
                while((temp=br.readLine())!=null){
                    tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    bw[Integer.valueOf(tem[0])].write(temp);bw[Integer.valueOf(tem[0])].newLine();
                }
                br.close();
                for(int i=0;i<bw.length;i++){
                    bw[i].flush();bw[i].close();
                }
            }
            catch(Exception ex){
                ex.getStackTrace();
            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }

    }
    public String tail( File file, int lines) {
        java.io.RandomAccessFile fileHandler = null;
        try {
            fileHandler =
                    new java.io.RandomAccessFile( file, "r" );
            long fileLength = fileHandler.length() - 1;
            StringBuilder sb = new StringBuilder();
            int line = 0;
            for(long filePointer = fileLength; filePointer != -1; filePointer--){
                fileHandler.seek( filePointer );
                int readByte = fileHandler.readByte();
                if( readByte == 0xA ) {
                    if (filePointer < fileLength) {
                        line = line + 1;
                    }
                } else if( readByte == 0xD ) {
                    if (filePointer < fileLength-1) {
                        line = line + 1;
                    }
                }
                if (line >= lines) {
                    break;
                }
                sb.append( ( char ) readByte );
            }
            String lastLine = sb.reverse().toString();
            return lastLine;
        } catch( java.io.FileNotFoundException e ) {
            e.printStackTrace();
            return null;
        } catch( java.io.IOException e ) {
            e.printStackTrace();
            return null;
        }
        finally {
            if (fileHandler != null )
                try {
                    fileHandler.close();
                } catch (IOException e) {
                }
        }
    }
    public void parseBed (){
        String inputFile="/data1/home/junxu/eQTL/7FastQTL/7_nor_box_sorted.bed.gz";
        String outputDirS="/data1/home/junxu/eQTL/7FastQTL/7ParseBed/";
        try{
            BufferedReader br = IOUtils.getTextGzipReader(inputFile);
            String temp=null;
            String [] tem = null;List<String> tList=new ArrayList();
            try{
                temp=br.readLine();
                BufferedWriter[] bw = new BufferedWriter[45];
                for(int i=0;i<45;i++){
                    bw[i] = IOUtils.getTextWriter(new File(outputDirS,i+".sorted.bed").getAbsolutePath());
                    bw[i].write(temp);bw[i].newLine();
                }
                while((temp=br.readLine())!=null){
                    tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    bw[Integer.valueOf(tem[0])].write(temp);bw[Integer.valueOf(tem[0])].newLine();
                }
                br.close();
                for(int i=0;i<bw.length;i++){
                    bw[i].flush();bw[i].close();
                }
            }
            catch(Exception ex){
                ex.getStackTrace();
            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void changeVCF(){
        String inputDirS="/data1/home/junxu/eQTL/7FastQTL/5nominals/92DS/";
//        String inputDirS="/data1/home/junxu/eQTL/FastQTL2/FastQTL/0.snp.DS.vcf";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "DS.vcf");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        try{
            fList.stream().forEach(f -> {
                try{
                    BufferedWriter bw =IOUtils.getTextWriter(inputDirS+f.getName().replace("DS.vcf", "DS.right.vcf"));
                    BufferedReader br =IOUtils.getTextReader(f.getAbsolutePath());
                    String temp=null;
                    String [] tem = null;List<String> tList=new ArrayList();
                    bw.write("##fileformat=VCFv4.1");bw.newLine();
                    bw.write("##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage from MaCH/Thunder\">");bw.newLine();
                    while((temp=br.readLine())!=null){
                        tList= PStringUtils.fastSplit(temp);
                        tem = tList.toArray(new String[tList.size()]);
                        char s = temp.charAt(0);
                        if ((int)s < 48 || (int)s > 57){
                            StringBuilder header = new StringBuilder();
                            for(int i=0;i<tem.length;i++){
                                if(i==2){
                                    header.append("ID"+"\t");
                                }
                                if(i==4){
                                    header.append("QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t");
                                }
                                header.append(tem[i]+"\t");
                            }
                            bw.write(header.toString().replaceAll("\\s+$", ""));
                            bw.newLine();
                            continue;
                        }
                        StringBuilder out = new StringBuilder();
                        for(int i=0;i<tem.length;i++){
                            if(i==2){
                                out.append("snp_"+tem[0]+"_"+tem[1]+"\t");
                            }
                            if(i==4){
                                out.append("100"+"\t"+"PASS"+"\t"+"INFO"+"\t"+"DS"+"\t");
                            }
                            out.append(tem[i]+"\t");
                        }
                        bw.write(out.toString().replaceAll("\\s+$", ""));
                        bw.newLine();
                    }
                    br.close();
                    bw.flush();bw.close();
                }
                catch(Exception ex){
                    ex.getStackTrace();
                }

            });
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void nominalsMode(){
        String inputDirS="/data1/home/junxu/eQTL/FastQTL2/genotype/";
        String outputDirS="/data1/home/junxu/eQTL/FastQTL2/vst/10nominalsThred/command/";
//        String inputDirS="/data1/home/junxu/eQTL/FastQTL2/FastQTL/0.snp.DS.vcf";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".snp.maf001.mis01.DS.right.vcf.gz");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        try{
            fList.stream().forEach(f ->{
                BufferedWriter bw = IOUtils.getTextWriter(outputDirS+"/"+f.getName().split("\\.")[0]+".sh");
                try{
                    for(int j=1;j<11;j++){
                        StringBuilder sb = new StringBuilder();
                        sb.append("/data1/home/junxu/software/FastQTL/bin/fastQTL");
                        sb.append(" --vcf "+f.getAbsolutePath());
                        sb.append(" --bed /data1/home/junxu/eQTL/FastQTL2/vst/P7ParseBed/"+f.getName().split("\\.")[0]+".sorted.bed.gz");
                        sb.append(" --chunk "+j+" 10");//--permute 100 100000
                        sb.append(" --cov /data1/home/junxu/eQTL/FastQTL2/P7_10_covariance.txt.gz --threshold 0.00001");
                        sb.append(" --out "+"/data1/home/junxu/eQTL/FastQTL2/vst/10nominalsThred/result/");//.append(" --permute 100 100000")
                        sb.append(f.getName().split("\\.")[0]+"_"+j+"nominals.txt.gz ");
                        String command = sb.toString();
                        //                            System.out.println(command);
                        bw.write(command);bw.newLine();
                    }
                    bw.flush();bw.close();
                }
                catch(Exception ex){
                    ex.getStackTrace();
                }

            });
//                }
//            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void fastQTL(){
//        String inputDirS="/data1/home/junxu/eQTL/FastQTL2/vst/0.05/0.05subGenotype/";
        String outputDirS="/data1/home/junxu/eQTL/7FastQTL/5nominals/D/command/";
        String inputDirS="/data1/home/junxu/eQTL/genotype/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".snp.maf001.mis01.DS.right.vcf.gz");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        try{
//            File dir = new File(new File ("/data1/home/junxu/eWTL/FastQTL2/").getAbsolutePath());
//            int count=0;List<File> fL1 = new ArrayList(Arrays.asList());
//            for(int i=0;i<fList.size();i++){
//                if((i+1)/10==count && i!=fList.size()-1){//实现了同时跑10个HTSeq
//                    fL1.add(fList.get(i));
//                }else{
//                    fL1.add(fList.get(i));
            fList.stream().forEach(f ->{
//                        String chr=f.getName().split("\\.")[0];
//                        System.out.println(f.getName().split("\\.")[0]);
                BufferedWriter bw = IOUtils.getTextWriter(outputDirS+"/"+f.getName().split("\\.")[0]+".sh");
                try{
                    for(int j=1;j<11;j++){
                        StringBuilder sb = new StringBuilder();
                        sb.append("/data1/home/junxu/software/FastQTL/bin/fastQTL");
                        sb.append(" --vcf "+f.getAbsolutePath());
                        sb.append(" --bed /data1/home/junxu/eQTL/7FastQTL/7ParseBed/"+f.getName().split("\\.")[0]+".sorted.bed.gz");
                        sb.append(" --chunk "+j+" 10");//--permute 100 100000
                        sb.append(" --cov /data1/home/junxu/eQTL/7FastQTL/7_5_covariance.txt.gz");
                        sb.append(" --threshold 0.00001 --out "+"/data1/home/junxu/eQTL/7FastQTL/5nominals/D/result/");//.append(" --permute 100 100000")
                        sb.append(f.getName().split("\\.")[0]+"_"+j+"nominals.txt.gz --window 5000000 --exclude-samples /data1/home/junxu/eQTL/7FastQTL/file.exc ");
                        String command = sb.toString();
                        //                            System.out.println(command);
                        bw.write(command);bw.newLine();
                    }
                    bw.flush();bw.close();
                }
                catch(Exception ex){
                    ex.getStackTrace();
                }

            });
//                }
//            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void geneRegion(){
        String inputDirS="/data1/home/junxu/eQTL/7check/";
//        String inputDirS="/data1/home/junxu/eQTL/FastQTL2/genotype/41.snp.maf001.mis01.DS.vcf.gz";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "recode.vcf.gz");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        GeneFeature gf = new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
//        GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/eQTL/N344/wheat_v1.1_Lulab.gff3");
        HashIntIntMap[] posGeneMaps = new HashIntIntMap[45];
        for(int i=0;i<45;i++){
            posGeneMaps[i] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        }
        for(int i=0;i<gf.genes.length;i++){
            int index=gf.getGeneIndex(gf.genes[i].geneName);
            int chr=gf.genes[i].geneRange.chr;
            int start= gf.genes[i].geneRange.start;
            int end = gf.genes[i].geneRange.end;
            for(int j=start;j<end;j++){
                posGeneMaps[chr].put(j,index);
            }
        }
        int index =posGeneMaps[10].get(54303);
        String gene =gf.getGeneName(index);
        try{
            fList.stream().forEach(f -> {
                try{
                    String temp=null;String [] tem = null;List<String> tList=new ArrayList();
                    BufferedReader br =IOUtils.getTextGzipReader(f.getAbsolutePath());
                    BufferedWriter bw =IOUtils.getTextGzipWriter("/data1/home/junxu/eQTL/7check/"+f.getName().replace("recode.vcf.gz","geneRegion.vcf.gz"));
                    while((temp=br.readLine())!=null){
                        if(temp.startsWith("##") || temp.startsWith("#")){
                            bw.write(temp);bw.newLine();
                            continue;
                        }
                        tList= PStringUtils.fastSplit(temp);
                        tem = tList.toArray(new String[tList.size()]);
                        int pos=Integer.valueOf(Integer.valueOf(tem[1]));
                        if(posGeneMaps[Integer.valueOf(tem[0])].get(pos)!=-1){
                            bw.write(temp);bw.newLine();
                        }
                    }
                    br.close();
                    bw.flush();bw.close();
                }
                catch(Exception ex){
                    ex.getStackTrace();
                }

            });
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void sort(){
//        String inputDirS="/data1/home/junxu/eQTL/FastQTL2/genotype/";
        String inputDirS="/data1/home/junxu/eQTL/FastQTL2/FastQTL/0.snp.DS.vcf";
//        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "DS.vcf.gz");
//        List<File> fList = Arrays.asList(fs);
//        Collections.sort(fList);
        RowTable rt = new RowTable("/data1/home/junxu/eQTL/FastQTL2/snpNum.txt");
        HashMap hm = new HashMap ();
        for(int i=0;i<rt.getRowNumber();i++){
            hm.put(i+1, rt.getCell(i, 0));
        }
        RowTable rt1 = new RowTable("/data1/home/junxu/eQTL/FastQTL2/nameList.txt");
        HashMap hp = new HashMap ();
        for(int i=0;i<rt1.getRowNumber();i++){
            hp.put(rt1.getCell(i, 0), i);
        }
        try{
//            fList.stream().forEach(f -> {
            try{
//                    BufferedWriter bw =IOUtils.getTextGzipWriter("/data1/home/junxu/eQTL/FastQTL2/FastQTL/"+f.getName());
//                    String b=hm.get(Integer.valueOf(f.getName().split("\\.")[0])).toString();
                BufferedWriter bw =IOUtils.getTextGzipWriter("/data1/home/junxu/eQTL/FastQTL2/FastQTL/0.snp.maf001.mis01.DS.vcf.gz");
//                    String b=hm.get(Integer.valueOf(inputDirS.split("/")[inputDirS.split("/").length-1].split("\\.")[0])).toString();
//                    int snpNumber=Integer.valueOf(b)-45;
                int snpNumber=1739540;
                System.out.println(snpNumber);
                if(snpNumber<0){
                    snpNumber=0;
                }
                String temp=null; int row=0;
                double [][] count=new double [snpNumber-1][96];
                String [][] news= new String [snpNumber-1][4];
                String [] header=new String[100];
//                    BufferedReader br =IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedReader br =IOUtils.getTextReader(inputDirS);
                HashMap olderNew = new HashMap();
                String [] tem = null;List<String> tList=new ArrayList();
                while((temp=br.readLine())!=null){
                    tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    char s = temp.charAt(0);
                    if ((int)s < 48 || (int)s > 57){
                        for(int i=0;i<tem.length;i++){
                            if(i<4){
                                header[i]=tem[i];
                            }else{
                                int col=Integer.valueOf(hp.get(tem[i]).toString());
                                header[col+4]=tem[i];
                                olderNew.put(i, col);
                            }
                        }
                        continue;
                    }
                    for(int i=0;i<tem.length;i++){
                        if(i<4){
                            news[row][i]=tem[i];
                        }else{
                            count[row][Integer.valueOf(olderNew.get(i).toString())]=Double.valueOf(tem[i]);
                        }

                    }
                    row++;
                }
                br.close();
                StringBuilder head= new StringBuilder();
                for(int i=0;i<header.length;i++){
                    head.append(header[i]+"\t");
                }
                bw.write(head.toString().replaceAll("\\s+$", ""));
                bw.newLine();
                for(int i=0;i<count.length;i++){
                    StringBuilder sb= new StringBuilder();
                    for(int a=0;a<news[i].length;a++){
                        sb.append(news[i][a]+"\t");
                    }
                    for(int j=0;j<count[i].length;j++){
                        sb.append(count[i][j]+"\t");
                    }
                    bw.write(sb.toString().replaceAll("\\s+$", ""));
                    bw.newLine();
                }
                bw.flush();bw.close();
            }
            catch(Exception ex){
                ex.getStackTrace();
            }

//            });
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void filterSample(){
        String sampleInfor="/data1/home/junxu/eQTL/7FastQTL/7_sampleName.txt";
        String inputFiles="/data2/junxu/SNP/";
        String outputFiles="/data1/home/junxu/eQTL/7FastQTL/5nominals/92genotype";
        RowTable rt = new RowTable(sampleInfor);
        HashSet nameSet=new HashSet(); HashMap hp =new HashMap();
        for (int i=0;i<rt.getRowNumber();i++){
            nameSet.add(rt.getCell(i, 0));
            hp.put(rt.getCell(i, 0), rt.getCell(i, 1));
        }
        File[] fs = new File(inputFiles).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        fList.stream().forEach(f -> {
            try{
//             BufferedReader br=IOUtils.getTextReader(f.getAbsolutePath());
//             BufferedWriter bw = IOUtils.getTextWriter(new File(outputFiles,f.getName().replace("tempory.vcf","92.vcf")).getAbsolutePath());
                BufferedReader br=IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(new File(outputFiles,f.getName()).getAbsolutePath());
                String temp=null; String [] tem = null;List<String> tList=new ArrayList();
                StringBuilder header = new StringBuilder();
                HashSet colSet=new HashSet();
                while((temp=br.readLine())!=null){
                    if(temp.startsWith("##")){
                        bw.write(temp);bw.newLine();
                        continue;
                    }
                    if(temp.startsWith("#")){
                        StringBuilder out =new StringBuilder();
                        tList= PStringUtils.fastSplit(temp);
                        tem = tList.toArray(new String[tList.size()]);
                        for(int i=0;i<tem.length;i++){
                            if(i<9){
                                out.append(tem[i]+"\t");
                            }else{
                                if(tem[i].contains("AT")){
                                    String name= tem[i].substring(0,tem[i].length()-1);
                                    if(nameSet.contains(name)){
                                        out.append(hp.get(name)+"\t");
                                        colSet.add(i);
                                    }
                                }else{
                                    if(nameSet.contains(tem[i])){
                                        out.append(hp.get(tem[i])+"\t");
                                        colSet.add(i);
                                    }
                                }
                            }
                        }
                        bw.write(out.toString().replaceAll("\\s+$", ""));
                        bw.newLine();
                        continue;
                    }
                    StringBuilder out = new StringBuilder();
                    tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    for(int i=0;i<tem.length;i++){
                        if(i<9){
                            out.append(tem[i]+"\t");
                        }else{
                            if(colSet.contains(i)){
                                out.append(tem[i]+"\t");
                            }
                        }
                    }
                    bw.write(out.toString().replaceAll("\\s+$", ""));
                    bw.newLine();
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        });
    }
    public void subSetVCF(){
        String inputDirS="/data1/home/junxu/eQTL/7FastQTL/5nominals/";
//        String inputDirS="/Users/xujun/Desktop/eQTL/N344/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "homoGene.nominals");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
//        RowTable rt = new RowTable("/data2/junxu/geneSNP/SNPNum");
//        HashMap hm = new HashMap ();
//        for(int i=0;i<rt.getRowNumber();i++){
//            hm.put(i, rt.getCell(i, 0));
////            System.out.println(hm.get(i+1));
//        }
        try{
            BufferedWriter bw =IOUtils.getTextGzipWriter("/data1/home/junxu/eQTL/7FastQTL/5nominals/001homoGene.nominals.txt.gz");
//            BufferedWriter bw =IOUtils.getTextGzipWriter("/Users/xujun/Desktop/eQTL/N344/001.honoGene.nominals.txt.gz");
            fList.stream().forEach(f -> {
                try{
//                    System.out.println(f.getName());
//                    String b=hm.get(Integer.valueOf(f.getName().split("\\.")[0])).toString();
//                    int snpNumber=Integer.valueOf(b)-46;
//                    System.out.println(snpNumber);
//                    if(snpNumber<0){
//                        snpNumber=0;
//                    }
                    int snpNumber=7569599;
                    List<Integer> total =new ArrayList();
                    for(int i=0;i<snpNumber+1;i++){
                        total.add(i);
                    }
                    Collections.shuffle(total);
                    int sub = (int)Math.round(snpNumber*0.01);
                    System.out.println(f.getName()+"\t"+snpNumber+"\t"+sub);
                    List subList = total.subList(0, sub+1);
                    String temp=null; int row=0;
                    BufferedReader br =IOUtils.getTextGzipReader(f.getAbsolutePath());
                    while((temp=br.readLine())!=null){
//                        char s = temp.charAt(0);
//                        if ((int)s < 48 || (int)s > 57) continue;
                        if(subList.contains(row)){
                            bw.write(temp);bw.newLine();
                        }
                        row++;
                    }
                    br.close();
                }
                catch(Exception ex){
                    ex.getStackTrace();
                }

            });
            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.getStackTrace();
        }

    }
    public void change(){
        String sampleInfor="/Users/xujun/Desktop/eQTL/N343/P7_ZX-B18.txt";
        String inputFiles="/Users/xujun/Desktop/eQTL/workflow/P7_nor_countResult.txt";
        String outputFiles="/Users/xujun/Desktop/eQTL/N343/P7_nor_countResult.txt";
        RowTable rt = new RowTable(sampleInfor);
        HashMap ht = new HashMap();
        for(int i=0;i<rt.getRowNumber();i++){
            ht.put(rt.getCell(i, 1), rt.getCell(i, 0));
        }
        try{
            BufferedReader br = IOUtils.getTextReader(inputFiles);
            BufferedWriter bw = IOUtils.getTextWriter(outputFiles);
            String temp=null; String [] tem = null;
            StringBuilder header = new StringBuilder();
            List<String> tList=new ArrayList();HashSet nameList=new HashSet();
            temp=br.readLine();tList= PStringUtils.fastSplit(temp);
            tem = tList.toArray(new String[tList.size()]);
            header.append("Gene"+"\t");
            for(int i=1;i<tem.length;i++){
                header.append(ht.get(tem[i])+"\t");
            }
            bw.write(header.toString().replaceAll("\\s+$", ""));
            bw.newLine();
            while((temp=br.readLine())!=null){
                bw.write(temp);bw.newLine();
            }
            br.close();
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    public void filterNominal(){//根据GTEx的首席的过滤标准
        BufferedReader br1 =IOUtils.getTextGzipReader("/data1/home/junxu/eQTL/FastQTL/fastqtl.df.txt.gz");
        BufferedReader br2 =IOUtils.getTextReader("/data1/home/junxu/eQTL/FastQTL/15nominals.all.txt");
        BufferedWriter bw =IOUtils.getTextWriter(new File("/data1/home/junxu/eQTL/FastQTL/15nominals.all.filter.txt").getAbsolutePath());
        GeneFeature gf= new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
//        BufferedReader br =IOUtils.getTextReader("/Users/xujun/Desktop/15nominals.all.storey.sorted.txt");
//        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/15nominals.all.storey.sorted.pos.txt").getAbsolutePath());
//        GeneFeature gf= new GeneFeature("/Users/xujun/Desktop/eQTL/workflow/wheat_v1.1_Lulab.gff3");
        HashMap hs = new HashMap();
        String temp=null;
        try{
            br1.readLine();
            while((temp=br1.readLine())!=null){
                hs.put(temp.split("\t")[0], temp.split("\t")[12]);
            }
            while((temp=br2.readLine())!=null){
                if(Double.valueOf(temp.split(" ")[3].toString()) < Double.valueOf(hs.get(temp.split(" ")[0]).toString())){
                    bw.write(temp);bw.newLine();
                }
            }
            br1.close();br2.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void SNPLocus(){
//        BufferedReader br =IOUtils.getTextReader("/data1/home/junxu/eQTL/FastQTL/15nominals.all.storey.sorted.txt");
//        BufferedWriter bw =IOUtils.getTextWriter(new File("/data1/home/junxu/eQTL/FastQTL/15nominals.all.storey.sorted.pos.txt").getAbsolutePath());
//        GeneFeature gf= new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
        BufferedReader br =IOUtils.getTextReader("/Users/xujun/Desktop/15nominals.all.storey.sorted.txt");
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/15nominals.all.storey.sorted.pos.txt").getAbsolutePath());
        GeneFeature gf= new GeneFeature("/Users/xujun/Desktop/eQTL/workflow/wheat_v1.1_Lulab.gff3");
        String temp=null; int a=0;int b=0;int c=0;
        try{
            br.readLine();
            while((temp=br.readLine())!=null){
                int pos=Integer.valueOf(temp.split(" ")[2]);
                int index = gf.getGeneIndex(temp.split(" ")[0]);
                int start =gf.getGeneStart(index);
                int end = gf.getGeneEnd(index);
                if(pos>=start && pos< end){
                    bw.write(temp+"\t"+0);//0表示在基因内部 -1上游 1下游
                    bw.newLine();
                    b++;
                }else if(pos< start){
                    bw.write(temp+"\t"+ -1);bw.newLine();
                    a++;
                }else{
                    bw.write(temp+"\t"+ 1);bw.newLine();
                    c++;
                }
            }
            br.close();
            bw.flush();bw.close();
            System.out.println(a+"\t"+b+"\t"+c);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void merge1(){
        String inputfsK="/data2/yafei/Bam_Pileup/Pileup.out/";
        File[] fsK = new File(inputfsK).listFiles();
        List<File> fListK = new ArrayList(Arrays.asList());
        fsK = IOUtils.listFilesEndsWith(fsK, ".pileup");
        fListK=Arrays.asList(fsK);
        HashIntIntMap[] posNum=new HashIntIntMap[42];
        for(int i=0;i<posNum.length;i++){
            posNum[i] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        }
        BufferedWriter bw =IOUtils.getTextWriter(new File("/data2/yafei/Bam_Pileup/Pileup.out/merge1.txt").getAbsolutePath());
        fListK.stream().forEach(f -> {
            BufferedReader br =IOUtils.getTextReader(f.getAbsolutePath());
            String temp =null;
            try{
                while((temp=br.readLine())!=null){
                    int chr=Integer.valueOf(temp.split(" ")[0])-1;
                    int pos = Integer.valueOf(temp.split(" ")[1]);
                    if(posNum[chr].get(pos)!=-1){
                        int num=posNum[chr].get(pos);
                        posNum[chr].put(pos, num+Integer.valueOf(temp.split(" ")[2]));
                    }else{
                        int num=Integer.valueOf(temp.split(" ")[2]);
                        posNum[chr].put(pos,num);
                    }

                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try{
            for(int i=0;i<posNum.length;i++){
                for(int a: posNum[i].keySet()){
                    bw.write(i+1+"\t"+a+"\t"+posNum[i].get(a));
                    bw.newLine();
                }
            }

            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void merge(){
//        String inputfsK="/data2/yafei/Bam_Pileup/Pileup.out/";
//        String inputfsK="/Users/xujun/Desktop";
//        String inputfsK="/data1/home/junxu/";
        String inputfsK="/Users/xujun/Desktop/";
        File[] fsK = new File(inputfsK).listFiles();
        List<File> fListK = new ArrayList(Arrays.asList());
        fsK = IOUtils.listFilesEndsWith(fsK, ".pileup");
        fListK=Arrays.asList(fsK);
        HashDoubleIntMap posNum = HashDoubleIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        HashDoubleIntMap posPos = HashDoubleIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
//        BufferedWriter bw =IOUtils.getTextWriter(new File("/data2/yafei/Bam_Pileup/Pileup.out/merge.txt").getAbsolutePath());
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/merge.txt").getAbsolutePath());
        fListK.stream().forEach(f -> {
            BufferedReader br =IOUtils.getTextReader(f.getAbsolutePath());
            String temp =null; double pos=0;
            try{
                while((temp=br.readLine())!=null){
                    pos= Double.valueOf(temp.split(" ")[0]+"."+temp.split(" ")[1]);
                    if(posNum.get(pos)!=-1){
                        int num=posNum.get(pos);
                        posNum.put(pos, num+Integer.valueOf(temp.split(" ")[2]));
                    }else{
                        int num=Integer.valueOf(temp.split(" ")[2]);
                        posNum.put(pos,num);
                        posPos.put(pos,temp.split(" ")[1].length() );
                    }

                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try{
            for (double a: posNum.keySet()){
                bw.write(String.valueOf(a).split("\\.")[0]+"\t");
                bw.write(this.addZeroForNum(String.valueOf(a).split("\\.")[1],posPos.get(a))+"\t"+posNum.get(a));
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch (Exception e) {

            e.printStackTrace();
        }
    }

    public static String addZeroForNum(String str, int strLength) {
        int strLen = str.length();
        if (strLen < strLength) {
            while (strLen < strLength) {
                StringBuffer sb = new StringBuffer();
//            sb.append("0").append(str);// 左补0
                sb.append(str).append("0");//右补0
                str = sb.toString();
                strLen = str.length();
            }
        }

        return str;
    }
    public void NPos(){
        BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/tree/Axiom_wheat660.na34.annot.35N35.fa").getAbsolutePath());
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/tree/pos.txt").getAbsolutePath());
        try{
            String tem =null;String temp=null;
            while((temp=br.readLine())!=null){
                temp=br.readLine();
                for(int i=0;i<temp.length();i++){
                    if(temp.charAt(i)=='N'){
                        System.out.println(i);
                        bw.write(i);
                        bw.newLine();
                    }
                }

            }
            br.close();bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void changPos(){
        BufferedReader br =IOUtils.getTextReader(new File("/data1/home/junxu/tree/chr1D.filt.vcf").getAbsolutePath());
        BufferedWriter bw =IOUtils.getTextWriter(new File("/data1/home/junxu/tree/changechr1D.vcf").getAbsolutePath());
//         BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/tree/AA358.unique.mapped.ok.hmp.txt").getAbsolutePath());
//         BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/tree/changePosAA358.unique.mapped.ok.hmp.txt").getAbsolutePath());
        try{
            String [] tem =null;String temp=null;
            while((temp=br.readLine())!=null){
                if(!temp.startsWith("#")){
                    tem=temp.split("\t");
                    for(int i=0;i<tem.length;i++){
                        if(i!=0){
                            bw.write(tem[i]+"\t");
                        }else{
                            bw.write("15"+"\t");
                        }
                    }
                    bw.newLine();
                }else{
                    bw.write(temp);bw.newLine();
                }
            }
            br.close();bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }
    public void overlapSNP(){
        HashIntIntMap posline = null;
        posline = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        ArrayList infor = new ArrayList();List<Integer> h1= new ArrayList();
        BufferedReader br1 =IOUtils.getTextReader(new File("/data1/home/junxu/tree/changechr1B.vcf").getAbsolutePath());
        BufferedReader br2 =IOUtils.getTextReader(new File("/data1/home/junxu/tree/zhengchr1B.vcf").getAbsolutePath());
        BufferedWriter bw1 =IOUtils.getTextWriter(new File("/data1/home/junxu/tree/overjiaochr1B.vcf").getAbsolutePath());
        BufferedWriter bw2 =IOUtils.getTextWriter(new File("/data1/home/junxu/tree/overzhengchr1B.vcf").getAbsolutePath());
//         BufferedReader br1 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/tree/zhengchr1A.vcf").getAbsolutePath());
//         BufferedReader br2 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/tree/overzhengchr1A.vcf").getAbsolutePath());
//         BufferedWriter bw1 =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/tree/overjiaochr1A.vcf").getAbsolutePath());
//         BufferedWriter bw2 =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/tree/overzhengchr1A.vcf").getAbsolutePath());

        try{
//            for (int i=0;i<rt.getRowNumber();i++){
//                if(rt.getCell(i,0).equals("1B")){
//                    pos.add(rt.getCell(i, 1));
//                    posline.put(Integer.valueOf(rt.getCell(i, 1).toString()).intValue(),i);
//                }
//            }
            String [] temt =null;String tempt=null;int i=0;
            while((tempt=br1.readLine())!=null){
                temt=tempt.split("\t");
                if(temt[0].equals("8")){
                    i++;
                    infor.add(tempt);
                    posline.put(i,Integer.valueOf(temt[1].toString()).intValue());
                }
            }
            String [] tem =null;String temp=null;
            while((temp=br2.readLine())!=null){
                tem=temp.split("\t");
                if(!tem[0].startsWith("#")){
                    if(posline.containsValue(Integer.valueOf(tem[1]).intValue())){
                        h1=getKeyList(posline,Integer.valueOf(tem[1]).intValue());
                        for(int j=0;j<h1.size();j++){
                            bw1.write(infor.get(h1.get(j).intValue()-1).toString());
                            bw1.newLine();
                        }
                        bw2.write(temp);bw2.newLine();
                    }

                }

            }
            bw1.flush();bw1.close();br2.close();bw2.flush();bw2.close();
//            for(int i=1;i<geneNumber.length;i=i+2){
//                for(int j=0;j<geneNumber[i].length;j++){
//                    System.out.println(i+"\t"+geneNumber[i][j]);
//                }
//            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }
    public static List<Integer> getKeyList(HashIntIntMap map,int value){//the method to get key from value
        List<Integer> keyList = new ArrayList();
        for(int getKey: map.keySet()){
            if(map.get(getKey)==value){
                keyList.add(getKey);
            }
        }
        return keyList;
    }
    public void geneNumber(){
        RowTable rt = new RowTable("/Users/xujun/Desktop/tempory/readme.txt");
        HashMap subChrLength=new HashMap();
        for( int i=0;i< rt.getRowNumber();i++){
            subChrLength.put(rt.getCellAsInteger(i, 0), rt.getCellAsInteger(i, 2));
        }
        BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/IGVmaterial/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3").getAbsolutePath());
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/geneDis.txt").getAbsolutePath());
        int [][] geneNumber = new int [45][10];
//         GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/IGVmaterial/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        int erro=0;
        try{
            String [] tem =null;String temp=null;
            while((temp=br.readLine())!=null){
                tem=temp.split("\t");
                if(tem[2].startsWith("gene")){
                    if(Integer.valueOf(tem[0])%2!=0 && Integer.valueOf(tem[0])!=0){
                        int chrLength=Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])).toString())+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])+1).toString());
                        double pos = (double)Integer.valueOf(tem[3])/chrLength;
                        bw.write(pos+"");bw.newLine();
                    }else{
                        if(Integer.valueOf(tem[0])!=0){
                            int chrLength=Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])).toString())+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])-1).toString());
                            double pos =(double) (Integer.valueOf(tem[3])+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])-1).toString()))/chrLength;
//                            geneNumber[Integer.valueOf(tem[0])-1][pos]++;
                            bw.write(pos+"");bw.newLine();
                        }
                    }
                }

            }
            bw.flush();bw.close();br.close();
//            for(int i=1;i<geneNumber.length;i=i+2){
//                for(int j=0;j<geneNumber[i].length;j++){
//                    System.out.println(i+"\t"+geneNumber[i][j]);
//                }
//            }
        }
        catch (Exception e) {
            System.out.println(erro);
            e.printStackTrace();
        }
    }

    public void isSame(){
        BufferedReader br1 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/homoAB.txt").getAbsolutePath());
        BufferedReader br2 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/homoAD.txt").getAbsolutePath());
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/samA.txt").getAbsolutePath());
        List<String> nameList1 = new ArrayList();
        List<String> nameList2 = new ArrayList();
        try{
            String [] tem =null;String temp=null;
            while((temp=br1.readLine())!=null){
                if(temp.split(" ").length>1){
                    tem=temp.split(" ");
                    nameList1.add(tem[1].split("\t")[1]);
                }
            }
            while((temp=br2.readLine())!=null){
                if(temp.split(" ").length>1){
                    tem=temp.split(" ");
                    nameList2.add(tem[1].split("\t")[1]);
                }
            }
            for(int i=0;i<nameList2.size();i++){
                if(nameList1.contains(nameList2.get(i))){
                    bw.write(nameList2.get(i));bw.newLine();
                }
            }
            br1.close();br2.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void highDiverganceGene(){
        RowTable rt = new RowTable("/Users/xujun/Desktop/tempory/readme.txt");
        HashMap chrStart=new HashMap();
        HashMap chrEnd = new HashMap();
        GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/IGVmaterial/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        int geneNumber =gf.getGeneNumber();
        for( int i=0;i< rt.getRowNumber();i++){
            chrStart.put(rt.getCellAsInteger(i, 0), rt.getCellAsInteger(i, 4));
            chrEnd.put(rt.getCellAsInteger(i, 0), rt.getCellAsInteger(i, 5));
        }
        String temp=null;String position []= null;
        BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/homoAB.txt").getAbsolutePath());
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/AB.txt").getAbsolutePath());
        try{
            String [] tem =null;
            while((temp=br.readLine())!=null){
                if(temp.split(" ").length>1){
                    tem=temp.split(" ");
                    if(Integer.valueOf(tem[2])%2==0){
                        int pos = Integer.valueOf(tem[3])+Integer.valueOf(chrStart.get(Integer.valueOf(tem[2])).toString());
                        double pro=(double)pos/Integer.valueOf(chrEnd.get(Integer.valueOf(tem[2])).toString());
//                        System.out.println(pro);
                        bw.write(pro+"\t");
                    }else{
                        int pos=Integer.valueOf(tem[3]);
                        double pro=(double)pos/Integer.valueOf(chrEnd.get(Integer.valueOf(tem[2])+1).toString());
//                        System.out.println(pro);
                        bw.write(pro+"\t");
                    }
                    if(Integer.valueOf(tem[6])%2==0){
                        int pos = Integer.valueOf(tem[7])+Integer.valueOf(chrStart.get(Integer.valueOf(tem[6])).toString());
                        double pro=(double)pos/Integer.valueOf(chrEnd.get(Integer.valueOf(tem[6])).toString());
                        bw.write(pro+"");bw.newLine();
                    }else{
                        int pos=Integer.valueOf(tem[7]);
//                        int end = Integer.valueOf(chrEnd.get(Integer.valueOf(tem[6])+1).toString());
                        double pro=(double)pos/Integer.valueOf(chrEnd.get(Integer.valueOf(tem[6])+1).toString());
                        bw.write(pro+" ");bw.newLine();
                    }
//                    bw.write(tem[2]+"\t"+tem[6]);bw.newLine();
                }
            }
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void getList(){//to get the data from the distance matrix
        ColumnTable ct = new ColumnTable("/Users/xujun/Desktop/Tempory/g.txt");
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/geneticList.txt").getAbsolutePath());
        try{
            for(int i =0 ; i < ct.getColumnNumber() ; i++){
                for(int j =i ; j < ct.getRowNumber();j++){
                    bw.write(ct.getCellAsString(i, j));bw.newLine();
                }
            }
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void shScoreParse(){
        String inputFileS = "/Users/xujun/summarys.txt";
        String outputFileS = "/Users/xujun/shScore.txt";
        String temp=null;String shtemp=null;
        double sh=0;double shM=0;int k=0;int tempK=0;String num=null;
        try {
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            while((temp=br.readLine())!=null){
                if(temp.startsWith("logical")){
                    bw.write(String.valueOf(k));bw.newLine();
                    sh=0;shM=0;k=0;tempK=0;
                }
                if(temp.startsWith("Silhouette")){
                    tempK=Integer.valueOf(temp.split(" ")[5]);
                    br.readLine();br.readLine();br.readLine();
                    br.readLine();br.readLine();
                    shtemp=br.readLine();
                    if(!shtemp.contains("Individual")){
                        if(shtemp.contains("  ")){
                            sh=Double.valueOf(shtemp.split("  ")[3]);
                        }else{
                            sh=Double.valueOf(shtemp.split(" ")[3]);
                        }
                    }else{
                        br.readLine();shtemp=br.readLine();
                        if(shtemp.contains("  ")){
                            sh=Double.valueOf(shtemp.split("  ")[3]);
                        }else{
                            sh=Double.valueOf(shtemp.split(" ")[3]);
                        }
                    }

                    if(sh>shM){
                        shM=sh;k=tempK;
                        if(k==6){

                        }
                    }
                }
            }
            br.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            System.out.println(sh);
            e.printStackTrace();
        }
    }
    public void shScore (){
        String inputfsK="/Users/xujun/Desktop/eQTL/total/Homology/shScore";
        File[] fsK = new File(inputfsK).listFiles();
        List<File> fListK = new ArrayList(Arrays.asList());
        fsK = IOUtils.listFilesEndsWith(fsK, ".txt");
        for(int i=0;i<fsK.length;i++){
            if(!fsK[i].getName().contains("temp")){
                fListK.add(fsK[i]);
            }
        }
        RCaller callerL = new RCaller();
        callerL.setRscriptExecutable("/usr/local/bin/R");
        RCode rCodeL = new RCode();
        rCodeL.addRCode("library(cluster);library(factoextra);library(vegan);install.packages(\"Runiversal\");");
        callerL.setRCode(rCodeL);
        callerL.runOnly();
        List KPool =new ArrayList();
        fListK.stream().forEach(f -> {
            double sh=0;double shM=0;int k=0;
            for(int i=1;i<11;i++){
                RCaller caller = new RCaller();
                caller.setRscriptExecutable("//usr//local//bin//R");
                RCode rCode = new RCode();
                rCode.addRCode("a<-read.csv(\""+f+"\",sep=\"\\t\",header = F)");
                rCode.addRCode("results=kmeans(a,");
                rCode.addRCode(i+");dis = dist(a^2);");
                rCode.addRCode("sil = silhouette (results$cluster, dis);sum=summary(sil);sh<-sum$avg.width");
                caller.setRCode(rCode);
                caller.runAndReturnResult("sh");
                sh=caller.getParser().getAsDoubleArray("sh")[0];
                if(sh>shM){
                    shM=sh;
                }else{
                    k=i;
                    break;
                }
            }
            KPool.add(k);
        });
        for(int i=0;i<KPool.size();i++){
            System.out.println(KPool.get(i));
        }
    }
    public void countTPM(){
        RCaller caller = new RCaller();
        caller.setRscriptExecutable("/usr/local/bin/R");
        RCode rCode = new RCode();
        rCode.addRCode("library(cluster);library(factoextra);library(vegan);results=kmeans(");
        rCode.addRCode(""+",");
        for(int i=1;i<11;i++){
            rCode.addRCode(i+");dis = dist(a");
            rCode.addRCode(""+")^2;");
        }
        rCode.addRCode("sil = silhouette (results$cluster, dis);sum=summary(sil);sh<-sum$avg.width");
        caller.setRCode(rCode);
        caller.runAndReturnResult("sh");
    }
    public void findTriad(){
        String inputFileS="/Users/xujun/Desktop/eQTL/total/DESeq/TaS-D/750MSiPAS-TPM.txt";
        RowTable rt = new RowTable(inputFileS);
        HashMap AExpre = new HashMap();
        HashMap BExpre = new HashMap();
        HashMap DExpre = new HashMap();
        for (int i=0;i<rt.getRowNumber();i++){
            if(rt.getCellAsDouble(i, 1)>0.3 ){
                if(rt.getCellAsString(i, 0).charAt(8)=='A'){
                    AExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                }else{
                    if(rt.getCellAsString(i, 0).charAt(8)=='B'){
                        BExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                    }else{
                        DExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                    }
                }
            }
        }
        Iterator iterA = AExpre.entrySet().iterator();
        while (iterA.hasNext()) {
            HashMap.Entry entryA = (HashMap.Entry)iterA.next();
            Object keyA = entryA.getKey();
            Object valA = entryA.getValue();
            Iterator iterB = BExpre.entrySet().iterator();
            while(iterB.hasNext()){
                HashMap.Entry entryB = (HashMap.Entry)iterB.next();
                Object keyB = entryB.getKey();
                if(keyB.toString().charAt(7)==keyA.toString().charAt(7)){
                    Object valB = entryB.getValue();
                    if(isEqual(Double.parseDouble(valA.toString()),Double.parseDouble(valB.toString()))){
                        Iterator iterD = DExpre.entrySet().iterator();
                        HashMap.Entry entryD = (HashMap.Entry)iterD.next();
                        Object keyD = entryD.getKey();
                        if(keyA.toString().charAt(7)==keyD.toString().charAt(7)){
                            Object valD = entryD.getValue();
                            while(iterD.hasNext()){
                                if(isEqual(Double.parseDouble(valA.toString()),Double.parseDouble(valD.toString()))){
                                    System.out.println(keyA.toString()+"\n"+keyB.toString()+"\n"+keyD.toString());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    public boolean isEqual(double a,double b){
        double eps = 1E-1;
        if (Math.abs(a - b) / a <= eps) {
            return true;
        }else{
            return false;
        }
    }
    public void splitSample(){
//         String inputDirS="/data2/junxu/SiPASData/20200324/Rawdata/20190118aptest/";
//         String outputDirS ="/data1/home/junxu/wheat/SiPAS1.1/testR/subFastqs";
//         RowTable rt = new RowTable("/data1/home/junxu/wheat/SiPAS1.1/testR/SampleInformation_test.txt");
        String inputDirS="/data2/junxu/SiPASData/20200324/Rawdata/20190118S1-7/";
        String outputDirS ="/data1/home/junxu/eQTL/7phenotype/subFastqs";
        RowTable rt = new RowTable("/data1/home/junxu/eQTL/7phenotype/SampleInformation.txt");
        HashMap barcodeName = new HashMap();
        List<String> nameList = new ArrayList();
        for (int i=0;i<rt.getRowNumber();i++){
//             barcodeName.put(rt.getCellAsString(i, 2), rt.getCellAsString(i, 0)+"-"+rt.getCell(i,3));
//             nameList.add(rt.getCellAsString(i, 0)+"-"+rt.getCell(i,3));
            barcodeName.put(rt.getCellAsString(i, 2), rt.getCellAsString(i, 0)+"_"+rt.getCell(i,1));
            nameList.add(rt.getCellAsString(i, 0)+"_"+rt.getCell(i,1));
        }
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        List<File> fList = Arrays.asList(fs);
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS, p+"_R1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"_R2.fq.gz").getAbsolutePath();
                BufferedWriter [] bw = new BufferedWriter[nameList.size()];
                BufferedWriter [] bw1 = new BufferedWriter[nameList.size()];
                BufferedWriter bwE =IOUtils.getTextGzipWriter(new File(outputDirS,"error.txt.gz").getAbsolutePath());
                for (int i = 0; i < nameList.size(); i++) {
                    bw[i]=IOUtils.getTextGzipWriter(new File(outputDirS, nameList.get(i)+"_R1.fq.gz").getAbsolutePath());
                    bw1[i]=IOUtils.getTextGzipWriter(new File(outputDirS, nameList.get(i)+"_R2.fq.gz").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                String temp = null;String seq=null;
                String currentBarcode=null;
                while ((temp = br.readLine()) != null) {
                    seq=br.readLine();
                    if(barcodeName.get(seq.substring(0, 8))!=null){
                        int index = nameList.indexOf(barcodeName.get(seq.substring(0, 8)));
                        bw[index].write(temp);bw[index].newLine();
                        bw[index].write(seq);bw[index].newLine();
                        bw[index].write(br.readLine());bw[index].newLine();
                        bw[index].write(br.readLine());bw[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                    }else{
                        br.readLine();br.readLine();
                        bwE.write(br1.readLine());bwE.newLine();
                        bwE.write(br1.readLine());bwE.newLine();
                        bwE.write(br1.readLine());bwE.newLine();
                        bwE.write(br1.readLine());bwE.newLine();
                    }
                }
                br.close();br1.close();
                for(int i=0;i<bw.length;i++){
                    bw[i].flush();bw[i].close();
                    bw1[i].flush();bw1[i].close();
                }

            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("finished");
    }
    public void avergeTranscriptLengthInWheat(){
        String inputFile="/Users/xujun/Desktop/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3";
        GeneFeature gf = new GeneFeature(inputFile);
        long totalTranscriptLength = 0;
        int transcriptNumber = 0;
        for(int i=0;i<gf.genes.length;i++){
            for(int j=0;j<gf.genes[i].ts.size();j++){
                totalTranscriptLength+=gf.genes[i].ts.get(j).getTranscriptEnd()-gf.genes[i].ts.get(j).getTranscriptStart()+1;
                transcriptNumber++;
            }
        }
        System.out.println(totalTranscriptLength/transcriptNumber);
    }
    private void removeUnmappedStartEnd () {
        List<String> diffValumnMethodList = new ArrayList<String>();
        String inputDirS = new File("/data1/home/junxu/analysis0215/test-result/test-withoutERCC").getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        for(int i =0;i<fList.size();i++){
            diffValumnMethodList.add(fList.get(i).getName().replace(".sam", ""));
        }
        Collections.sort(diffValumnMethodList);
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(inputDirS+"/change"+f.getName());
                String temp = null;
                String [] tem = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("[")) continue;
                    if(temp.startsWith("@")){
                        bw.write(temp);bw.newLine();
                    }else{
                        List<String> tList= FStringUtils.fastSplit(temp);
                        tem = tList.toArray(new String[tList.size()]);
                        if(!tem[2].equals("*")){
                            bw.write(temp);bw.newLine();
                        }
                    }
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    public void findIndex(){
        String inputFile="/data1/home/junxu/wheat/splitByIndex/indexSiPAS-all.fq";
        List <String> indexList=new ArrayList<>();
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw =IOUtils.getTextGzipWriter(inputFile.replace(".fq", ".txt"));
            String temp=null;
            while((temp = br.readLine()) != null){
                String index=temp.split(" ")[1].split(":")[3];
                if(!indexList.contains(index)){
                    indexList.add(index);
                }
            }
            BufferedReader br1 = IOUtils.getTextReader(inputFile);
            int [] indexNumber = new int[indexList.size()];
            while((temp = br1.readLine()) != null){
                String index=temp.split(" ")[1].split(":")[3];
                indexNumber[indexList.indexOf(index)]++;
            }
            for(int i=0;i<indexList.size();i++){
                bw.write(indexList.get(i)+"\t"+indexNumber[i]);
                bw.newLine();
            }
            br.close();br1.close();
            bw.flush();bw.close();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
