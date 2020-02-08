/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;



import static com.sun.javafx.animation.TickCalculation.sub;
import com.sun.javafx.geom.Point2D;
import pgl.infra.table.RowTable;
//import static htsjdk.samtools.util.SequenceUtil.a;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import org.apache.commons.lang.StringUtils;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.FractionalSimilarityScorer;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.phylo.DistanceMatrixCalculator;
import rcaller.RCaller;
import rcaller.RCode;
import pgl.infra.utils.IOUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;
/**
 *
 * @author kanglipeng
 /逢10无法count*/
public class SiPASHomo {

    public SiPASHomo() {
//        System.out.println(this.getGeneSequenceNew(29, 206100, 206160));
//        System.out.println(this.getGeneSequence("/Users/xujun/Desktop/IGVmaterial/wheat/abd_iwgscV1.fa",1,98019574,98020000));
//        System.out.println(this.getGeneSequenceNew(1,98019574,98020000));
//        this.homoAllignmentOld();
//        this.withinAlignment();
//        this.pairwiseDistance();
//        this.KPattern();
//        this.getFQ2();
//        this.geneDensity();
        this.homoAllignmentNew();

    }
    public void geneDensity(){
        String geneFile="/data1/home/junxu/SiPAS-homo/homologyGene.txt";
        String outputDir="/data1/home/junxu/SiPAS-homo/";
        BufferedWriter[] bw =new BufferedWriter[3];
        for(int i=0;i<3;i++){
            bw[i]=IOUtils.getTextWriter(outputDir+i+".txt");
        }
        BufferedWriter bwB =IOUtils.getTextWriter("/data1/home/junxu/SiPAS-homo/homoBDensity.txt");
        BufferedWriter bwD =IOUtils.getTextWriter("/data1/home/junxu/SiPAS-homo/homoDDensity.txt");
        RowTable rt = new RowTable("/data1/home/junxu/SiPAS-homo/readme.txt");
         HashMap subChrLength=new HashMap();
         for( int i=0;i< rt.getRowNumber();i++){
             subChrLength.put(rt.getCellAsInteger(i, 0), rt.getCellAsInteger(i, 2));
         }
//        BufferedWriter bw1 =IOUtils.getTextWriter(outputFile1);
//        GeneFeature gf =new GeneFeature("/Users/xujun/Desktop/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        GeneFeature gf =new GeneFeature("/data1/home/junxu/SiPAS-homo/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        String tempF = null;String tempG=null;String[] tem = null;
        try {
            for(int i =0;i<3;i++){
                BufferedReader brG=IOUtils.getTextReader(geneFile);
                while((tempG=brG.readLine())!=null){
                    tem=tempG.split("\t");
                    if(Integer.valueOf(gf.getGeneChromosome(gf.getGeneIndex(tem[i])))%2!=0 ){
                        int chrLength=Integer.valueOf(subChrLength.get(gf.getGeneChromosome(gf.getGeneIndex(tem[i]))).toString())+Integer.valueOf(subChrLength.get(gf.getGeneChromosome(gf.getGeneIndex(tem[i]))+1).toString());
                        double pos = (double)Integer.valueOf(gf.getGeneStart(gf.getGeneIndex(tem[i])))/chrLength;
                        bw[i].write(pos+"");bw[i].newLine();
                    }else{
                        int chrLength=Integer.valueOf(subChrLength.get(gf.getGeneChromosome(gf.getGeneIndex(tem[i]))).toString())+Integer.valueOf(subChrLength.get(gf.getGeneChromosome(gf.getGeneIndex(tem[i]))-1).toString());    
                        double pos =(double) (Integer.valueOf(gf.getGeneStart(gf.getGeneIndex(tem[i])))+Integer.valueOf(subChrLength.get(gf.getGeneChromosome(gf.getGeneIndex(tem[i]))-1).toString()))/chrLength;
    //                        geneNumber[Integer.valueOf(tem[0])-1][pos]++;
                            bw[i].write(pos+"");bw[i].newLine();
                    }
                } 
                brG.close();bw[i].flush();bw[i].close();
            }
            
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
//        
    }
    public void getFQ2(){
        String inputDirS="/data1/home/junxu/SiPAS-homo/A/FQ"; 
        String outputDirS="/data1/home/junxu/SiPAS-homo/A/FQ2";
        File[] fs= new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".fq");
        for(int i=0;i<fs.length;i++){
            if(!fs[i].getName().contains("temp")){
                fList.add(fs[i]);
            }
        }
        fList.stream().forEach(f -> {
            try{
                BufferedReader br=IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw=IOUtils.getTextWriter(outputDirS+"/"+f.getName());
                String temp=null;
                while((temp=br.readLine())!=null){
                    if(temp.endsWith("2")){
                        bw.write(temp);bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                    }
                }
                br.close();bw.flush();bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }     
        });
        
    }
    public void KPattern(){
        String[] inputDirS=new String [3];
        inputDirS[0]="/data1/home/junxu/SiPAS-homo/A/FQ2";
        inputDirS[1]="/data1/home/junxu/SiPAS-homo/B/FQ2";
        inputDirS[2]="/data1/home/junxu/SiPAS-homo/D/FQ2";
        String outputDirS ="/data1/home/junxu/SiPAS-homo/kValue/3002";
        String geneFile="/data1/home/junxu/SiPAS-homo/homologyGene.txt";
        File[] fs1 = new File(inputDirS[0]).listFiles();
        File[] fs2 = new File(inputDirS[1]).listFiles();
        File[] fs3 = new File(inputDirS[2]).listFiles();
        HashSet<String> nameSet1 = new HashSet();
        HashSet<String> nameSet2 = new HashSet();
        HashSet<String> nameSet3 = new HashSet();
        fs1 = IOUtils.listFilesEndsWith(fs1, ".fq");
        fs2 = IOUtils.listFilesEndsWith(fs2, ".fq");
        fs3 = IOUtils.listFilesEndsWith(fs3, ".fq");
        for(int i=0;i<fs1.length;i++){
           if (fs1[i].length()<=34900) continue;
           nameSet1.add(fs1[i].getName().split("\\.")[0]);
        }
        for (int i = 0; i < fs2.length; i++) {
//            if (fs2[i].length()<=70000) continue;//changable
            if (fs2[i].length()<=34900) continue;
            nameSet2.add(fs2[i].getName().split("\\.")[0]);
        }
        for (int i = 0; i < fs3.length; i++) {
            if (fs3[i].length()<=34900) continue;
            nameSet3.add(fs3[i].getName().split("\\.")[0]);
        }
        String temp =null;int n=0;
        try{
            BufferedReader br=IOUtils.getTextReader(geneFile);
            while((temp=br.readLine())!=null){
                n++;
                if(nameSet1.contains(temp.split("\t")[0]) && nameSet2.contains(temp.split("\t")[1]) && nameSet3.contains(temp.split("\t")[2])){
                    for(int i=0;i<3;i++){
                        File dir = new File(new File (outputDirS).getAbsolutePath());
                        StringBuilder sb = new StringBuilder();
                        sb.append("sed -n '1,400p' ").append(new File(inputDirS[i],temp.split("\t")[i]+".fq").getAbsolutePath());
                        sb.append(" > ").append( "temp"+i+".fq");
                        String command = sb.toString();
        //                System.out.println(command);  
                        String [] cmdarry ={"/bin/bash","-c",command};
                        Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                        p.waitFor();p.destroy();
                    }
                    File dir = new File(new File (outputDirS).getAbsolutePath());
                    StringBuilder sb1 = new StringBuilder();
                    sb1.append("cat ").append(new File(outputDirS,"temp0.fq").getAbsolutePath()+" ");
                    sb1.append(new File(outputDirS,"temp1.fq").getAbsolutePath()+" ").append(new File(outputDirS,"temp2.fq").getAbsolutePath());
                    sb1.append(" > ").append( n+".fq");
                    String command1 = sb1.toString();
            //        System.out.println(command);  
                    String [] cmdarry1 ={"/bin/bash","-c",command1};
                    Process p1=Runtime.getRuntime().exec(cmdarry1,null,dir);
                    p1.waitFor();p1.destroy();
                }  
            }
            
        }catch (Exception e) {
            e.printStackTrace();
        }
        File[] fs = new File(outputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".fq");
        for(int i=0;i<fs.length;i++){
            if(!fs[i].getName().contains("temp")){
                fList.add(fs[i]);
            }
        }
        fList.stream().forEach(f -> {
            try {
                BufferedReader br=IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw=IOUtils.getTextWriter(outputDirS+"/pairwiseValue/"+f.getName().replace(".fq", ".txt"));
                String tempV=null;String temdV=null;int j=0;
                SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
                SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
                SubstitutionMatrix<NucleotideCompound> M = SubstitutionMatrixHelper.getNuc4_2();
                SequencePair<DNASequence, NucleotideCompound> psa = null;
                DNASequence query=new DNASequence("");
                DNASequence hit=new DNASequence("");
                while((tempV=br.readLine())!=null){
                    j++;
                    query=new DNASequence(br.readLine());
                    br.readLine();br.readLine();
                    int m=0;
                    BufferedReader br2=IOUtils.getTextReader(f.getAbsolutePath());
                    while((temdV=br2.readLine())!=null){
                        m++;
                        if(m>j){
                            hit=new DNASequence(br2.readLine());
                            br2.readLine();br2.readLine();
                            psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                            if(psa.getNumSimilars()>=30){
                                int value=0;
                                for(int i=0;i<psa.getLength();i++){
                                    if(psa.getQuery().toString().charAt(i)!=psa.getTarget().toString().charAt(i)){
                                        value++;
                                    }
                                }
                                bw.write(String.valueOf(value));bw.newLine();
                            }
                        }else{
                            br2.readLine();br2.readLine();br2.readLine();
                        }
                    }
                    br2.close();
                }
                bw.flush();bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
//        String inputfsK=outputDirS+"/pairwiseValue";
//        File[] fsK = new File(inputfsK).listFiles();
//        List<File> fListK = new ArrayList(Arrays.asList());
//        fsK = IOUtils.listFilesEndsWith(fsK, ".txt");
//        for(int i=0;i<fsK.length;i++){
//            if(!fsK[i].getName().contains("temp")){
//                fListK.add(fsK[i]);
//            }
//        }
//        RCaller callerL = new RCaller();
//        callerL.setRscriptExecutable("/usr/local/bin/R");
//        RCode rCodeL = new RCode();
//        rCodeL.addRCode("library(cluster);library(factoextra);library(vegan);");
//        callerL.setRCode(rCodeL);
//        callerL.runOnly();
//        List KPool =new ArrayList();
//        fListK.stream().forEach(f -> {
//            double sh=0;double shM=0;int k=0;
//            for(int i=1;i<11;i++){
//                RCaller caller = new RCaller();
//                caller.setRscriptExecutable("/usr/local/bin/R");
//                RCode rCode = new RCode();
//                rCode.addRCode("results=kmeans("+f+",");
//                rCode.addRCode(i+");dis = dist(");
//                rCode.addRCode(f+")^2;");
//                rCode.addRCode("sil = silhouette (results$cluster, dis);sum=summary(sil);sh<-sum$avg.width");
//                caller.setRCode(rCode);
//                caller.runAndReturnResult("sh");
//                sh=caller.getParser().getAsDoubleArray("sh")[0];
//                if(sh>shM){
//                    shM=sh;
//                }else{
//                    k=i;
//                    break;
//                }
//            }
//            KPool.add(k);
//        });
     }
    public int distancesum(int x[], int y[], int n) { 
        int sum = 0; 
        for (int i = 0; i < n; i++) 
            for (int j = i + 1; j < n; j++) 
                sum += (Math.abs(x[i] - x[j]) +  Math.abs(y[i] - y[j])); 
        return sum; 
    } 
    public void pairwiseDistance(){
        String inputFile="/data1/home/junxu/SiPAS-homo/cluster/pool.fq";
//        String inputFile="/Users/xujun/Desktop/eQTL/total/Homology/pool.fq";
        BufferedReader br1=IOUtils.getTextReader(inputFile);
        BufferedWriter bw=IOUtils.getTextWriter("/data1/home/junxu/SiPAS-homo/cluster/parsepool.fq");
        BufferedWriter bww=IOUtils.getTextWriter("/data1/home/junxu/SiPAS-homo/cluster/pairwiseValue.txt");
//        BufferedWriter bw=IOUtils.getTextWriter("/Users/xujun/Desktop/eQTL/total/Homology/parsepool.fq");
        String temp=null;String temd=null;int n=0;
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        SubstitutionMatrix<NucleotideCompound> M = SubstitutionMatrixHelper.getNuc4_2();
        SequencePair<DNASequence, NucleotideCompound> psa = null;
        try{
            DNASequence query=new DNASequence("");
            DNASequence hit=new DNASequence("");
            while((temp=br1.readLine())!=null){
                n++;
                query=new DNASequence(br1.readLine());
                br1.readLine();br1.readLine();
                int m=0;
                BufferedReader br2=IOUtils.getTextReader(inputFile);
                while((temd=br2.readLine())!=null){
                    m++;
                    if(m>n){
                        hit=new DNASequence(br2.readLine());
                        br2.readLine();br2.readLine();
                        psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                        if(psa.getNumSimilars()>=50){
//                            bw.write(psa.getQuery().toString());bw.newLine();
//                            bw.write(psa.getTarget().toString());bw.newLine();
//                            int [] x=new int [psa.getLength()];int [] y=new int [psa.getLength()];
//                            String te=null;
//                            te=psa.getQuery().toString().replace("_", "0");
//                            te=te.replace("A", "1");te=te.replace("T", "2");
//                            te=te.replace("C", "3");te=te.replace("G", "4");
//                            for(int i=0;i<te.length();i++){
//                                x[i]=te.charAt(i);
//                            }
//                            te=psa.getTarget().toString().replace("_", "0");
//                            te=te.replace("A", "1");te=te.replace("T", "2");
//                            te=te.replace("C", "3");te=te.replace("G", "4");
//                            for(int i=0;i<te.length();i++){
//                                y[i]=te.charAt(i);
//                            }
//                            bww.write(String.valueOf(this.distancesum(x, y, psa.getLength())));bww.newLine();
                            int value=0;
                            for(int i=0;i<psa.getLength();i++){
                                if(psa.getQuery().toString().charAt(i)!=psa.getTarget().toString().charAt(i)){
                                    value++;
                                }
                            }
                            bww.write(String.valueOf(value));bww.newLine();
                        }
                    }else{
                        br2.readLine();br2.readLine();br2.readLine();
                    }
                }
                br2.close();
            }
            bw.flush();bw.close();bww.flush();bww.close();
            br1.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        

    }

    public void homoAllignmentOld(){
        String faFile="/Users/xujun/Desktop/IGVmaterial/wheat/abd_iwgscV1.fa";
//        System.out.println(this.getGeneSequence(faFile, 37,218264266,218265670));
//        String geneFile="/Users/xujun/Desktop/eQTL/total/Homology/homologyGene.txt";
//        String outputFile="/Users/xujun/Desktop/eQTL/total/homology/homoAB.txt"; 
//        String outputFile1="/Users/xujun/Desktop/eQTL/total/homo1.txt";
//        String faFile="/data1/home/junxu/SiPAS-homo/abd_iwgscV1.fa";
        String geneFile="/data1/home/junxu/SiPAS-homo/homologyGene.txt";
//        String outputFile="/data1/home/junxu/SiPAS-homo/AB/homoAB.txt";
        String outputFile="/data1/home/junxu/SiPAS-homo/BD/Local/homoBDlocal.txt";
//        String outputFile1="/data1/home/junxu/SiPAS-homo/homo1.txt";
//        System.out.println(this.getGeneSequence(faFile,60,62 ));
        BufferedReader brG=IOUtils.getTextReader(geneFile);
        BufferedWriter bw =IOUtils.getTextWriter(outputFile);
//        BufferedWriter bw1 =IOUtils.getTextWriter(outputFile1);
//        GeneFeature gf =new GeneFeature("/Users/xujun/Desktop/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        GeneFeature gf =new GeneFeature("/data1/home/junxu/SiPAS-homo/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        String tempF = null;String tempG=null;String[] tem = null;
        int[] m=new int[2];
        int[] n=new int[2];
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        SequencePair<DNASequence, NucleotideCompound> psa = null;
        int line=0;
        try {
            DNASequence query=new DNASequence("");
            DNASequence hit=new DNASequence("");
            DNASequence queryTemp=new DNASequence("");
            while((tempG=brG.readLine())!=null){
                line++;
                tem=tempG.split("\t");
                m[0]=gf.getGeneStart(gf.getGeneIndex(tem[1]));
                n[0]=gf.getGeneEnd(gf.getGeneIndex(tem[1]));
                m[1]=gf.getGeneStart(gf.getGeneIndex(tem[2]));
                n[1]=gf.getGeneEnd(gf.getGeneIndex(tem[2]));  
                
                BufferedReader[] brF=new BufferedReader[2];
                for(int i=0;i<brF.length;i++){
                    brF[i]=IOUtils.getTextReader(faFile);
                }
                for(int i=0;i<m.length;i++){
                    int pos=0;
                    while ((tempF = brF[i].readLine()) != null) {
                        if (tempF.startsWith(">")) {
                            continue;
                        }
                        if (pos==(int)(60*Math.floor(m[i]/60))) {
                            if(n[i]-m[i]>60){
                                String sequence = tempF.substring((int) (m[i]-60*Math.floor(m[i]/60)));
                                queryTemp=new DNASequence(sequence);
    //                            System.out.println(StringUtils.repeat("N",60-sequence.length())+sequence);
                                pos=pos+60;
                                while (pos < n[i]) {
                                    queryTemp=new DNASequence(queryTemp.toString()+brF[i].readLine());
    //                                System.out.println(brF.readLine());
                                    pos=pos+60;
                                }
                                if(pos>=n[i]){
                                    queryTemp=new DNASequence(queryTemp.toString()+brF[i].readLine().substring(0, (int) (n[i]-60*Math.floor(n[i]/60))+1));
    //                                System.out.println(brF.readLine().substring(0, (int) (n-60*Math.floor(n/60))+1));
                                    if(i==0){
                                        query=queryTemp;
                                        brF[0].close();
                                        break;
                                    }else{
                                        hit=queryTemp;
                                        brF[1].close();
                                        break;
                                    }
                                }
                            }else{
                                queryTemp=new DNASequence(queryTemp.toString()+tempF.substring((int) (m[i]-60*Math.floor(m[i]/60)),n[i]-pos));
    //                            String sequence = tempF.substring((int) (m-60*Math.floor(m/60)),n-pos);
    //                            System.out.println(sequence);
                                if(i==0){
                                    query=queryTemp;
                                    brF[0].close();
                                    break;
                                }else{
                                    hit=queryTemp;
                                    brF[1].close();
                                    break;
                                }
                            }

                        }else{pos=pos+60;}
                    }
                }  
                if(query.getLength()>=25000 || hit.getLength()>=25000){continue;}
//                Alignments.getAllPairsAlignments(sequences, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix)
                psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
//                System.out.println(psa.toString(60));
                int gapNumber=0;
                for(int i=1;i<psa.getLength()+1;i++){
                    if(psa.hasGap(i)){
                        gapNumber++;
                    }
                } 
//                System.out.println(psa.getNumSimilars()+" "+psa.getLength());
                System.out.println(1-(double)(psa.getNumSimilars()+gapNumber)/(double)(psa.getLength()));
                bw.write(String.format("%f",1-(double)(psa.getNumSimilars()+gapNumber)/(double)(psa.getLength())));
                bw.newLine();
//                if(line>50){break;}
            }  
            brG.close();
            bw.flush();bw.close();
//            bw1.flush();bw1.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(psa.getQuery().getLength());
            System.out.println(psa.getTarget().getLength());
            System.exit(1);
        }
//        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
//        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
//        try {
//        DNASequence query = new DNASequence ("AACATTAACATCCGAGCAATTGAAAAGGACATTTCC");
//        DNASequence hit = new DNASequence ("AACCACAGATATACTCAAATATGAGCCATATT");
////        DNASequence hit2 = new DNASequence ("AACATTAACATCCGAGCAATTGAAAAGGACATTT");
//        SequencePair<DNASequence, NucleotideCompound> psa = null;
//        psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.GLOBAL, gapPen, subMatrix);
////        System.out.println(psa.toString(60));
//        int gapNumber=0;
//        for(int i=1;i<psa.getLength()+1;i++){
//            if(psa.hasGap(i)){
//                gapNumber++;
//            }
//        }
//        System.out.println((double)(psa.getLength()-psa.getNumSimilars()-gapNumber)/(psa.getLength()-gapNumber));
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
    }
    public void homoAllignmentNew(){
//        String faFile="/Users/xujun/Desktop/IGVmaterial/wheat/abd_iwgscV1.fa";
//        String geneFile="/Users/xujun/Desktop/eQTL/total/SiPAS-homo/homologyGene.txt";
//        String outputFile="/Users/xujun/Desktop/eQTL/total/SiPAS-homo/homo.txt"; 
//        String outputFile1="/Users/xujun/Desktop/eQTL/total/homo1.txt";
        String faFile="/data1/home/junxu/SiPAS-homo/abd_iwgscV1.fa";
        String geneFile="/data1/home/junxu/SiPAS-homo/homologyGene.txt";
        String outputFile="/data1/home/junxu/SiPAS-homo/AD/homoAD-seq.txt";
//        String outputFile="/data1/home/junxu/SiPAS-homo/AD/Local/homoABlocal.txt";
        BufferedReader brG=IOUtils.getTextReader(geneFile);
        BufferedWriter bw =IOUtils.getTextWriter(outputFile);
//        BufferedWriter bw1 =IOUtils.getTextWriter(outputFile1);
//        GeneFeature gf =new GeneFeature("/Users/xujun/Desktop/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        GeneFeature gf =new GeneFeature("/data1/home/junxu/SiPAS-homo/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        String tempF = null;String tempG=null;String[] tem = null;
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        SequencePair<DNASequence, NucleotideCompound> psa = null;
        int queryIndex=0;int hitIndex=2;//need change
        try {
            while((tempG=brG.readLine())!=null){
                tem=tempG.split("\t");
                String queryS=getGeneSequenceNew(gf.getGeneChromosome(gf.getGeneIndex(tem[queryIndex])),gf.getGeneStart(gf.getGeneIndex(tem[queryIndex])),gf.getGeneEnd(gf.getGeneIndex(tem[queryIndex])));
                String hitS=getGeneSequenceNew(gf.getGeneChromosome(gf.getGeneIndex(tem[hitIndex])),gf.getGeneStart(gf.getGeneIndex(tem[hitIndex])),gf.getGeneEnd(gf.getGeneIndex(tem[hitIndex]))); 
                DNASequence query=new DNASequence(queryS);
                DNASequence hit=new DNASequence(hitS);  
                if(query.getLength()>=25000 || hit.getLength()>=25000){continue;}
//                Alignments.getAllPairsAlignments(sequences, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix)
                psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.GLOBAL, gapPen, subMatrix);
//                psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
//                System.out.println(psa.toString(60));
                int gapNumber=0;
                for(int i=1;i<psa.getLength()+1;i++){
                    if(psa.hasGap(i)){
                        gapNumber++;
                    }
                }
                FractionalSimilarityScorer<DNASequence, NucleotideCompound> scorer = new FractionalSimilarityScorer<>(psa);
                Double score = scorer.getScore();
//                if((psa.getLength()-psa.getNumSimilars()-gapNumber)/(psa.getLength()-gapNumber)==1){
//                    System.out.println(psa.getQuery().toString());
//                    System.out.println(psa.getTarget().toString());
//                    bw1.write(psa.getQuery().toString());
//                    bw1.newLine();
//                }
//                System.out.println(String.format("%f",1-(double)(psa.getNumSimilars()+gapNumber)/(double)(psa.getLength())));
                bw.write(String.format("%f",1-(double)(psa.getNumSimilars()+gapNumber)/(double)(psa.getLength()))+"\t"+score);
                if(1-(double)(psa.getNumSimilars()+gapNumber)/(double)(psa.getLength())> 0.1){
                    bw.write(" "+tempG+" ");
                    bw.write(gf.getGeneChromosome(gf.getGeneIndex(tem[queryIndex]))+" "+gf.getGeneStart(gf.getGeneIndex(tem[queryIndex]))+" "+gf.getGeneEnd(gf.getGeneIndex(tem[queryIndex]))+" ");
                    bw.write(String.format("%d",gf.getGeneEnd(gf.getGeneIndex(tem[queryIndex]))-gf.getGeneStart(gf.getGeneIndex(tem[queryIndex])))+" ");
                    bw.write(gf.getGeneChromosome(gf.getGeneIndex(tem[hitIndex]))+" "+gf.getGeneStart(gf.getGeneIndex(tem[hitIndex]))+" "+gf.getGeneEnd(gf.getGeneIndex(tem[hitIndex]))+" ");
                    bw.write(String.format("%d",gf.getGeneEnd(gf.getGeneIndex(tem[hitIndex]))-gf.getGeneStart(gf.getGeneIndex(tem[hitIndex])))+"");
                }
                bw.newLine();
                bw.write(psa.toString());
                bw.newLine();
//                System.out.println((double)(psa.getLength()-psa.getNumSimilars()-gapNumber)/(psa.getLength()-gapNumber));
            }  
            brG.close();
            bw.flush();bw.close();
//            bw1.flush();bw1.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
//        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
//        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
//        try {
//        DNASequence query = new DNASequence ("AACATTAACATCCGAGCAATTGAAAAGGACATTTCC");
//        DNASequence hit = new DNASequence ("AACCACAGATATACTCAAATATGAGCCATATT");
////        DNASequence hit2 = new DNASequence ("AACATTAACATCCGAGCAATTGAAAAGGACATTT");
//        SequencePair<DNASequence, NucleotideCompound> psa = null;
//        psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.GLOBAL, gapPen, subMatrix);
////        System.out.println(psa.toString(60));
//        int gapNumber=0;
//        for(int i=1;i<psa.getLength()+1;i++){
//            if(psa.hasGap(i)){
//                gapNumber++;
//            }
//        }
//        System.out.println((double)(psa.getLength()-psa.getNumSimilars()-gapNumber)/(psa.getLength()-gapNumber));
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
    }
    public void withinAlignment(){
//        String faFile="/Users/xujun/Desktop/IGVmaterial/wheat/abd_iwgscV1.fa";
//        String geneFile="/Users/xujun/Desktop/eQTL/total/SiPAS-homo/homologyGene.txt";
//        String outputFile="/Users/xujun/Desktop/eQTL/total/SiPAS-homo/A.txt"; 
//        String outputFile1="/Users/xujun/Desktop/eQTL/total/homo1.txt";
        String faFile="/data1/home/junxu/SiPAS-homo/abd_iwgscV1.fa";
        String geneFile="/data1/home/junxu/SiPAS-homo/homologyGene.txt";
//        String outputDir="/data1/home/junxu/SiPAS-homo/D";//need change
        String outputDir="/data1/home/junxu/SiPAS-homo/D/Local";//need change
//        String outputFile1="/data1/home/junxu/SiPAS-homo/homo1.txt";
        BufferedReader brG=IOUtils.getTextReader(geneFile);
//        BufferedWriter bw =IOUtils.getTextWriter(outputDir+"/D.txt");//need change
        BufferedWriter bw =IOUtils.getTextWriter(outputDir+"/Dlocal.txt");
//        BufferedWriter bw1 =IOUtils.getTextWriter(outputFile1);
//        GeneFeature gf =new GeneFeature("/Users/xujun/Desktop/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        GeneFeature gf =new GeneFeature("/data1/home/junxu/SiPAS-homo/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        String tempF = null;String tempG=null;String[] tem = null;String tempFQ=null;
        int gh=2;//hit
        int a=0;int b=0;
        int chr=-1;
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        SequencePair<DNASequence, NucleotideCompound> psa = null; 
        try {
            DNASequence query=new DNASequence("");
            DNASequence hit=new DNASequence("");
            DNASequence queryTemp =new DNASequence("");
            while((tempG=brG.readLine())!=null){
                tem=tempG.split("\t");
                BufferedReader brF=IOUtils.getTextReader(faFile);
                chr=gf.getGeneChromosome(gf.getGeneIndex(tem[gh]));
                if(gf.getGeneStrand(gf.getGeneIndex(tem[gh]))==1){//1是正链，0是负链
                    a=gf.getGeneEnd(gf.getGeneIndex(tem[gh]))-500;
                    b=gf.getGeneEnd(gf.getGeneIndex(tem[gh]))+500;
                }else{
                    a=gf.getGeneStart(gf.getGeneIndex(tem[gh]))-500;
                    b=gf.getGeneStart(gf.getGeneIndex(tem[gh]))+500;
                }
                File dir = new File(new File (outputDir+"/BAM/").getAbsolutePath());
                StringBuilder sb = new StringBuilder();
                sb.append("samtools view -h ").append("/data1/home/junxu/wheat/doubleAll/sams/5000000SiPAS-10am1sorted.bam ");
                sb.append(chr+":"+a+"-"+b+" -o ").append(tem[gh]+".sorted.bam");
                String command = sb.toString();
//                System.out.println(command);  
                String [] cmdarry ={"/bin/bash","-c",command};
                Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                p.waitFor();p.destroy();
                File dir1 = new File(new File (outputDir+"/FQ/").getAbsolutePath());
                StringBuilder sb1 = new StringBuilder();
                sb1.append("samtools bam2fq ");
                sb1.append(dir.getAbsolutePath()+"/"+tem[gh]+".sorted.bam ").append("> ");
                sb1.append(tem[gh]+ ".fq");
                String command1 = sb1.toString();
//                System.out.println(command1);  
                String [] cmdarry1 ={"/bin/bash","-c",command1};
                try{
                    Process p1=Runtime.getRuntime().exec(cmdarry1,null,dir1);
                    p1.waitFor();p1.destroy();
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                int m=0;int n=0;
                int gq=2;//query 
                chr=gf.getGeneChromosome(gf.getGeneIndex(tem[gq]));
                if(gf.getGeneStrand(gf.getGeneIndex(tem[gq]))==1){//1是正链，0是负链
                    m=gf.getGeneEnd(gf.getGeneIndex(tem[gq]))-500;
                    n=gf.getGeneEnd(gf.getGeneIndex(tem[gq]))+500;
                }else{
                    m=gf.getGeneStart(gf.getGeneIndex(tem[gq]))-500;
                    n=gf.getGeneStart(gf.getGeneIndex(tem[gq]))+500;
                }
                query=new DNASequence(this.getGeneSequenceNew(chr, m, n));
                if(new File(dir1.getAbsolutePath()+"/"+tem[gh]+".fq").length()==0){
                    continue;
                }else{
                    BufferedReader brFQ=IOUtils.getTextReader(dir1.getAbsolutePath()+"/"+tem[gh]+".fq");
                    while((tempFQ=brFQ.readLine())!=null){
                        hit=new DNASequence(brFQ.readLine());
                        brFQ.readLine();brFQ.readLine();
//                        psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.GLOBAL, gapPen, subMatrix);
                        psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
        //                System.out.println(psa.toString(60));
                        int gapNumber=0;
                        for(int i=1;i<psa.getLength()+1;i++){
                            if(psa.hasGap(i)){
                                gapNumber++;
                            }
                        }
//                        if((double)(psa.getLength()-psa.getNumSimilars()-gapNumber)/(psa.getLength()-gapNumber)!=1){
                            bw.write(String.format("%f",1-(double)(psa.getNumSimilars()+gapNumber)/(double)(psa.getLength())));
//                            System.out.println(String.format("%f",1-(double)(psa.getNumSimilars()+gapNumber)/(double)(psa.getLength())));
//                        }
                        bw.newLine();
                    }
                    brFQ.close();
                }
                
//                System.out.println((double)(psa.getLength()-psa.getNumSimilars()-gapNumber)/(psa.getLength()-gapNumber));
                
            }  
            brG.close();
            bw.flush();bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    public String getGeneSequence(String genomeDir,int chr,int m,int n) {
        BufferedReader br;
        if (genomeDir.endsWith("gz")) {
            br = IOUtils.getTextGzipReader(genomeDir);
        } else {
            br = IOUtils.getTextReader(genomeDir);
        }
        String temp = null;String temp1=null;
        int pos = 0; String queryTemp="";
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    if(Integer.parseInt(temp.substring(1))==chr){
                         while((temp = br.readLine()) != null){
                             if (pos==(int)(60*Math.floor(m/60))) {
                                if(n-m>60){
                                    String sequence = temp.substring((int) (m-60*Math.floor(m/60)));
                                    queryTemp=sequence;
            //                        System.out.println(sequence);//StringUtils.repeat("N",60-sequence.length())+
                                    pos=pos+60;
                                    while (pos < n) {
                                        queryTemp=queryTemp+br.readLine();
                                        pos=pos+60;
                                    }
                                    if(pos>=n){
            //                            System.out.println(br.readLine().substring(0, (int) (n-60*Math.floor(n/60))+1));
                                        queryTemp+=br.readLine().substring(0, (int) (n-60*Math.floor(n/60))+1);
                                        break;
                                    }
                                }else{
                                    queryTemp+=temp.substring((int)(m-60*Math.floor(m/60)),n-pos);
                                    break;
                                }
                            }else{pos=pos+60;}
                         }
                    }else{
                        continue;
                    }
                }    
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return queryTemp;
    }
    public String getGeneSequenceNew(int chr,int m,int n) {
//        String readFileDir="/data1/home/junxu/SiPAS-homo/FA";
        String readFileDir="/Users/xujun/Desktop";
//        String readFileDir="/Users/xujun/";
        BufferedReader br;
        br = IOUtils.getTextReader(readFileDir+"/"+chr+".fa");
        String temp = null;String temp1=null;
        int pos = 0; String queryTemp="";
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    continue;
                }
                if (pos==(int)(60*Math.floor(m/60))) {
                    if(n-m>60){
                        String sequence = temp.substring((int) (m-60*Math.floor(m/60)));
                        queryTemp=sequence;
            //                        System.out.println(sequence);//StringUtils.repeat("N",60-sequence.length())+
                        pos=pos+60;
                        while (pos < n) {
                            queryTemp=queryTemp+br.readLine();
                            pos=pos+60;
                        }
                        if(pos>=n){
            //                            System.out.println(br.readLine().substring(0, (int) (n-60*Math.floor(n/60))+1));
                            queryTemp+=br.readLine().substring(0, (int) (n-60*Math.floor(n/60))+1);
                            break;
                        }
                    }else{
                        queryTemp+=temp.substring((int)(m-60*Math.floor(m/60)),n-pos);
                        break;
                    }
                }else{pos=pos+60;}
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return queryTemp;
    }
//    public static void main(String[] args){
//        int x[] = { -1, 1, 3, 2 }; 
//        int y[] = { 5, 6, 5, 3 }; 
//        int n = x.length;   
//        System.out.println(distancesum(x, y, n)); 
//        new SiPASHomo();
//    }

}

