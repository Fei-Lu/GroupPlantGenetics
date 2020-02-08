/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.table.RowTable;

/**
 *
 * @author Aoyue
 */
public class TaxaInfo {

    public TaxaInfo() {
        //this.checkTaxaCategory();
        //this.taxaGroupSummary();
        //this.removeSpacing();
        //this.checkCategoryandAddGermplasm();
//        this.checktaxaList();
//        this.reversechecktaxaList();
        //this.mergeInfo();
        //this.addGroupToMDS();
        //this.selected();
        //this.selectUseful();
        //this.removemixed();
        //new B73Distance();
        //this.countTaxaInGroup();
        //this.MDSwithHighDepth();
        //this.testPCA();
        this.addGrouptoMDSmatrix();
        
        
        
        
        
        
        
    }
    
    /**
     * 将自己抽取的9236 SNPs
     */
    
    public void addGrouptoMDSmatrix (){
        String infileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/012_PCA/004_mds_matrix.txt";
        String groupFileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/000_group/geneticGroup.manual.txt";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/012_PCA/004_mds_matrix.addgroup.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(groupFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            HashMap<String,String> hm = new HashMap<>();
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                String taxa = l.get(0); String group = l.get(4);
                hm.put(taxa, group);
            }
            br.close();
            System.out.println(hm.size() + "    hashMap");
            
            br=IOUtils.getTextReader(infileS);
            temp = br.readLine();
            bw.write(  "Taxa\tGroup\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\n");
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                String taxa = l.get(0); String group = hm.get(taxa);
                StringBuilder sb = new StringBuilder();
                sb.append(taxa).append("\t").append(group);
                for(int i =2; i<l.size(); i++){
                    sb.append("\t").append(l.get(i));
                }
                bw.write(sb.toString());bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 将IBS matrix矩阵添加分组信息，进行PCA练习,按照网页上教学的方法进行测试。
     */
    public void testPCA(){
        String infileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/012_PCA/003_distance_matrix_fei.txt";
        String groupFileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/000_group/geneticGroup.manual.txt";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/012_PCA/003_distance_matrix_fei.addgroup.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(groupFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            HashMap<String,String> hm = new HashMap<>();
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                String taxa = l.get(0); String group = l.get(4);
                hm.put(taxa, group);
            }
            br.close();
            System.out.println(hm.size() + "    hashMap");
            
            br=IOUtils.getTextReader(infileS);
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                String taxa = l.get(0); String group = hm.get(taxa);
                StringBuilder sb = new StringBuilder();
                sb.append(taxa).append("\t").append(group);
                for(int i =1; i<l.size(); i++){
                    sb.append("\t").append(l.get(i));
                }
                bw.write(sb.toString());bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    public void MDSwithHighDepth(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/006_MDS_addGroupInfo.txt";
        String highDeptFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/reccesiveDeleterious_hmp321_highDepth.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/010_MDS_highDepth.txt";
        RowTable<String> t = new RowTable<>(highDeptFileS);
        Set<String> s = new HashSet<>(t.getColumn(0));
        String[] sa = s.toArray(new String[s.size()]);
        Arrays.sort(sa);
 
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String taxa = PStringUtils.fastSplit(temp).get(0);
                String group = PStringUtils.fastSplit(temp).get(7);
                int index = Arrays.binarySearch(sa, taxa);
                if(index < 0) continue;
                if(group.equals("unknown")) continue;
                if(group.equals("mixed"))continue;
                if(group.equals("popcorn"))continue;
                if(group.startsWith("sweet"))continue;
                
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    public void countTaxaInGroup(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/005_loadVSB73Distance/004_reccesiveDeleterious_merge_notraitMixedUnknow.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(3);
        Set<String> s = new HashSet<>(l);
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "    " + Collections.frequency(l, a));
        }
    }
    
    /**
     * I only want  6 groups Landrace Teosinte Tripsacum ss nss ts 
     */
     public void removemixed(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/006_MDS_addGroupInfo.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/009_MDS_selectLanTeoTrissnss_ts.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String taxon = PStringUtils.fastSplit(temp).get(0);
                String germplasm = PStringUtils.fastSplit(temp).get(6);
                String feiGroup = PStringUtils.fastSplit(temp).get(7);
                if(germplasm.equals("unknown")) continue;
                if(germplasm.equals("Improved")) {
//                    if(feiGroup.equals("teo")){
//                        bw.write(PStringUtils.fastSplit(temp).get(0));bw.write("\t");
//                        bw.write(PStringUtils.fastSplit(temp).get(1));bw.write("\t");
//                        bw.write(PStringUtils.fastSplit(temp).get(2));bw.write("\t");
//                        bw.write(PStringUtils.fastSplit(temp).get(3));bw.write("\t");
//                        bw.write(PStringUtils.fastSplit(temp).get(4));bw.write("\t");
//                        bw.write(PStringUtils.fastSplit(temp).get(5));bw.write("\t");
//                        bw.write("Teo");bw.newLine();
//                    }
                    if(feiGroup.equals("unknown"))continue;
                    if(feiGroup.equals("mixed"))continue;
                    if(feiGroup.equals("popcorn"))continue;
                    if(feiGroup.startsWith("sweet"))continue;
                    bw.write(temp);bw.newLine();
                }
                else{
                    bw.write(temp);bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
     
    public void selectUseful(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/006_MDS_addGroupInfo.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/008_MDS_selectLanTeoTrissnsstsmixed.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String taxon = PStringUtils.fastSplit(temp).get(0);
                String germplasm = PStringUtils.fastSplit(temp).get(6);
                String feiGroup = PStringUtils.fastSplit(temp).get(7);
                if(germplasm.equals("unknown")) continue;
                if(germplasm.equals("Improved")) {
                    if(feiGroup.equals("unknown"))continue;
                    bw.write(temp);bw.newLine();
                }
                else{
                     bw.write(temp);bw.newLine();
                }
                
               
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    /**
     * Germplasm有5个分类：Landrace Teosinte Tripsacum Improved unknown
     */
    public void selected(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/005_taxa1210GroupInfo_mergeFeiGroupInfo_removeTIP498.2.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/007_MDS_selectLanTeoTri.txt";
        String a = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/MDS_PCs_Matrix_68300sites.txt";
        RowTable<String> t = new RowTable<>(infileS);
        HashMap<String, String> hm = new HashMap<>();
        for(int i = 0; i< t.getRowNumber(); i++){
            hm.put(t.getCellAsString(i, 0), t.getCellAsString(i, 3));
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(a);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.write("\tGroup");bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String taxon = PStringUtils.fastSplit(temp).get(0);
                if(hm.get(taxon).equals("Improved")) continue;
                if(hm.get(taxon).equals("unknown")) continue;
                bw.write(temp);bw.write("\t");bw.write(hm.get(taxon));bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void addGroupToMDS(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/005_taxa1210GroupInfo_mergeFeiGroupInfo_removeTIP498.2.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/006_MDS_addGroupInfo.txt";
        String a = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/MDS_PCs_Matrix_68300sites.txt";
        RowTable<String> t = new RowTable<>(infileS);
        HashMap<String, String> hm = new HashMap<>();
        HashMap<String, String> hm1 = new HashMap<>();
        
        int cntLandrace =0; int cntTeosinte =0; int cntTripsacum =0; int cntImproved =0; int cntunknown =0;
        for(int i =0; i< t.getRowNumber(); i++){
            hm.put(t.getCellAsString(i, 0), t.getCellAsString(i, 3));
            hm1.put(t.getCellAsString(i, 0), t.getCellAsString(i, 4));
            
            String germplasm = t.getCellAsString(i, 3);
            if(germplasm.equals("Landrace")) cntLandrace++;
            if(germplasm.equals("Teosinte")) cntTeosinte++;
            if(germplasm.equals("Tripsacum")) cntTripsacum++;
            if(germplasm.equals("Improved")) cntImproved++;
            if(germplasm.equals("unknown")) cntunknown++;
        }
        System.out.println(cntLandrace + "  Landrace\n" + cntTeosinte + "   Teosinte\n" + cntTripsacum + "  Tripsacum\n" + cntImproved + "  Improved\n" + cntunknown + "    unknown\n");
        /**
         * 25  Landrace
            21   Teosinte
            1  Tripsacum
            1111  Improved
            52    unknown
         */
        try{
            BufferedReader br = IOUtils.getTextReader(a);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.write("\tGroup");bw.write("\tFeiGroup");bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String taxon = PStringUtils.fastSplit(temp).get(0);
                bw.write(temp);bw.write("\t");bw.write(hm.get(taxon));bw.write("\t");bw.write(hm1.get(taxon));bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    public void mergeInfo(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/source/hmp321_taxaGroup_version3.txt";
        String finalfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/003_taxa1210GroupInfo_addGermplasm.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/004_taxa1210GroupInfo_mergeFeiGroupInfo.txt";
        RowTable<String> t = new RowTable<>(infileS);
        HashMap<String, String> hm = new HashMap<>();
        for(int i =0; i< t.getRowNumber(); i++){
            hm.put(t.getCellAsString(i, 0), t.getCellAsString(i, 1));
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(finalfileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.write("\tFeiGroup");bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String taxon = PStringUtils.fastSplit(temp).get(0);
                if (taxon.equals("TIP-498.2")){
                    bw.write(temp);bw.write("\tTeo");bw.newLine();
                }
                else{
                    bw.write(temp);bw.write("\t");bw.write(hm.get(taxon));bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void reversechecktaxaList(){
        String infile1S = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/source/Taxa1210List.txt";
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/003_taxa1210GroupInfo_addGermplasm.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(0);
        String[] taxa = l.toArray(new String[l.size()]);
        Arrays.sort(taxa);
        try{
            int cnt = 0;
            BufferedReader br = IOUtils.getTextReader(infile1S);
            String header = br.readLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String taxon = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(taxa, taxon);
                if (index >= 0) {
                    cnt++;
                }
                else{
                    System.out.println(temp);
                }
            }
            System.out.println(cnt + "\n");
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        /**
         * run:
            TIP-498.2	HapMap2 extra	Parviglumis	Teosinte
            1209

            W64a
            1209
         */
    }
    
    public void checktaxaList(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/source/Taxa1210List.txt";
        String infile1S = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/003_taxa1210GroupInfo_addGermplasm.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(0);
        String[] taxa = l.toArray(new String[l.size()]);
        Arrays.sort(taxa);
        try{
            int cnt = 0;
            BufferedReader br = IOUtils.getTextReader(infile1S);
            String header = br.readLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String taxon = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(taxa, taxon);
                if (index >= 0) {
                    cnt++;
                }
                else{
                    System.out.println(temp + "\n");
                    
                }
            }
            System.out.println(cnt + "\n");
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void checkCategoryandAddGermplasm(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/002_taxa1210GroupInfo_withoutSpacing.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/003_taxa1210GroupInfo_addGermplasm.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(2);
        Set<String> s = new HashSet<>(l);
        System.out.println(s + "\n");
        List<String> k = new ArrayList<>();
        for(int i= 0; i<t.getRowNumber(); i++){
            String category = t.getCellAsString(i, 2);
            String germp = null;
            if(category.equals("Landrace")){
                germp =  "Landrace";
            }
            if(category.equals("Mexicana") || category.equals("Parviglumis")){
                germp =  "Teosinte";
            }
            if(category.equals("Gamma Grass")){
                germp =  "Tripsacum";
            }
            if(category.equals("Improved")){
                germp =  "Improved";
                
            }
            if(category.equals("unknown")){
                germp =  "unknown";
            }
            k.add(germp);
        }
        t.insertColumn("Germplasm", 3, k);
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
    public void removeSpacing(){
        String a = "Improved";
        String b = a.trim();
        System.out.println(a);
        System.out.println(b);
        
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/001_taxa1210GroupInfo.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/002_taxa1210GroupInfo_withoutSpacing.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String Taxon = PStringUtils.fastSplit(temp).get(0);
                String Dataset = PStringUtils.fastSplit(temp).get(1);
                String Category = PStringUtils.fastSplit(temp).get(2);
                String category = Category.trim();
                StringBuilder sb = new StringBuilder();
                sb.append(Taxon).append("\t").append(Dataset).append("\t").append(category);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void taxaGroupSummary(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/source/hmp321taxaInfoSummary.txt";
        String CAUFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/002_SingleSummary/001_CAU.txt";
        String HapMap2extraFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/002_SingleSummary/001_HapMap2_extra.txt";
        String HapMap2FileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/002_SingleSummary/001_HapMap2.txt";
        
        String taxa9FileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/003_extractInfo/9taxaInHapMap2extra.txt";
        String taxa30FileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/003_extractInfo/30taxaHapMap2_HapMap2extra.txt";
        String taxa63FileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/003_extractInfo/63taxaInHapMap2.txt";
        
        String outFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/004_taxaGroup/001_taxa1210GroupInfo.txt";
        
        BufferedReader br = IOUtils.getTextReader(CAUFileS);
        BufferedWriter bw = IOUtils.getTextWriter(outFileS);
        
        RowTable<String> t1 = new RowTable<>(taxa9FileS);
        List<String> taxa9 = t1.getColumn(0);
        String[] taxa9a = taxa9.toArray(new String[taxa9.size()]);
        Arrays.sort(taxa9a);
        System.out.println(taxa9);
        BufferedReader br1 = IOUtils.getTextReader(HapMap2extraFileS);
        
        RowTable<String> t2 = new RowTable<>(taxa30FileS);
        List<String> taxa30 = t2.getColumn(0);
        String[] taxa30a = taxa30.toArray(new String[taxa30.size()]);
        Arrays.sort(taxa30a);
        
        RowTable<String> t3 = new RowTable<>(taxa63FileS);
        List<String> taxa63 = t3.getColumn(0);
        String[] taxa63a = taxa63.toArray(new String[taxa63.size()]);
        Arrays.sort(taxa63a);
        BufferedReader br2 = IOUtils.getTextReader(HapMap2FileS);
        
        try{
            String header = br.readLine();
            bw.write("Taxon\tDataSet\tCategory");bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String Taxon = PStringUtils.fastSplit(temp).get(1);
                String Dataset = PStringUtils.fastSplit(temp).get(2);
                String Category = PStringUtils.fastSplit(temp).get(4);
                StringBuilder sb = new StringBuilder();
                sb.append(Taxon).append("\t").append(Dataset).append("\t").append(Category);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            
            temp = br1.readLine();
            while((temp = br1.readLine()) != null){
                String Taxon = PStringUtils.fastSplit(temp).get(1);
                String Dataset = PStringUtils.fastSplit(temp).get(2);
                String Category = PStringUtils.fastSplit(temp).get(4);
                int index = Arrays.binarySearch(taxa9a, Taxon);
                if(index >= 0 ){
                    StringBuilder sb = new StringBuilder();
                    sb.append(Taxon).append("\t").append(Dataset).append("\t").append(Category);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                int index2 = Arrays.binarySearch(taxa30a, Taxon);
                if(index2 >= 0 ){
                    StringBuilder sb = new StringBuilder();
                    sb.append(Taxon).append("\t").append(Dataset).append("\t").append(Category);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            br1.close();
            
            temp = br2.readLine();
            while((temp = br2.readLine()) != null){
                String Taxon = PStringUtils.fastSplit(temp).get(1);
                String Dataset = PStringUtils.fastSplit(temp).get(2);
                String Category = PStringUtils.fastSplit(temp).get(4);
                int index = Arrays.binarySearch(taxa63a, Taxon);
                if(index >= 0 ){
                    StringBuilder sb = new StringBuilder();
                    sb.append(Taxon).append("\t").append(Dataset).append("\t").append(Category);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            br2.close();
            
            String taxa282setFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/003_extractInfo/282setTaxa.txt";
            RowTable<String> t4 = new RowTable<>(taxa282setFileS);
            List<String> taxa282set = t4.getColumn(0);
            String[] taxa282seta = taxa282set.toArray(new String[taxa282set.size()]);
            Arrays.sort(taxa282seta);        
            BufferedReader br3 = IOUtils.getTextReader(infileS);
            temp = br3.readLine();
            while((temp = br3.readLine()) != null){
                String Taxon = PStringUtils.fastSplit(temp).get(1);
                String Dataset = PStringUtils.fastSplit(temp).get(2);
                String Category = PStringUtils.fastSplit(temp).get(4);
                if(Dataset.equals("CIMMYT/BGI")){
                    StringBuilder sb = new StringBuilder();
                    sb.append(Taxon).append("\t").append(Dataset).append("\t").append(Category);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                if(Dataset.equals("282-4x"))continue;
                int index = Arrays.binarySearch(taxa282seta, Taxon);
                if(index >= 0 ){  
                    StringBuilder sb = new StringBuilder();
                    sb.append(Taxon).append("\t").append(Dataset).append("\t").append(Category);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                if(Dataset.equals("German")){
                    StringBuilder sb = new StringBuilder();
                    sb.append(Taxon).append("\t").append(Dataset).append("\t").append(Category);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            br3.close();
 
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void checkTaxaCategory(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/source/hmp321taxaInfoSummary.txt";
        RowTable<String> t = new RowTable<>(infileS);
        Set<String> DataSet = new HashSet<>();
        for(int i=0; i < t.getRowNumber(); i++){
            DataSet.add(t.getCellAsString(i, 2));
        }
        String[] datasets = DataSet.toArray(new String[DataSet.size()]);
        System.out.println( DataSet);
        for(int i=0; i< datasets.length; i++){
            try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            String header = br.readLine();
            String temp = null;
            Set<String> categorySet = new HashSet<>();
            while((temp = br.readLine()) != null){
                String Dataset = PStringUtils.fastSplit(temp).get(2);
                String category = PStringUtils.fastSplit(temp).get(4);
                if(Dataset.equals(datasets[i])){
                    categorySet.add(category);
                }
            }
            System.out.println(datasets[i] + categorySet);
            br.close();
            
            }
            catch(Exception e){
                //System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }
        }
        
        /**
        [CAU, HapMap2 extra, HapMap2, 282-2x, German, 282-4x, CIMMYT/BGI]
        CAU[ Improved]
        HapMap2 extra[Parviglumis,  Improved, landrace]
        HapMap2[Parviglumis, Mexicana, Gamma Grass,  Improved, Landrace]
        282-2x[ Improved]
        German[Parviglumis,  Improved]
        282-4x[ Improved]
        CIMMYT/BGI[ Improved, Landrace, unknown]
        */
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            String header = br.readLine();
            String temp = null;
            Set<String> categorySet = new HashSet<>();
            int cntUnknow = 0;
            int cntLandrace = 0;
            int cntImproved = 0;
            while((temp = br.readLine()) != null){
                String Dataset = PStringUtils.fastSplit(temp).get(2);
                String category = PStringUtils.fastSplit(temp).get(4);
                if(Dataset.equals("CIMMYT/BGI")){
                    if(category.equals(" Improved")) cntImproved++;
                    if(category.equals("Landrace")) cntLandrace++;
                    if(category.equals("unknown")) cntUnknow++;
                }
            }
            br.close();
            System.out.println( "#################" );
            System.out.println(String.valueOf(cntImproved) + "  Improved lines;" );
            System.out.println(String.valueOf(cntLandrace) + "  Landrace;" );
            System.out.println(String.valueOf(cntUnknow) + "  unknown;" );
            System.out.println( "#################" + "\n");
        }
        catch(Exception e){
            //System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
        }
        
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            String header = br.readLine();
            String temp = null;
            Set<String> taxaSet = new HashSet<>();
            //Boolean iftrue = false;
            int cntOverlap =0;
            while((temp = br.readLine()) != null){
                String Dataset = PStringUtils.fastSplit(temp).get(2);
                String category = PStringUtils.fastSplit(temp).get(4);
                String taxa = PStringUtils.fastSplit(temp).get(1);
                if(Dataset.equals("CAU")) taxaSet.add(taxa);
                
                if(Dataset.equals("HapMap2")) {
                    
                    if(!(taxaSet.add(taxa))){
                        cntOverlap++;
                        System.out.println(temp);
                    }
                    //taxaSet.add(taxa);
                    
                }
                
                //if(Dataset.equals("CAU")) taxaSet.add(taxa);
            }
            String[] taxa = taxaSet.toArray(new String[taxaSet.size()]);
            System.out.println(cntOverlap + "   HapMap2 and HapMap2 extra Overlap" );
            System.out.println( "#################" );
            System.out.println(taxaSet.size());
            
            br.close();
        }
        catch(Exception e){
            //System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static class B73Distance {

        public B73Distance() {
            //this.mkDistanceToB73();
            //this.mergeRecDelHighAndDsitanceToB73();
            this.filterGroup();
            
        }
        
        public void filterGroup(){
            String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/005_loadVSB73Distance/003_reccesiveDeleterious_merge.txt";
            String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/005_loadVSB73Distance/004_reccesiveDeleterious_merge_notraitMixedUnknow.txt";
            try{
                BufferedReader br = IOUtils.getTextReader(infileS);
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String header = br.readLine();
                bw.write(header);bw.newLine();
                String temp = null;
                while((temp = br.readLine()) != null){
                    String group = PStringUtils.fastSplit(temp).get(3);
                    if (group.equals("popcorn")) continue;
                    if (group.equals("mixed")) continue;
                    if (group.equals("sweet")) continue;
                    if (group.equals("unknown")) continue;
                    bw.write(temp);
                    bw.newLine();
                }
                br.close();bw.flush();bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        }
        
        private void mergeRecDelHighAndDsitanceToB73(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/reccesiveDeleterious_hmp321_highDepth.txt";
        String distanceFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/005_loadVSB73Distance/002_1210taxaDistanceToB73.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/005_loadVSB73Distance/003_reccesiveDeleterious_merge.txt";
        RowTable<String> t = new RowTable<>(distanceFileS);
        HashMap<String,String> taxaDistance = new HashMap<>();
        for(int i= 0; i< t.getRowNumber(); i++){
            taxaDistance.put(t.getCellAsString(i, 1),t.getCellAsString(i, 2));
        }
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.write("\tIBS_DistancetoB73");bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String key = PStringUtils.fastSplit(temp).get(0);
                bw.write(temp);bw.write("\t");bw.write(taxaDistance.get(key));
                bw.newLine();
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
        
        private void mkDistanceToB73(){
        String matrixFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/005_loadVSB73Distance/001_B73_Matrix_68300sites.txt";
        String taxaListFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/source/Taxa1210List.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/005_loadVSB73Distance/002_1210taxaDistanceToB73.txt";
        RowTable<String> t = new RowTable<>(taxaListFileS);
        List<String> taxaList = new ArrayList<>();
        taxaList = t.getColumn(0);
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        try{
            BufferedReader br = IOUtils.getTextReader(matrixFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxon1\tTaxon2\tIBS_Distance");
            bw.newLine();
            String temp = br.readLine();
            List<String> matrixList = PStringUtils.fastSplit(temp);
            for(int i = 0; i < taxaList.size(); i++){
                bw.write("B73\t");
                bw.write(taxaList.get(i));bw.write("\t");
                int j = i+1;
                bw.write(matrixList.get(j));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
        
        
        
    }
}
