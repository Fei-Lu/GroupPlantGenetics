/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Treetest {

    public Treetest() {
        //this.colbyType();
        //this.colbyContinent();
        this.binarybyType();
        
        
    }
    
    
    
    /**
     * 根据育种状态信息进行 binary标记。
     * 3种标签：Traditional cultivar/Landrace  Landrace   Advanced/improved cultivar   Cultivar  Breeding/Research Material 
     * 3种颜色：#FFD700   #87CEFA  #FF6A6A
     * 3种标记：HV  2   3
     * 状态：hm1 为 landrace的状态标记；如果是 Traditional cultivar/Landrace，为1，如果是 Landrace，为0，如果是其他，为-1；
     *      hm2 为 cultivar的状态标记；如果是 Cultivar，为1，如果是Advanced/improved cultivar，为0；如果是其他，为-1；
     *      hm3 为 breeding material的状态标记；如果是 Breeding/Research Material，为1，如果不是，为-1.
     */
    
    public void binarybyType(){
        String[] type = {"Traditional cultivar/Landrace","Landrace","Advanced/improved cultivar","Cultivar","Breeding/Research Material","NA","Other"};
        HashMap<String,String> hm1 = new HashMap<>();
        HashMap<String,String> hm2 = new HashMap<>();
        HashMap<String,String> hm3 = new HashMap<>();
        //System.out.println(type);
        for(int i=0;i<type.length;i++){
            if(type[i].equals("Traditional cultivar/Landrace")){
                hm1.put(type[i], "1"); hm2.put(type[i], "-1"); hm3.put(type[i], "-1");
            }
            if(type[i].equals("Landrace")){
                hm1.put(type[i], "0"); hm2.put(type[i], "-1"); hm3.put(type[i], "-1");
            }
            if(type[i].equals("Cultivar")){
                hm1.put(type[i], "-1"); hm2.put(type[i], "1"); hm3.put(type[i], "-1");
            }
            if(type[i].equals("Advanced/improved cultivar")){
                hm1.put(type[i], "-1"); hm2.put(type[i], "0"); hm3.put(type[i], "-1");
            }
            if(type[i].equals("Breeding/Research Material")){
                hm1.put(type[i], "-1"); hm2.put(type[i], "-1"); hm3.put(type[i], "1");
            }
            if(type[i].equals("Other")){
                hm1.put(type[i], "-1"); hm2.put(type[i], "-1"); hm3.put(type[i], "-1");
            }
            if(type[i].equals("NA")){
                hm1.put(type[i], "-1"); hm2.put(type[i], "-1"); hm3.put(type[i], "-1");
            }
        }
        
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/source/002_merge/All373wheat_ABD_CountryBreedingStatus.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/002_addcolor/002_binarybyType.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            bw.write("DATASET_BINARY\n");
            bw.write("SEPARATOR TAB\n");
            bw.write("DATASET_LABEL\tbinary_data\n");
            bw.write("COLOR\t#ff0000\n\n");
            bw.write("FIELD_LABELS\tLandrace\tCultivar\tBreeding/Research Material\n");
            bw.write("FIELD_COLORS\t#FFD700\t#87CEFA\t#FF6A6A\n");
            bw.write("FIELD_SHAPES\t1\t2\t3\n\n");
            bw.write("DATA\n");
            while((temp = br.readLine()) != null){
                l =PStringUtils.fastSplit(temp);
                String DatabaseID = l.get(2);
                String status = l.get(13);
                bw.write(DatabaseID + "\t" + hm1.get(status) + "\t"+ hm2.get(status) +"\t"+ hm3.get(status) + "\n");
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    
    
    /**
     * 根据亚洲部洲和其他大洲信息，分组为10组进行着色显示
     */
    
    public void colbyContinent(){
        String[] group = {"Oceania","Africa","North America","South America","Europe","Central Asia","South Asia","Western Asia","East Asia","NA"};
        //String[] col = {"#F1E1FF","#F4D03F","#F1948A","#5DADE2","#ABEBC6","#239B56","#CD6155","#FF6347","#7B241C","#EBEDEF"};
        String[] col = {"#5DADE2","#7B241C","#F1E1FF","#F4D03F","#FF9900","#006600","#389038","#82C782","#CCFFCC","#EBEDEF"};
        HashMap<String,String> hm = new HashMap<>();
        for(int i=0;i<group.length;i++){
            hm.put(group[i], col[i]);
        }
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/source/002_merge/All373wheat_ABD_CountryBreedingStatus.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/002_addcolor/";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter[] bw = new BufferedWriter[2];
            bw[0] = IOUtils.getTextWriter(new File(outfileDirS,"001_addcol_branch_byContient.txt").getAbsolutePath());
            bw[1] = IOUtils.getTextWriter(new File(outfileDirS,"001_addcol_range_byContient.txt").getAbsolutePath());
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            bw[0].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            bw[1].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            while((temp = br.readLine()) != null){
                l =PStringUtils.fastSplit(temp);
                String DatabaseID = l.get(2);
                String partContinent = l.get(14);
                bw[0].write( DatabaseID + "\tbranch\t" + hm.get(partContinent) + "\tnormal");
                bw[1].write( DatabaseID + "\trange\t" + hm.get(partContinent) + "\tnormal");
                bw[0].newLine(); bw[1].newLine();
            }
            br.close();bw[0].flush();bw[1].flush();bw[0].close();bw[1].close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 根据育种状态，将373样品分为7类，并赋予其一定颜色：
     * Traditional cultivar/Landrace      PaleGreen3	124 205 124	#7CCD7C
     * Advanced/improved cultivar     	  LightGoldenrod1	255 236 139	#FFEC8B
     * Breeding/Research Material   DarkGreen	0 100 0	#006400
     * Landrace   Green3	0 205 0	#00CD00
     * Cultivar     Goldenrod1	255 193 37	#FFC125
     * Other LightBlue	173 216 230	#ADD8E6    
     * NA Black	0 0 0	#000000
     */
    
    public void colbyType(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/source/002_merge/All373wheat_ABD_CountryBreedingStatus.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/002_addcolor/";
        BufferedWriter[] bw = new BufferedWriter[2];
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            bw[0] = IOUtils.getTextWriter(new File(outfileDirS,"001_addcol_branch.txt").getAbsolutePath());
            bw[1] = IOUtils.getTextWriter(new File(outfileDirS,"001_addcol_range.txt").getAbsolutePath());
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while((temp = br.readLine()) != null){
                l =PStringUtils.fastSplit(temp);
                String DatabaseID = l.get(2);
                String type = l.get(13);
                if(type.equals("Traditional cultivar/Landrace")){
                    bw[0].write( DatabaseID + " branch #7CCD7C normal");
                    bw[1].write(DatabaseID + " range #7CCD7C Traditional cultivar/Landrace");
                    bw[0].newLine(); bw[1].newLine();
                }
                if(type.equals("Advanced/improved cultivar")){
                    bw[0].write( DatabaseID + " branch #FFEC8B normal");
                    bw[1].write(DatabaseID + " range #FFEC8B Advanced/improved cultivar");
                    bw[0].newLine(); bw[1].newLine();
                }
                if(type.equals("Breeding/Research Material")){
                    bw[0].write( DatabaseID + " branch #006400 normal");
                    bw[1].write(DatabaseID + " range #006400 Breeding/Research Material");
                    bw[0].newLine(); bw[1].newLine();
                }
                if(type.equals("Landrace")){
                    bw[0].write( DatabaseID + " branch #00CD00 normal");
                    bw[1].write(DatabaseID + " range #00CD00 Landrace");
                    bw[0].newLine(); bw[1].newLine();
                }
                if(type.equals("Cultivar")){
                    bw[0].write( DatabaseID + " branch #FFC125 normal");
                    bw[1].write(DatabaseID + " range #FFC125 Cultivar");
                    bw[0].newLine(); bw[1].newLine();
                }
                if(type.equals("Other")){
                    bw[0].write( DatabaseID + " branch #ADD8E6 normal");
                    bw[1].write(DatabaseID + " range #ADD8E6 Other");
                    bw[0].newLine(); bw[1].newLine();
                }
                if(type.equals("NA")){
                    bw[0].write( DatabaseID + " branch #000000 normal");
                    bw[1].write(DatabaseID + " range #000000 NA");
                    bw[0].newLine(); bw[1].newLine();
                }
            }
            br.close();bw[0].flush();bw[1].flush();bw[0].close();bw[1].close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
}
