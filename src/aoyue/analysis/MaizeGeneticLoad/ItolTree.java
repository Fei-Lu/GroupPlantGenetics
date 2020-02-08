/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class ItolTree {

    public ItolTree() {
        this.addcolor();
    }
    
    public void addcolor (){
        String infileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/000_group/geneticGroup.manual.txt";
        //String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/014_tree/000_addcolorbranch.txt";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/014_tree/000_addcolorRange.txt";
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while((temp = br.readLine()) != null){
                l =PStringUtils.fastSplit(temp);
                String taxa = l.get(0);
                String geneticIndex = l.get(3);
                if(geneticIndex.equals("0")){
                    //bw.write( taxa + " branch #0072B2 normal");
                    bw.write( taxa + " range #0072B2 Teosinte");
                    bw.newLine();
                }
                if(geneticIndex.equals("1")){
                    //bw.write( taxa + " branch #D55E00 normal");
                    bw.write( taxa + " range #D55E00 Tropcal-subtropical");
                    bw.newLine();
                }
                if(geneticIndex.equals("2")){
                    //bw.write( taxa + " branch #56B4E9 normal");
                    bw.write( taxa + " range #56B4E9 Non-stiff_stalk");
                    bw.newLine();
                }
                if(geneticIndex.equals("3")){
                    //bw.write( taxa + " branch #009E73 normal");
                    bw.write( taxa + " range #009E73 Stiff_stalk");
                    bw.newLine();
                }
                if(geneticIndex.equals("4")){
                    //bw.write( taxa + " branch #E69F00 normal");
                    bw.write( taxa + " range #E69F00 Mixed");
                    bw.newLine();
                }
                if(geneticIndex.equals("5")){
                    //bw.write( taxa + " branch #999999 normal");
                    bw.write( taxa + " range #999999 China_specific");
                    bw.newLine();
                }
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
}
