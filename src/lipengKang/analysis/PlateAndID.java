/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

/**
 *
 * @author kanglipeng
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Collections;
import java.util.List;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class PlateAndID {
    public PlateAndID(){
        this.plateAndID();
    }
    public void plateAndID(){
        String PlateinfileS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/实验记录/Re2extractDNA_Plate.txt";
        String OriginalInfileS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/实验记录/GBS_ID2.txt";
        String outfileS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/实验记录/Re2-2extractDNA_ID.txt";
        RowTable<String> t = new RowTable<>(PlateinfileS);
        List<String> plateList = t.getColumn(2);
        Collections.sort(plateList);
        try {
            BufferedReader br = IOUtils.getTextReader(OriginalInfileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("SerialNumber\tPlate\tID");
            bw.newLine();
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                if (Collections.binarySearch(plateList, tem[1]) < 0) continue;
                cnt++;
                bw.write(String.valueOf(cnt)+"\t"+temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
        
    
}
