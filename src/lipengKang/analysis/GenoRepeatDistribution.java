/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

/**
 *
 * @author kanglipeng
 * //a simplified counter to count unmasked region in soft-masked genome
 * //a bug : when the last base of seq is uppercase, it will always lose the last unmasked seq-length, but the numbers of unmasked segments is right
 *///need to skip each chr header(need modified)
// another bug when the unmasked seq in twoline, the length will -2
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

public class GenoRepeatDistribution {

    public GenoRepeatDistribution(String inFile, String outFile) {

        this.callNoRepeatLength(inFile, outFile);

    }

    private void callNoRepeatLength(String inFile, String outFile) {
        BufferedReader br;
       BufferedWriter bw;
        if (inFile.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(inFile);
            bw = YaoIOUtils.getTextGzipWriter(outFile);
        } else {
            br = YaoIOUtils.getTextReader(inFile);
            bw = YaoIOUtils.getTextWriter(outFile);
        }

        String temp = null;
        int noRepeatLength = 0;
        int noRepeatNumber = 1;
        TIntArrayList temp2 = new TIntArrayList();
        try {
            while ((temp = br.readLine()) != null) {
                int lineLength = temp.length();
                for (int i = 0; i <= lineLength - 2; i++) {
                    char x = temp.charAt(i);
                    char y = temp.charAt(i + 1);
                    if (Character.isUpperCase(x)) {
                        if (Character.isUpperCase(y)) {
                            noRepeatLength++;

                        } else {

                            temp2.add(noRepeatLength + 1);
                            noRepeatLength = 0;

                        }
                    } else {
                        if (Character.isUpperCase(y)) {
                            
                            noRepeatNumber++;
                        }

                    }

                }
            }
            System.out.print("There are" + " " + noRepeatNumber + " " + "unmasked segments in wheat genome");
         int repeatCount[]= temp2.toArray();
          StringBuilder sb = new StringBuilder();
            for (int i = 0; i <= repeatCount.length-1; i++) {
                sb.append(repeatCount[i]+"\n");
              
            }
            bw.write(sb.toString());
               
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Error in Reading: " + inFile);
        }
    }

    public static void main(String[] args) {
        String inFile = "/data1/home/lipeng/database/wheat/smRef/D/wheatD.fa";
        String outFile = "/data1/home/lipeng/result/whrepeatditr2.txt";
        new GenoRepeatDistribution(inFile, outFile);
    }

}
