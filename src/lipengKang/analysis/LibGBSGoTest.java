/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import lipengKang.analysis.LibraryInfoTest;
import lipengKang.analysis.TagParserTest;
import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import utils.IOUtils;

/**
 *
 * @author kanglipeng
 */
public class LibGBSGoTest {

    String workingDirS = null;
    String barcodeFileS = null;
    String libraryFastqMapFileS = null;
    String[] subDirS = {"tagsBySample", "tagsLibrary", "alignment", "rawGenotype", "filteredGenotype"};
    LibraryInfoTest li = null;
    String tagBySampleDirS = null;

    public LibGBSGoTest(String parameterFileS) {
        this.initializeParameter(parameterFileS);
    }

    public void initializeParameter(String parameterFileS) {
        ArrayList<String> paList = new ArrayList();
        try {
            boolean check = false;
            BufferedReader br = IOUtils.getTextReader(parameterFileS);
            if (!br.readLine().equals("Author: Fei Lu")) {
                check = true;
            }
            if (!br.readLine().equals("Email: flu@genetics.ac.cn; dr.lufei@gmail.com")) {
                check = true;
            }
            if (!br.readLine().equals("Homepage: http://plantgeneticslab.weebly.com/")) {
                check = true;
            }
            if (check) {
                System.out.println("Please keep the author information, or the program quits.");
            }
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("!Parameter")) {
                    paList.add(br.readLine());
                }
            }
            br.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
        this.workingDirS = paList.get(0);
        this.barcodeFileS = paList.get(1);
        this.libraryFastqMapFileS = paList.get(2);
        
        File workingDir = new File(this.workingDirS);
        workingDir.mkdir();
        for (int i = 0; i < this.subDirS.length; i++) {
            File f = new File(this.workingDirS, subDirS[i]);
            f.mkdir();
            if (i == 0) {
                File tagsLibraryFile = new File(this.workingDirS, subDirS[0]);
                tagBySampleDirS = tagsLibraryFile.getAbsolutePath();
                continue;
            }

        }

    }

    public static void main(String[] args) {
        LibGBSGoTest liGBSGoTest=new LibGBSGoTest(args[0]);
        LibraryInfoTest li =  new LibraryInfoTest(liGBSGoTest.barcodeFileS,  liGBSGoTest.libraryFastqMapFileS);
        new LibraryInfoTest(liGBSGoTest.barcodeFileS, liGBSGoTest.libraryFastqMapFileS);
        new TagParserTest(li, liGBSGoTest.tagBySampleDirS);
    }
}
