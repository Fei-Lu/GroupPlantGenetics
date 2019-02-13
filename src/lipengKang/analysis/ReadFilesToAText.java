/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author kanglipeng
 */
public class ReadFilesToAText {

    String outputDirS = "/Users/kanglipeng/Desktop/result";

    String dir = "/Users/kanglipeng/Desktop/untitled folder";

    File libs[] = null;
    File aimlibs[] = null;
    File nameArray[] = null;
    File baseFileArray[] = null;

    public ReadFilesToAText() {
        this.ReadFilesAbsolutePath(dir, outputDirS);

    }

    public void ReadFilesAbsolutePath(String dir, String outputDirS) {

        File[] subDirs = new File(dir).listFiles();

        List<File> aimFileDirList = new ArrayList<>();
        List<File> baseFileList = new ArrayList<>();
        List<File> subbaseFileList = new ArrayList<>();
        Set<File> nameSet = new HashSet<>();

        for (int i = 0; i < subDirs.length; i++) {
            if (subDirs[i].isHidden()) {
                continue;
            }
            nameSet.add(subDirs[i]);
        }

        nameArray = nameSet.toArray(new File[nameSet.size()]);

        for (int i = 0; i < nameArray.length; i++) {
            if (!nameArray[i].isDirectory()) {
                continue;
            }
            baseFileList.add(nameArray[i]);
        }
        baseFileArray = baseFileList.toArray(new File[baseFileList.size()]);
        File[] subsubDirs = new File[baseFileList.size()];
        for (int i = 0; i < baseFileList.size(); i++) {
            subsubDirs = baseFileArray[i].listFiles();
            for (int j = 0; j < subsubDirs.length; j++) {
                subbaseFileList.add(subsubDirs[j]);
            }
        }
        try {
            Set<File> s = new HashSet<>(subbaseFileList);
            libs = s.toArray(new File[s.size()]);
            for (int i = 0; i < libs.length; i++) {
                if (libs[i].getName().endsWith(".txt")) {
                    aimFileDirList.add(libs[i]);
                }
            }

            Set<File> m = new HashSet(aimFileDirList);
            aimlibs = m.toArray(new File[m.size()]);

            String line = null;

            BufferedReader[] brs = new BufferedReader[aimlibs.length];
            String outfileS = new File(outputDirS, "Allmd5.txt").getAbsolutePath();
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfileS)));

            for (int i = 0; i < brs.length; i++) {
                String infileS = aimlibs[i].getAbsolutePath();
                brs[i] = utils.IOUtils.getTextReader(infileS);

                while ((line = brs[i].readLine()) != null) {
                    bw.write(line);
                }
                bw.newLine();
                bw.flush();
            }
            for (int i = 0; i < brs.length; i++) {

                brs[i].close();
            }
            bw.close();
            File workingDir = new File(this.outputDirS);
            workingDir.mkdir();
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Something wrong while reading ");

        }
    }

    public static void main(String[] args) {

        new ReadFilesToAText();
    }

}
