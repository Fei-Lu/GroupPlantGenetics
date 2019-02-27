package daxing;

import utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class PCA {
    public static void extractRandomRowFromFile(String inputFile, String outFile, Integer numberOfRowsForExtract,
                                                boolean head){
        try{
            long numberOfRows =Files.lines(Paths.get(inputFile)).count();
            List<String> l= new ArrayList<>((int)numberOfRowsForExtract);
            BufferedReader br= IOUtils.getTextReader(inputFile);
            BufferedWriter bw =IOUtils.getTextWriter(outFile);
            String line;
            while ((line=br.readLine())!=null){
                if (head) {br.readLine();}
                line=br.readLine();
                l.add(line);
            }
            br.close();
            Collections.shuffle(l, new Random());
            for(int i=0;i<numberOfRowsForExtract;i++){
                bw.write(l.get(i));
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
}
