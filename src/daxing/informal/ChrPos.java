package daxing.informal;

import daxing.common.StringTool;
import gnu.trove.list.array.TIntArrayList;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class ChrPos {
    String chr=null;
    TIntArrayList pos=null;

    private  void initialize(File filePath){
        pos=new TIntArrayList();
        String chromosome=null;
        TIntArrayList position =new TIntArrayList();
        try (BufferedReader br= IOUtils.getTextGzipReader(filePath.getAbsolutePath())){
            br.readLine();
            String temp;
            HashSet<String> chrSet=new HashSet<>();
            while ((temp=br.readLine())!=null){
                List<String> list= PStringUtils.fastSplit(temp);
                chrSet.add(list.get(0));
                pos.add(Integer.parseInt(list.get(1)));
            }
            if (chrSet.size()==1){
                chr=chrSet.iterator().next();
            }else {
                System.out.println(filePath+" has "+chrSet.size()+" chromosomes");
                System.exit(1);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    ChrPos(File filePath){
        this.initialize(filePath);
    }

    public void add(ChrPos c){
        if (this.chr.equals(c.chr)){
            this.pos.addAll(c.pos);
            this.pos.sort();
        }else {
            System.out.println("This two files are different chromosomes");
            System.exit(1);
        }
    }

    public void write(File outPath){
        try(BufferedWriter bw=IOUtils.getTextGzipWriter(outPath.getAbsolutePath())) {
            bw.write("Chr"+"\t"+"Pos");
            bw.newLine();
            for (int i = 0; i < this.pos.size() ; i++) {
                bw.write(this.chr+"\t"+String.valueOf(this.pos.get(i)));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public void write(File parentDir, String child ){
        this.write(new File(parentDir, child));
    }

    /**
     *
     * @param chrPosPathDir inputDir
     * @param outDir
     */
    public static void merge(String chrPosPathDir, String outDir){
        File[] files=IOUtils.listRecursiveFiles(new File(chrPosPathDir));
        Predicate<File> p= file -> {
            int chrNum=StringTool.getNumFromString(file.getName());
            if(chrNum==0 || chrNum==43 || chrNum==44) return true;
            return false;
        };
        List<File> filesList=Arrays.stream(files).filter(e->(!e.getName().contains(".DS_Store"))).filter(p.negate()).collect(Collectors.toList());
        List<ChrPos> chrPoslist= new ArrayList<>();
        for (int i = 0; i < filesList.size(); i++) {
            chrPoslist.add(new ChrPos(filesList.get(i)));
        }
        for (int i = 0; i < chrPoslist.size(); i=i+2) {
            chrPoslist.get(i).add(chrPoslist.get(i+1));
        }
        for (int i = 0; i < chrPoslist.size(); i=i+2) {
            chrPoslist.get(i).write(new File(outDir), StringTool.changeNumToChr(filesList.get(i).getName()));
        }
    }
}
