package daxing.ancestralAllele;

import daxing.common.ChrConvertionRule;
import utils.IOUtils;

import java.io.BufferedReader;
import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;

public class Start {
    private String workingDir;
    private String chrConvertionRule;
    private String chrallvcfFile;
    private String outgroup1InputDir;
    private int indexOfWheatInOutGroup1;
    private String outgroup2InputDir;
    private int indexOfWheatInOutGroup2;
    private String[] subDir={"refOutgroupAllele", "merge", "sorted", "ancestralAllele"};

    Start(String parameterFileS){
        this.initializeParameter(parameterFileS);
        this.getOutgroupAllele();
//        this.sort();
//        this.merge();
//        this.getAncestrallAllele();
    }

    private void initializeParameter(String parameterFileS){
        ArrayList<String> paList = new ArrayList();
        try {
            boolean check = false;
            BufferedReader br = IOUtils.getTextReader(parameterFileS);
            if (!br.readLine().equals("Author: Daxing Xu")) check = true;
            if (!br.readLine().equals("Email: dxxu@genetics.ac.cn; xiaoxiaoms218@gmail.com")) check = true;
            if (!br.readLine().equals("Homepage: http://plantgeneticslab.weebly.com/")) check = true;
            if (check) {
                System.out.println("Please keep the author information, or the program quits.");
            }
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("!Parameter")) {
                    paList.add(br.readLine());
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.workingDir = paList.get(0);
        this.chrConvertionRule=paList.get(1);
        this.chrallvcfFile=paList.get(2);
        this.outgroup1InputDir=paList.get(3);
        this.indexOfWheatInOutGroup1=Integer.parseInt(paList.get(4));
        this.outgroup2InputDir=paList.get(5);
        this.indexOfWheatInOutGroup2=Integer.parseInt(paList.get(6));

        File workingDir = new File(this.workingDir);
        workingDir.mkdir();
        File f;
        for (int i = 0; i < subDir.length; i++) {
            f=new File(this.workingDir, subDir[i]);
            f.mkdir();
        }
        System.out.println("Pipeline parameters initialized");
    }

    private void getOutgroupAllele(){
        ChrConvertionRule chrConvertionRule=new ChrConvertionRule(Paths.get(this.chrConvertionRule));
        MAF maf1=new MAF(this.indexOfWheatInOutGroup1, chrConvertionRule, Paths.get(this.outgroup1InputDir));
        MAF maf2=new MAF(this.indexOfWheatInOutGroup2, chrConvertionRule, Paths.get(this.outgroup2InputDir));
        maf1.getAllele(new File(this.workingDir, this.subDir[0]).getAbsolutePath());
        maf2.getAllele(new File(this.workingDir, this.subDir[0]).getAbsolutePath());
    }

    private void merge(){
        File files= new File(this.workingDir, this.subDir[2]);
        MAF.merge(files.getAbsolutePath(), new File(this.workingDir, this.subDir[1]).getAbsolutePath());
    }

    private void sort(){
        String inDir=new File(this.workingDir, this.subDir[0]).getAbsolutePath();
        String outDir=new File(this.workingDir, this.subDir[2]).getAbsolutePath();
        MAF.sort(inDir, outDir);
    }

    private void getAncestrallAllele(){
        String input=new File(this.workingDir, this.subDir[1]).getAbsolutePath();
        String out=new File(this.workingDir, this.subDir[3]).getAbsolutePath();
        MAF.getAncestralAlleleParallel(input, out);
    }

    public static void main(String[] args) {
//        new Start(args[0]);
    }


}
