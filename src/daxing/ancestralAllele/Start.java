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
    private String[] subDir={"refOutgroupAllele0904", "merge0904", "ancestralAllele0904"};

    Start(String parameterFileS){
        this.initializeParameter(parameterFileS);
        this.getOutgroupAllele();
//        this.merge();
//        this.getAncestralAllele();

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
//        AllelesInfor allelesInfor=new AllelesInfor(Paths.get(this.chrallvcfFile));
        maf1.getAllele(new File(this.workingDir, this.subDir[0]).getAbsolutePath());
        maf2.getAllele(new File(this.workingDir, this.subDir[0]).getAbsolutePath());
    }

    private void merge(){
        File[] files= new File(this.workingDir, this.subDir[0]).listFiles();
//        String outfile=new File(this.workingDir, "ancestralAllele.txt").getAbsolutePath();
//        MAF.mergeTwoFiles(files[0].getAbsolutePath(), files[1].getAbsolutePath(), outfile);
        MAF.merge(this.workingDir, new File(this.workingDir, this.subDir[1]).getAbsolutePath());
    }

    public void getAncestralAllele(){
        File input=new File(this.workingDir, this.subDir[1]);
        File output=new File(this.workingDir, this.subDir[2]);
        MAF.getAncestralAlleleParallel(input.getAbsolutePath(), output.getAbsolutePath());
    }


    public static void main(String[] args) {
        new Start(args[0]);
//        MD5.getMD5FromDir("/Users/xudaxing/Desktop/work");
//        SeqByte seqByte=new SeqByte("TCTTCCCCTA");
//        System.out.println(seqByte.getSequence(0,3));
//        Script.getLastz(args[0], args[1], args[2], args[3]);
    }
}
