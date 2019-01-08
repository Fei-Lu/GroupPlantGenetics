package daxing;

import daxing.md5.PCA;
import format.table.RowTable;
import utils.IOUtils;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;

/**
 *
 * @author xudaxing
 */
public class Go {
    private String workingDirS;
    private String sampleIDFastqMapFileS;
    private String wellBarcodeFileS;
    private String[] subDirS = {"allFlowcellLaneIndexFastq","barcodeFile","barcodeFastq"};

    Go(String parameterFileS){
        this.initializeParameter(parameterFileS);
        this.writeFlowcellLaneIndexFastq();
        //this.addBarcodeAndgetLibraries();
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
        this.workingDirS = paList.get(0);
        this.sampleIDFastqMapFileS=paList.get(1);
        this.wellBarcodeFileS=paList.get(2);
        File workingDir = new File(this.workingDirS);
        workingDir.mkdir();
        for (String subDir : this.subDirS) {
            File f = new File(this.workingDirS, subDir);
            f.mkdir();
        }
        System.out.println("Pipeline parameters initialized");
    }

    private void writeFlowcellLaneIndexFastq(){
        Fastqs fqs = new Fastqs(this.sampleIDFastqMapFileS);
        String allFlowcellLaneIndexFastq=new File(this.workingDirS,this.subDirS[0]).getAbsolutePath();
        fqs.writeFlowcellLaneIndexFastq(allFlowcellLaneIndexFastq);
    }

    private void addBarcodeAndgetLibraries(){
        Fastqs fqs=new Fastqs(new File(this.workingDirS, this.subDirS[0]));
        MergeFlowcellLaneIndexFastq mf=new MergeFlowcellLaneIndexFastq(fqs);
        String barcodeFileDirS=new File(this.workingDirS,this.subDirS[1]).getAbsolutePath();
        String bacodeFile=new File(barcodeFileDirS, "BarcodeFile.txt").getAbsolutePath();
        RowTable<String> t=new RowTable<>(this.wellBarcodeFileS);
        String librariesOutputDirS=new File(this.workingDirS, this.subDirS[2]).getAbsolutePath();
        String librariesFastqMap=new File(barcodeFileDirS, "LibrariesFastqMap.txt").getAbsolutePath();
        mf.addBarcodeAndgetLibraries(t, bacodeFile,librariesOutputDirS,librariesFastqMap);
    }

    public static void main(String[] args) {
        //new Go("/Users/xudaxing/Desktop/AM/parameter.txt");
        //MD5.getMD5FromDir((new File("/Users/xudaxing/Desktop/out")),new File("/Users/xudaxing/Desktop/md5.txt"));
        //md5.checkMD5ForDir(new File("/Users/xudaxing/Desktop/md5.txt"));
        PCA.extractRandomRowFromFile("/Users/xudaxing/Desktop/mkTagsBySampleLog.txt",
                "/Users/xudaxing/Desktop/mkTagsBySampleLog2.txt",1000, true);
    }
}

