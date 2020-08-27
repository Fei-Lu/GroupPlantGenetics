package daxing;

import daxing.common.PlotTools;
import daxing.common.RowTableTool;
import daxing.common.VCF;
import pgl.infra.utils.IOUtils;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * 个体load来源于毕傲月(包括杂合基因型load)
 */
public class Temp20200818 {

    public static void main(String[] args){
        String inputDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/005_test/003_hexaploid_test_pos_forSlidingWindow";
        String outDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/005_test/004_hexaploid_test_pos_slidingWindow";
        String synNonDel="syn";
        splitChrAndSlidWindow(inputDir,outDir,"syn");
        splitChrAndSlidWindow(inputDir,outDir,"non");
        splitChrAndSlidWindow(inputDir,outDir,"del");
    }

    public static void splitChrAndSlidWindow(String inputDir, String outDir, String synNonDel){
        String[] synNonDelArray={"syn","non","del"};
        int columnIndex= Arrays.asList(synNonDelArray).indexOf(synNonDel)+2;
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        String[] outFileNames= files.stream().map(File::getName).
                map(s->s.replaceAll("txt.gz",synNonDel+".10M_slidWindow_1M_step.txt.gz")).toArray(String[]::new);
        String[] outDirNames=files.stream().map(File::getName).map(s->s.replaceAll(".txt.gz",".splitChr")).toArray(String[]::new);
        File[] subDir=new File[outDirNames.length];
        File[] subSlidWindowDir=new File[outDirNames.length];
        for (int i = 0; i < outDirNames.length; i++) {
            subDir[i]=new File(outDir, outDirNames[i]);
            subDir[i].mkdir();
            subSlidWindowDir[i]=new File(outDir, outDirNames[i]+"_slidWindow");
            subSlidWindowDir[i].mkdir();
        }
        VCF vcf;
        for (int i = 0; i < files.size(); i++) {
            vcf=new VCF(files.get(i));
            vcf.writeVcfToSplitedChr(subDir[i].getAbsolutePath());
        }
        for (int i = 0; i < subDir.length; i++) {
            PlotTools.slidingWindow(subDir[i].getAbsolutePath(),subSlidWindowDir[i].getAbsolutePath(),columnIndex);
        }
        for (int i = 0; i < subSlidWindowDir.length; i++) {
            RowTableTool.meregMultipleTable(subSlidWindowDir[i].getAbsolutePath(), "Chr",new File(outDir,
                    outFileNames[i]).getAbsolutePath());
        }
        for (int i = 0; i < subDir.length; i++) {
            for(File file: subDir[i].listFiles()){
                if (!file.isDirectory()){
                    file.delete();
                }
            }
            subDir[i].delete();
            for (File file: subSlidWindowDir[i].listFiles()){
                if (!file.isDirectory()){
                    file.delete();
                }
            }
            subSlidWindowDir[i].delete();
        }

    }

}
