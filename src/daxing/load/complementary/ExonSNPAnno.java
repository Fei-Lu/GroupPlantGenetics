package daxing.load.complementary;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import daxing.load.ancestralSite.SNPAnnotation;
import pgl.infra.pos.ChrPos;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.File;
import java.util.*;


/**
 * this class is use for reference bias
 * reference is ancestral or derived
 * daf bins
 * in each daf bin, 计算reference is ancestral条件下, non-synonymous ratio
 *                      reference is derived条件下, synonymous ratio
 */
public class ExonSNPAnno {

    ChrSNPAnnoDB[] chrSNPAnnoDB; // by chrID

    public ExonSNPAnno(String exonSNPAnnoDir){
        chrSNPAnnoDB=new ChrSNPAnnoDB[42];
        List<File> fileList= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        for (int i = 0; i < fileList.size(); i++) {
            chrSNPAnnoDB[i]=new ChrSNPAnnoDB(fileList.get(i));
        }
    }

    public ChrSNPAnnoDB[] getChrSNPAnnoDB() {
        return chrSNPAnnoDB;
    }

    private static double[] initializeDafBin(){
        double[] dafBins={3.85,7.59,11.3,15.1,18.8,22.5,26.2,30,45,60,70,80,90,100};
        double[] res=new double[dafBins.length];
        for (int i = 0; i < dafBins.length; i++) {
            res[i]=dafBins[i]/100;
        }
        return res;
    }

    private static int getDafBinIndex(double[] dafBins, double daf){
        int hit= Arrays.binarySearch(dafBins, daf);
        return hit < 0 ? -hit-1 : hit;
    }

    public Table<Integer,Integer,Double> getSubCorrRatio(String sub){
        int[] chrID= RefV1Utils.getChrIDsOfSubgenome(sub);
        List<ChrSNPAnnoDB> subSNPAnnoDBList=new ArrayList<>();
        ChrSNPAnnoDB[] chrSNPAnnoDBS=this.getChrSNPAnnoDB();
        for (int k : chrID) {
            subSNPAnnoDBList.add(chrSNPAnnoDBS[k - 1]);
        }
        double[] dafBins=initializeDafBin();
        List<SNPAnnotation> snpAnnotationList;
        SNPAnnotation snpAnnotation;
        List<SNPAnnotation>[][] refIsAncDafBinSNPList=new List[2][];
        for (int i = 0; i < refIsAncDafBinSNPList.length; i++) {
            refIsAncDafBinSNPList[i]=new List[dafBins.length];
            for (int j = 0; j < refIsAncDafBinSNPList[i].length; j++) {
                refIsAncDafBinSNPList[i][j]=new ArrayList<>();
            }
        }
        for (ChrSNPAnnoDB chrSNPAnnoDB: subSNPAnnoDBList){
            snpAnnotationList=chrSNPAnnoDB.getSnpAnnotationList();
            for (SNPAnnotation annotation : snpAnnotationList) {
                snpAnnotation = annotation;
                if (!snpAnnotation.hasAncestral()) continue;
                if (snpAnnotation.getDaf() < 0) continue;
                int dafBinIndex = getDafBinIndex(dafBins, snpAnnotation.getDaf());
                if (annotation.isRefAlleleAncestral()) {
                    refIsAncDafBinSNPList[0][dafBinIndex].add(snpAnnotation);
                } else {
                    refIsAncDafBinSNPList[1][dafBinIndex].add(snpAnnotation);
                }
            }
        }
        List<SNPAnnotation> ancDafBinSnpList, derDafBinSnpList;
        Table<Integer,Integer,Double> chrPosCorrRatioTable= HashBasedTable.create();
        ChrPos chrPos;
        for (int i = 0; i < dafBins.length; i++) {
            ancDafBinSnpList=refIsAncDafBinSNPList[0][i];
            derDafBinSnpList=refIsAncDafBinSNPList[1][i];
            double ancNonCount=0, derNonCount=0;
            double ancNonRatio, derNonRatio;
            for (SNPAnnotation snpAnnotation1: ancDafBinSnpList){
                if (!snpAnnotation1.getVariant_type().equals("NONSYNONYMOUS")) continue;
                ancNonCount++;
            }
            for (SNPAnnotation snpAnnotation1: derDafBinSnpList){
                if (!snpAnnotation1.getVariant_type().equals("NONSYNONYMOUS")) continue;
                derNonCount++;
            }
            ancNonRatio=ancNonCount/ancDafBinSnpList.size();
            derNonRatio=derNonCount/derDafBinSnpList.size();
            double corrRatio=ancNonRatio/derNonRatio;
            for (SNPAnnotation snpAnnotation1: ancDafBinSnpList){
                chrPos=snpAnnotation1.getChrPos();
                chrPosCorrRatioTable.put((int)chrPos.getChromosome(), chrPos.getPosition(), corrRatio);
            }
            for (SNPAnnotation snpAnnotation1: derDafBinSnpList){
                chrPos=snpAnnotation1.getChrPos();
                chrPosCorrRatioTable.put((int)chrPos.getChromosome(), chrPos.getPosition(), 1.0);
            }
        }
        return chrPosCorrRatioTable;
    }

    public Table<Integer,Integer,Double> getCorrRatio(){
        Table<Integer,Integer,Double> chrPosDoubleMap=HashBasedTable.create();
        chrPosDoubleMap.putAll(this.getSubCorrRatio("A"));
        chrPosDoubleMap.putAll(this.getSubCorrRatio("B"));
        chrPosDoubleMap.putAll(this.getSubCorrRatio("D"));
        return chrPosDoubleMap;
    }

}
