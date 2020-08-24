package daxing.load.ancestralSite.complementary;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

/**
 * 该类用于存储每个个体在亚基因水平上syn non del在cds区域为derived的概率
 * exon vcf
 * 基因cds区域
 */
public class IndividualTaxonDerivedProbability {

    public static void getIndividualTaxonDerivedProbability(String inputDir, String taxa_InfoDBFile, String outFile){
        List<File> fileList= IOUtils.getVisibleFileListInDir(inputDir);
        BufferedReader br;
        Map<String,String> taxonGroupMap= RowTableTool.getMap(taxa_InfoDBFile,0, 1);
        try (BufferedWriter bw = IOTool.getTextWriter(outFile)) {
            bw.write("Taxon\tSub\tDerivedSynCount\tDerivedNonCount\tDerivedDelCount\tGroup\tSubspecies" +
                    "\tSynRatio\tNonRatio\tDelRatio\tRatio_NonVSyn\tRatio_DelVSyn\tRatio_NonVSyn_ByCount" +
                    "\tRatio_DelVSyn_ByCount");
            bw.newLine();
            for (int i = 0; i < fileList.size(); i++) {
                br=IOTool.getReader(fileList.get(i));

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
