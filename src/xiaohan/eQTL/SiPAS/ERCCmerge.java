package xiaohan.eQTL.SiPAS;

import pgl.infra.table.RowTable;
import xiaohan.eQTL.simulation.simulationData;
import xiaohan.utils.IOUtils;

import java.io.BufferedWriter;
import java.text.DecimalFormat;

/**
 * @ author: yxh
 * @ created: 2021-07-21 : 3:22 PM
 */
public class ERCCmerge {

    public ERCCmerge(){

    }

    public void ERCCRoc() {
        String prefix = "Truseq";
        int sample = 3;
        BufferedWriter bw = IOUtils.getTextWriter("/Users/yxh/Documents/eQTL/SiPAS/ERCC/20210318/5M/" + prefix + "/" + prefix + "PR.txt");
        RowTable<String> rt = new RowTable<>("/Users/yxh/Documents/eQTL/SiPAS/ERCC/20210318/5M/" + prefix + "/" + prefix + "Roc.txt");
        try {
            bw.write("Gene\tMix1_C\tMix1_precision\tMix1_recall\tMix2_C\tMix2_precision\tMix2_recall\n");
            for (int i = 0; i < rt.getRowNumber(); i++) {
                double[] TP = new double[sample * 2];
                double[] FP = new double[sample * 2];
                double[] FN = new double[sample * 2];
                for (int j = 1; j <= 2; j++) {
                    for (int k = 1; k <= sample; k++) {
                        String index = j + "_" + k;
                        String observed = rt.getCell(i, rt.getColumnIndex("log2mix" + index));
                        String expected = rt.getCell(i, rt.getColumnIndex("predictmix" + index));
                        double[] TPFPFN = simulationData.getPrecisionandRecall(expected, observed);
                        int ix = (j - 1) * sample + k - 1;
                        TP[ix] = TPFPFN[0];
                        FP[ix] = TPFPFN[1];
                        FN[ix] = TPFPFN[2];
                    }
                }
                StringBuilder sb = new StringBuilder();
                sb.append(rt.getCell(i, rt.getColumnIndex("Gene")) + "\t");
                sb.append(rt.getCell(i, rt.getColumnIndex("Mix1_C")) + "\t");
                DecimalFormat decfor = new DecimalFormat("0.00000000");
                double precision1 = 0.00000000;
                double recall1 = 0.00000000;
                double precision2 = 0.00000000;
                double recall2 = 0.00000000;
                double TP1 = 0.00000000;
                double TP2 = 0.00000000;
                double FP1 = 0.00000000;
                double FP2 = 0.00000000;
                double FN1 = 0.00000000;
                double FN2 = 0.00000000;
                for (int j = 0; j < sample; j++) {
                    TP1 += TP[j];
                    FP1 += FP[j];
                    FN1 += FN[j];
                }
                for (int j = sample; j < sample * 2; j++) {
                    TP2 += TP[j];
                    FP2 += FP[j];
                    FN2 += FN[j];
                }
                precision1 = (double) TP1 / (TP1 + FP1);
                recall1 = (double) TP1 / (TP1 + FN1);
                precision2 = (double) TP2 / (TP2 + FP2);
                recall2 = (double) TP2 / (TP2 + FN2);
                sb.append(precision1 + "\t" + recall1 + "\t");
                sb.append(rt.getCell(i, rt.getColumnIndex("Mix2_C")) + "\t");
                sb.append(precision2 + "\t" + recall2);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args){
        new ERCCmerge();
    }
}
