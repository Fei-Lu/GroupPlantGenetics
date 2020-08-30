package daxing.load.ancestralSite.complementary;

import daxing.common.IOTool;
import daxing.common.NumberTool;
import daxing.common.PGF;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class LoadHeter {

    /**
     * 使用滑窗计算syn non del 在个体水平上的的杂合分布
     * @param pgfFile
     * @param countMergeFile
     * @param windowSize
     * @param stepSize
     * @param outFile
     */
    public static void calculateHeterozygous(String pgfFile, String countMergeFile, int windowSize,
                                             int stepSize, String outFile){
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByName();
        try (BufferedReader br = IOTool.getReader(countMergeFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            br.readLine();
            String line, geneName;
            List<String> temp;
            int geneIndex, chrID, startOnChrID, endOnChrID, start, end, midPos, chrIndex;
            double[] heterSynNonDel, totalSynNonDel;
            double synInWindowSum, nonInWindowSum, delInWindowSum;
            double synHeterInWindowSum, nonHeterInWindowSum, delHeterInWindowSum;
            String chr;
            SlidingWindowForLoadComplement[] totalSimpleWindowArray=new SlidingWindowForLoadComplement[21];
            SlidingWindowForLoadComplement[] heterSimpleWindowArray=new SlidingWindowForLoadComplement[21]; // 第一维度为染色体，第二维度为syn non del
            String[] chrs=RefV1Utils.getChromosomes();
            for (int i = 0; i < chrs.length; i++) {
                totalSimpleWindowArray[i]=new SlidingWindowForLoadComplement(RefV1Utils.getChromosomeLength(chrs[i]), windowSize, stepSize);
                heterSimpleWindowArray[i]=new SlidingWindowForLoadComplement(RefV1Utils.getChromosomeLength(chrs[i]), windowSize, stepSize);
            }
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                geneName=temp.get(0);
                totalSynNonDel=new double[3];
                heterSynNonDel=new double[3];
                totalSynNonDel[0]=Integer.parseInt(temp.get(1));
                totalSynNonDel[1]=Integer.parseInt(temp.get(4));
                totalSynNonDel[2]=Integer.parseInt(temp.get(7));
                heterSynNonDel[0]=Integer.parseInt(temp.get(3));
                heterSynNonDel[1]=Integer.parseInt(temp.get(6));
                heterSynNonDel[2]=Integer.parseInt(temp.get(9));
                geneIndex=pgf.getGeneIndex(geneName);
                chrID=pgf.getGene(geneIndex).getGeneRange().getRangeChromosome();
                startOnChrID = pgf.getGene(geneIndex).getGeneRange().getRangeStart();
                endOnChrID=pgf.getGene(geneIndex).getGeneRange().getRangeEnd();
                chr= RefV1Utils.getChromosome(chrID, startOnChrID);
                start=RefV1Utils.getPosOnChromosome(chrID, startOnChrID);
                end=RefV1Utils.getPosOnChromosome(chrID, endOnChrID);
                midPos=(start+end)/2;
                chrIndex= Arrays.binarySearch(chrs, chr);
                totalSimpleWindowArray[chrIndex].addValue(midPos, totalSynNonDel);
                heterSimpleWindowArray[chrIndex].addValue(midPos, heterSynNonDel);
            }
            String header="Chr\tStart\tEnd\tGeneNum\tSynHeterRatio\tNonHeterRatio\tDelHeterRatio";
            bw.write(header);
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < totalSimpleWindowArray.length; i++) {
                for (int j = 0; j < totalSimpleWindowArray[i].getWindowNum(); j++) {
                    synInWindowSum=totalSimpleWindowArray[i].getSynWindowValue(j).sum();
                    synHeterInWindowSum= heterSimpleWindowArray[i].getSynWindowValue(j).sum();
                    nonInWindowSum=totalSimpleWindowArray[i].getNonWindowValue(j).sum();
                    nonHeterInWindowSum=heterSimpleWindowArray[i].getNonWindowValue(j).sum();
                    delInWindowSum=totalSimpleWindowArray[i].getDelWindowValue(j).sum();
                    delHeterInWindowSum=heterSimpleWindowArray[i].getDelWindowValue(j).sum();
                    sb.setLength(0);
                    sb.append(chrs[i]).append("\t").append(totalSimpleWindowArray[i].getWindowStarts()[j]).append("\t");
                    sb.append(totalSimpleWindowArray[i].getWindowEnds()[j]).append("\t").append(totalSimpleWindowArray[i].getGeneNum(j)).append("\t");
                    if (Double.isNaN(synHeterInWindowSum/synInWindowSum)){
                        sb.append("NaN").append("\t");
                    }else {
                        sb.append(NumberTool.format(synHeterInWindowSum/synInWindowSum, 5)).append("\t");
                    }
                    if (Double.isNaN(nonHeterInWindowSum/nonInWindowSum)){
                        sb.append("NaN").append("\t");
                    }else {
                        sb.append(NumberTool.format(nonHeterInWindowSum/nonInWindowSum, 5)).append("\t");
                    }
                    if (Double.isNaN(delHeterInWindowSum/delInWindowSum)){
                        sb.append("NaN");
                    }else {
                        sb.append(NumberTool.format(delHeterInWindowSum/delInWindowSum, 5));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
