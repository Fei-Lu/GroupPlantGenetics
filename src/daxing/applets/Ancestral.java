package daxing.applets;

import daxing.common.*;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.lang.StringUtils;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.IntPredicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class Ancestral {

    public static void getAncestral(String ancestralDir, String outDir, int outgroupNum, double thresh, String chrConvertionRule){
        List<File> files= IOUtils.getVisibleFileListInDir(ancestralDir);
        File refChr=new File(outDir, "001_byRefChr");
        File chrID=new File(outDir, "002_byChrID");
        refChr.mkdir();
        chrID.mkdir();
        if (outgroupNum==2){
            IntStream.range(0, files.size()).forEach(e->getAncestral_twoOutgroup(files.get(e), refChr, thresh));
            split(refChr.getAbsolutePath(), chrID.getAbsolutePath());
        }else if (outgroupNum==3){
            IntStream.range(0, files.size()).forEach(e->getAncestral_threeOutgroup(files.get(e), refChr, thresh));
            split(refChr.getAbsolutePath(), chrID.getAbsolutePath());
        }
    }

    private static String getAncestralAllele_twoOutGroup(String line, double thresh){
        String[] acgtArray={"A","C","G","T"};
        String[] temp= StringUtils.split(line, "\t ");
        double p_major_ancestral= Double.parseDouble(temp[6]);
        int[] countOfACGT= Stream.of(temp[3].split(",")).mapToInt(Integer::parseInt).toArray();
        int majorAlleleIndex= ArrayTool.getMaxIndex(countOfACGT);
        double[] p_alt= Arrays.asList(temp).subList(7, 10).stream().mapToDouble(Double::parseDouble).toArray();
        IntPredicate majorAllelePredicate= e->e==majorAlleleIndex;
        double[] altAlleleCount= IntStream.range(0, p_alt.length).filter(majorAllelePredicate.negate())
                .boxed().map(e->p_alt[e]).mapToDouble(Double::doubleValue).toArray();
        int maxIndexOfAltAllele=ArrayTool.getMaxIndex(altAlleleCount);
        double p_alt_ancestral=altAlleleCount[maxIndexOfAltAllele];
        List<String> acgtRemain= Arrays.stream(acgtArray).collect(Collectors.toList());
        acgtRemain.remove(majorAlleleIndex);
        if (p_major_ancestral>thresh) {
            return acgtArray[majorAlleleIndex];
        } else if (p_alt_ancestral>thresh) {
            return acgtRemain.get(maxIndexOfAltAllele);
        }
        return "NA";
    }

    private static String getAncestralAllele_threeOutgroup(String line, double thresh){
        String[] acgtBase={"A","C","G","T"};
        String[] temp= StringUtils.split(line, "\t ");
        int[] focal_ACGT_Count=Stream.of(temp[3].split(",")).mapToInt(Integer::parseInt).toArray();
        int majorAlleleIndex=ArrayTool.getMaxIndex(focal_ACGT_Count);
        double p_major_ancestral=Double.parseDouble(temp[6+3]);
        if (p_major_ancestral>thresh) return acgtBase[majorAlleleIndex];
        List<Double> node1States= Arrays.stream(temp).skip(7+3).mapToDouble(Double::parseDouble).boxed()
                .collect(Collectors.toList());
        TDoubleArrayList node1States_ACGT=new TDoubleArrayList();
        for (int i = 0; i < acgtBase.length; i++) {
            node1States_ACGT.add(node1States.subList(4*i, 4*i+4).stream().mapToDouble(Double::doubleValue).sum());
        }
        node1States_ACGT.removeAt(majorAlleleIndex);
        int indexOfMaxAF_Node1=ArrayTool.getMaxIndex(node1States_ACGT.toArray());
        List<String> acgt_remanin_Base=CollectionTool.changeToList(acgtBase);
        acgt_remanin_Base.remove(majorAlleleIndex);
        if (node1States_ACGT.get(indexOfMaxAF_Node1)>thresh) return acgt_remanin_Base.get(indexOfMaxAF_Node1);
        return "NA";
    }

    private static void getAncestral_twoOutgroup(File inputFile, File outDir, double thresh){
        try (BufferedReader bufferedReader = IOTool.getReader(inputFile)) {
            String outFileName="chr"+inputFile.getName().substring(0,2)+".secer.hovul.exon.probs.ancestral.txt.gz";
            BufferedWriter bw=IOTool.getWriter(new File(outDir, outFileName));
            WheatLineage subgenome=WheatLineage.valueOf(inputFile.getName().substring(1,2));
            bw.write("chr\tpos\tancestralAllele\n");
            String line, chr, pos, ancestralAllele;
            StringBuilder sb;
            List<String> temp;
            while ((line= bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                chr=temp.get(0);
                pos=temp.get(2);
                ancestralAllele= getAncestralAllele_twoOutGroup(line, thresh);
                sb=new StringBuilder(32);
                sb.append(chr, 3, 4).append(subgenome).append("\t").append(pos).append("\t").append(ancestralAllele);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void getAncestral_threeOutgroup(File inputFile, File outDir, double thresh){
        try (BufferedReader bufferedReader = IOTool.getReader(inputFile)) {
            Set<String> chrs= RowTableTool.getColumnSet(inputFile.getAbsolutePath(), 0);
            chrs.remove("0");
            List<String> chrsList=new ArrayList<>(chrs);
            Collections.sort(chrsList);
            Map<String, BufferedWriter> chrBufferWriterMap=new HashMap<>();
            BufferedWriter bw;
            WheatLineage subgenome=WheatLineage.valueOf(inputFile.getName().substring(5,6));
            String outgroupName=inputFile.getName().replaceAll("traes.", "").replaceAll("\\.gz$", "");
            for (int i = 0; i < chrsList.size(); i++) {
                bw=IOTool.getWriter(new File(outDir, chrsList.get(i)+subgenome+outgroupName+".ancestral.txt"));
                bw.write("chr\tpos\tancestralAllele\n");
                chrBufferWriterMap.put(chrsList.get(i), bw);
            }
            String line, chr, pos, ancestralAllele;
            List<String> temp;
            StringBuilder sb=new StringBuilder(32);
            while ((line=bufferedReader.readLine()).startsWith("0")){}
            temp= PStringUtils.fastSplit(line);
            chr=temp.get(0);
            pos=temp.get(2);
            ancestralAllele= getAncestralAllele_threeOutgroup(line, thresh);
            sb.append(chr.substring(3,4)+subgenome).append("\t").append(pos).append("\t").append(ancestralAllele);
            bw=chrBufferWriterMap.get(chr);
            bw.write(sb.toString());
            bw.newLine();
            while ((line= bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                chr=temp.get(0);
                pos=temp.get(2);
                ancestralAllele= getAncestralAllele_threeOutgroup(line, thresh);
                sb=new StringBuilder(32);
                sb.append(chr.substring(3,4)+subgenome).append("\t").append(pos).append("\t").append(ancestralAllele);
                bw=chrBufferWriterMap.get(chr);
                bw.write(sb.toString());
                bw.newLine();
            }
            for (Map.Entry<String, BufferedWriter> entry: chrBufferWriterMap.entrySet()){
                entry.getValue().flush();
                entry.getValue().close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void split(String inputDir, String outDir){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        VCF vcf;
        for (int i = 0; i < files.size(); i++) {
            vcf=new VCF(files.get(i));
            vcf.changeToVCFChr();
            vcf.writeVcfToSplitedChrID(outDir);
        }
    }

    public static void getProbabilityDistributionMap_two(String inputDir, String outFile){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        BufferedReader bufferedReader;
        BufferedWriter bufferedWriter=IOTool.getWriter(outFile);
        String line, chr, pos;
        TDoubleArrayList node2_AF;
        String[] temp;
        StringBuilder sb=new StringBuilder(32);
        WheatLineage wheatLineage;
        int[] count_ACGT;
        try {
            bufferedWriter.write("chrpos\tp_major_ancestral\tp_alt_maxAF\tp_alt_maxAF_removedMajor\n");
            for (int i = 0; i < files.size(); i++) {
                bufferedReader=IOTool.getReader(files.get(i));
                wheatLineage=WheatLineage.valueOf(files.get(i).getName().substring(1,2));
                while ((line=bufferedReader.readLine())!=null){
                    temp=StringUtils.split(line, "\t ");
                    chr=temp[0].substring(0,4)+wheatLineage;
                    pos=temp[2];
                    node2_AF= new TDoubleArrayList(Arrays.stream(temp).skip(7).collect(Collectors.toList()).stream().mapToDouble(Double::parseDouble).toArray());
                    int maxIndexOfNode2AF=ArrayTool.getMaxIndex(node2_AF.toArray());
                    sb=new StringBuilder();
                    sb.append(chr).append("_").append(pos).append("\t").append(temp[6]).append("\t").append(node2_AF.get(maxIndexOfNode2AF));
                    count_ACGT= Arrays.stream(Arrays.stream(temp).skip(3).limit(1).collect(Collectors.toList()).get(0).split(","))
                            .mapToInt(Integer::parseInt).toArray();
                    int indexOfMajorAllele=ArrayTool.getMaxIndex(count_ACGT);
                    node2_AF.removeAt(indexOfMajorAllele);
                    maxIndexOfNode2AF=ArrayTool.getMaxIndex(node2_AF.toArray());
                    sb.append("\t").append(node2_AF.get(maxIndexOfNode2AF));
                    bufferedWriter.write(sb.toString());
                    bufferedWriter.newLine();
                }
                bufferedReader.close();
            }
            bufferedWriter.flush();
            bufferedWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void getProbabilityDistributionMap_three(String inputDir, String outFile){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        BufferedReader bufferedReader;
        BufferedWriter bufferedWriter=IOTool.getWriter(outFile);
        String line, chr, pos;
        TDoubleArrayList alt_AF;
        double[] node1_ACGT;
        String[] temp;
        StringBuilder sb=new StringBuilder(32);
        WheatLineage wheatLineage;
        try {
            bufferedWriter.write("chrpos\tp_major_ancestral\tp_alt_maxAF\n");
            for (int i = 0; i < files.size(); i++) {
                bufferedReader=IOTool.getReader(files.get(i));
                wheatLineage=WheatLineage.valueOf(files.get(i).getName().substring(5,6));
                while ((line=bufferedReader.readLine()).startsWith("0")){}
                temp=StringUtils.split(line, "\t ");
                chr=temp[0]+wheatLineage;
                pos=temp[2];
                node1_ACGT=new double[4];
                alt_AF=new TDoubleArrayList(Arrays.stream(temp).skip(10).collect(Collectors.toList()).stream().mapToDouble(Double::parseDouble).toArray());
                for (int j = 0; j < node1_ACGT.length; j++) {
                    node1_ACGT[j]= alt_AF.subList(4*j, 4*j+4).sum();
                }
                int maxIndexOfAltAF=ArrayTool.getMaxIndex(node1_ACGT);
                sb=new StringBuilder();
                sb.append(chr).append("_").append(pos).append("\t").append(temp[9]).append("\t").append(node1_ACGT[maxIndexOfAltAF]);
                bufferedWriter.write(sb.toString());
                bufferedWriter.newLine();
                while ((line=bufferedReader.readLine())!=null){
                    temp=StringUtils.split(line, "\t ");
                    chr=temp[0]+wheatLineage;
                    pos=temp[2];
                    node1_ACGT=new double[4];
                    alt_AF=new TDoubleArrayList(Arrays.stream(temp).skip(10).collect(Collectors.toList()).stream().mapToDouble(Double::parseDouble).toArray());
                    for (int j = 0; j < node1_ACGT.length; j++) {
                        node1_ACGT[j]= alt_AF.subList(4*j, 4*j+4).sum();
                    }
                    maxIndexOfAltAF=ArrayTool.getMaxIndex(node1_ACGT);
                    sb=new StringBuilder();
                    sb.append(chr).append("_").append(pos).append("\t").append(temp[9]).append("\t").append(node1_ACGT[maxIndexOfAltAF]);
                    bufferedWriter.write(sb.toString());
                    bufferedWriter.newLine();
                }
                bufferedReader.close();
            }
            bufferedWriter.flush();
            bufferedWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
