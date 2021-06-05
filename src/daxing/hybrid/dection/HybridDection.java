package daxing.hybrid.dection;

import daxing.common.DateTime;
import daxing.common.IOTool;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

public class HybridDection {

    public static void getDepthHeterozygosity(String inputDir, String outDir, int windowSize, int stepSize){
        List<File> fileList= IOTool.getVisibleFileListInDir(inputDir);
        String[] subdirs={"001_depth","002_heterozygosity","003_IndividualIntervalHeterozygosity"};
        List<File> subdirFileList=new ArrayList<>();
        File file;
        for (String subdir : subdirs) {
            file = new File(outDir, subdir);
            subdirFileList.add(file);
            file.mkdir();
        }
        String[] outFileNameDepth= fileList.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", ".depth.txt.gz")).toArray(String[]::new);
        String[] outFileNameHeterozygosity= fileList.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", ".heterozygosity.txt.gz")).toArray(String[]::new);
        String[] outFileNameIndividualIntervalHeterozygosity= fileList.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", ".IndividualIntervalHeterozygosity.txt.gz")).toArray(String[]::new);
        System.out.println(DateTime.getDateTimeOfNow());
        IntStream.range(0, fileList.size()).parallel().forEach(e->getDepth(fileList.get(e), new File(subdirFileList.get(0), outFileNameDepth[e])));
        IntStream.range(0, fileList.size()).parallel().forEach(e->getHeterozygosity(fileList.get(e), new File(subdirFileList.get(1), outFileNameHeterozygosity[e])));
        IntStream.range(0, fileList.size()).parallel().forEach(e->getIndividualIntervalHeterozygosity(e+1, windowSize
                ,stepSize, fileList.get(e), new File(subdirFileList.get(2), outFileNameIndividualIntervalHeterozygosity[e])));
    }

    private static void getDepth(File inputFile, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            List<String> temp, tem, te;
            List<String> taxaList;
            String line;
            int depth;
            while ((line=br.readLine()).startsWith("##")) {}
            temp=PStringUtils.fastSplit(line);
            taxaList=temp.subList(9,temp.size());
            StringBuilder sb=new StringBuilder();
            sb.append("Chr").append("\t").append("Pos").append("\t");
            sb.append(String.join("\t", taxaList));
            bw.write(sb.toString());
            bw.newLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t");
                for (int i = 9; i < temp.size(); i++) {
                    if (temp.get(i).startsWith("./.")){
                        sb.append("NA").append("\t");
                        continue;
                    }
                    tem=PStringUtils.fastSplit(temp.get(i), ":");
                    te=PStringUtils.fastSplit(tem.get(1), ",");
                    depth=Integer.parseInt(te.get(0))+Integer.parseInt(te.get(1));
                    sb.append(depth).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void getHeterozygosity(File inputFile, File outFile){
        GenotypeGrid genotypeGrid=new GenotypeGrid(inputFile.getAbsolutePath(), GenoIOFormat.VCF_GZ);
        int taxonNum=genotypeGrid.getTaxaNumber();
        short chr= genotypeGrid.getChromosome(0);
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            sb.append("Taxon").append("\t").append("Chr").append("\t").append("HeterozygoteNumber").append("\t");
            sb.append("HomozygoteNumber");
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i <taxonNum; i++) {
                sb.setLength(0);
                sb.append(genotypeGrid.getTaxonName(i)).append("\t").append(chr).append("\t");
                sb.append(genotypeGrid.getHeterozygoteNumberByTaxon(i)).append("\t");
                sb.append(genotypeGrid.getHomozygoteNumberByTaxon(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void getIndividualIntervalHeterozygosity(int chrID, int windowSize, int stepSize, File inputFile,
                                                           File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile)) {
            String line;
            List<String> temp, taxaList;
            while ((line=br.readLine()).startsWith("##")) {}
            temp=PStringUtils.fastSplit(line);
            taxaList=temp.subList(9,temp.size());
            List<IndividualIntervalHeterozygosity.IndividualHomoHeteroPosition> individualIntervalHeterozygosityList=new ArrayList<>();
            IndividualIntervalHeterozygosity.IndividualHomoHeteroPosition individualHomoHeteroPosition;
            for (String s : taxaList) {
                individualHomoHeteroPosition = new IndividualIntervalHeterozygosity.IndividualHomoHeteroPosition(s);
                individualIntervalHeterozygosityList.add(individualHomoHeteroPosition);
            }
            int position;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                position= Integer.parseInt(temp.get(1));
                for (int i = 9; i < temp.size(); i++) {
                    if (temp.get(i).startsWith("./.")) continue;
                    if (temp.get(i).startsWith("0/1")){
                        individualIntervalHeterozygosityList.get(i-9).addHeteroPos(position);
                    }else {
                        individualIntervalHeterozygosityList.get(i-9).addHomoPos(position);
                    }
                }
            }
            IndividualIntervalHeterozygosity individualIntervalHeterozygosity =
                    new IndividualIntervalHeterozygosity(taxaList, RefV1Utils.getChrIDLength(chrID), windowSize, stepSize);
            individualIntervalHeterozygosity.addTaxonHomoHeteroCount(individualIntervalHeterozygosityList);
            individualIntervalHeterozygosity.write(outFile.getAbsolutePath(), chrID);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
