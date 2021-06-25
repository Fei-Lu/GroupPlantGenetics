package daxing.hybrid.dection;

import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class HybridIdentification {

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
        String[] chrs=RefV1Utils.getChromosomes();
        IntStream.range(0, fileList.size()).forEach(e->getDepth(fileList.get(e), new File(subdirFileList.get(0), outFileNameDepth[e])));
//        IntStream.range(0, fileList.size()).parallel().forEach(e->getHeterozygosity(fileList.get(e), new File(subdirFileList.get(1), outFileNameHeterozygosity[e])));
        IntStream.range(0, fileList.size()).forEach(e->getIndividualIntervalHeterozygosity(chrs[e],
                windowSize,stepSize, fileList.get(e), new File(subdirFileList.get(2), outFileNameIndividualIntervalHeterozygosity[e])));
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
            Set<String> naSet;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t");
                naSet=new HashSet<>();
                for (int i = 9; i < temp.size(); i++) {
                    if (temp.get(i).startsWith("./.")){
                        sb.append("NA").append("\t");
                        naSet.add("NA");
                        continue;
                    }
                    tem=PStringUtils.fastSplit(temp.get(i), ":");
                    te=PStringUtils.fastSplit(tem.get(1), ",");
                    depth=Integer.parseInt(te.get(0))+Integer.parseInt(te.get(1));
                    naSet.add(String.valueOf(depth));
                    sb.append(depth).append("\t");
                }
                if (naSet.size() < 2) continue;
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

    public static void getIndividualIntervalHeterozygosity(String chr, int windowSize, int stepSize, File inputFile,
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
                    new IndividualIntervalHeterozygosity(taxaList, RefV1Utils.getChromosomeLength(chr), windowSize, stepSize);
            individualIntervalHeterozygosity.addTaxonHomoHeteroCount(individualIntervalHeterozygosityList);
            individualIntervalHeterozygosity.write(outFile.getAbsolutePath(), chr);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void calculateIntervalIBS(String hbVCFDir, String parentsDir, String df_indiIntervalHeterozigosityDir
            , String outDir){
        List<File> hbFiles= IOTool.getVisibleDir(hbVCFDir);
        List<File> paFiles = IOTool.getVisibleDir(parentsDir);
        List<File> intervalFiles= IOTool.getVisibleFileListInDir(df_indiIntervalHeterozigosityDir);
        GenotypeGrid hb, pa, genotypeGrid;
        BufferedReader br;
        BufferedWriter bw;
        String line, header;
        List<String> temp, p1List, p2List;
        Set<String> p1Set, p2Set;
        p1Set= RowTableTool.getColumnSet(intervalFiles.get(0).getAbsolutePath(), 5);
        p2Set= RowTableTool.getColumnSet(intervalFiles.get(0).getAbsolutePath(), 6);
        p1List = new ArrayList<>(p1Set);
        p2List = new ArrayList<>(p2Set);
        p1List.addAll(p2List);
        Collections.sort(p1List);
        try {
            String refChr, taxon;
            int refPosStart, refPosEnd, posStart, posEnd, chrIDStart, chrIDEnd;
            int site1Index, site2Index, taxon1Index, taxon2Index;
            double ibs;
            TaxonIBS taxonIBS;
            List<TaxonIBS> taxonIBSList;
            TDoubleList ibsList;
            List<String> taxaList;
            Set<String> taxaSet1, taxaSet2;
            StringBuilder sb;
            String[] outNames= hbFiles.stream().map(File::getName).map(str->str.replaceAll(".vcf.gz",".txt")).toArray(String[]::new);
            for (int i = 0; i < hbFiles.size(); i++) {
                hb = new GenotypeGrid(hbFiles.get(i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
                pa = new GenotypeGrid(paFiles.get(i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
                genotypeGrid = GenotypeOperation.mergeGenotypesByTaxon(hb, pa);
                br = IOTool.getReader(intervalFiles.get(i));
                bw = IOTool.getWriter(new File(outDir, outNames[i]));
                header=br.readLine();
                sb= new StringBuilder();
                sb.append(header).append("\t").append("IBS_Mini1").append("\t").append("IBS_Mini2").append("\t").append(
                        "IfP1P2Mini");
                bw.write(sb.toString());
                bw.newLine();
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    taxon = temp.get(3);
                    refChr = temp.get(0);
                    refPosStart = Integer.parseInt(temp.get(1));
                    refPosEnd = Integer.parseInt(temp.get(2));
                    chrIDStart = RefV1Utils.getChrID(refChr, refPosStart);
                    chrIDEnd = RefV1Utils.getChrID(refChr, refPosEnd);
                    posStart = RefV1Utils.getPosOnChrID(refChr, refPosStart);
                    posEnd = RefV1Utils.getPosOnChrID(refChr, refPosEnd);
                    site1Index = genotypeGrid.getSiteIndex( (short) chrIDStart, posStart);
                    site2Index = genotypeGrid.getSiteIndex( (short) chrIDEnd, posEnd);
                    site1Index = site1Index < 0 ? -site1Index-1 : site1Index;
                    site2Index = site2Index < 0 ? -site2Index-1 : site2Index;
                    taxonIBSList= new ArrayList<>();
                    taxon2Index = genotypeGrid.getTaxonIndex(taxon);
                    for (int j = 0; j < p1List.size(); j++) {
                        taxon1Index = genotypeGrid.getTaxonIndex(p1List.get(j));
                        ibs = genotypeGrid.getIBSDistance(taxon1Index, taxon2Index, site1Index, site2Index);
                        if (Double.isNaN(ibs)){
                            taxonIBS = new TaxonIBS(Double.MAX_VALUE, p1List.get(j));
                        }else {
                            taxonIBS = new TaxonIBS(ibs, p1List.get(j));
                        }
                        taxonIBSList.add(taxonIBS);
                    }
                    Collections.sort(taxonIBSList);
                    taxaList  = new ArrayList<>();
                    ibsList = new TDoubleArrayList();
                    for (int j = 0; j < 2; j++) {
                        taxaList.add(taxonIBSList.get(j).taxon);
                        ibsList.add(taxonIBSList.get(j).ibs);
                    }
                    taxaSet1 = new HashSet<>(taxaList);
                    taxaSet2 = new HashSet<>();
                    taxaSet2.add(temp.get(5));
                    taxaSet2.add(temp.get(6));
                    sb.setLength(0);
                    sb.append(line).append("\t").append(ibsList.get(0)).append("\t");
                    sb.append(ibsList.get(1)).append("\t");
                    if (taxaSet1.equals(taxaSet2) && ibsList.get(0)!=0 && ibsList.get(1)!=0){
                        sb.append("Yes");
                    }else {
                        sb.append("NO");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println(hbFiles.get(i).getName() + " completed");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static class TaxonIBS implements Comparable<TaxonIBS>{

        double ibs;
        String taxon;

        TaxonIBS(double ibs, String taxon){
            this.ibs=ibs;
            this.taxon=taxon;
        }

        @Override
        public int compareTo(TaxonIBS o) {
            if (this.ibs == o.ibs){
                return 0;
            }else if (this.ibs < o.ibs){
                return -1;
            }else return 1;
        }

        public String getTaxon() {
            return taxon;
        }
    }

}
