package daxing.common.grt;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import daxing.common.genotype.GenoGrid;
import daxing.common.utiles.IOTool;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class HybridIdentification {

    public static void calculateIBS(String genotypeDir, String nam2ParentsFile, String outDir){
        List<File> vcfFiles = IOTool.getFileListInDirEndsWith(genotypeDir, ".vcf.gz");
        Nam2Parensts nam2Parensts = Nam2Parensts.getInstance(nam2ParentsFile);
        String[] outNames = vcfFiles.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".hybrid.txt")).toArray(String[]::new);
        IntStream.range(0,vcfFiles.size()).forEach(e->{
            try (BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                bw.write("SampleID\tHybridID\tP1\tP2\tMinIBSP1\tMinIBSP2\tMinIBS1\tMinIBS2");
                bw.newLine();
                StringBuilder sb = new StringBuilder();
                GenoGrid genoGrid = new GenoGrid(vcfFiles.get(e).getAbsolutePath(), GenoIOFormat.VCF_GZ);
                String[] sampleID = nam2Parensts.getSampleID();
                List<Parent2IBS> parent1ToIBSList, parent2ToIBSList;
                String[] p1s = Arrays.stream(nam2Parensts.getP1s()).filter(l->!l.equals("NA")).toArray(String[]::new);
                String[] p2s = Arrays.stream(nam2Parensts.getP2s()).filter(l->!l.equals("NA")).toArray(String[]::new);
                for (int i = 0; i < sampleID.length; i++) {
                    int sampleIndex = genoGrid.getTaxonIndex(sampleID[i]);
                    parent1ToIBSList = new ArrayList<>();
                    parent2ToIBSList = new ArrayList<>();
                    for (String parent: p1s){
                        int parentIndex = genoGrid.getTaxonIndex(parent);
                        double ibs = genoGrid.getIBSDistance(sampleIndex, parentIndex);
                        parent1ToIBSList.add(new Parent2IBS(parent, ibs));
                    }
                    for (String parent: p2s){
                        int parentIndex = genoGrid.getTaxonIndex(parent);
                        if (sampleIndex==parentIndex) continue;
                        double ibs = genoGrid.getIBSDistance(sampleIndex, parentIndex);
                        parent2ToIBSList.add(new Parent2IBS(parent, ibs));
                    }
                    Collections.sort(parent1ToIBSList, Comparator.comparing(l->l.getIbs()));
                    Collections.sort(parent2ToIBSList, Comparator.comparing(l->l.getIbs()));
                    sb.setLength(0);
                    sb.append(sampleID[i]).append("\t");
                    sb.append(nam2Parensts.getNamID(i)).append("\t");
                    sb.append(nam2Parensts.getP1(i)).append("\t");
                    sb.append(nam2Parensts.getP2(i)).append("\t");
//                    sb.append(genoGrid.getChromosome(0)).append("\t");
                    sb.append(parent1ToIBSList.get(0).parents).append("\t");
                    sb.append(parent2ToIBSList.get(0).parents).append("\t");
                    sb.append(parent1ToIBSList.get(0).ibs).append("\t");
                    sb.append(parent2ToIBSList.get(0).ibs).append("\t");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static class Nam2Parensts implements Swapper, IntComparator {
        String[] sampleID;
        String[] namIDs;
        String[] p1s;
        String[] p2s;

        public Nam2Parensts(String[] sampleID, String[] namIDs, String[] p1s, String[] p2s){
            this.sampleID=sampleID;
            this.namIDs=namIDs;
            this.p1s=p1s;
            this.p2s=p2s;
            this.sortBySampleID();
        }

        public String[] getNamIDs() {
            return namIDs;
        }

        public String[] getP1s() {
            return p1s;
        }

        public String[] getP2s() {
            return p2s;
        }

        public String[] getSampleID() {
            return sampleID;
        }

        @Override
        public void swap(int a, int b) {
            String temp= sampleID[a];
            sampleID[a]=sampleID[b];
            sampleID[b]=temp;
            temp = namIDs[a];
            namIDs[a] = namIDs[b];
            namIDs[b]=temp;
            temp = p1s[a];
            p1s[a] = p1s[b];
            p1s[b] = temp;
            temp = p2s[a];
            p2s[a] = p2s[b];
            p2s[b]=temp;
        }

        @Override
        public int compare(int a, int b) {
            return sampleID[a].compareTo(sampleID[b]);
        }

        public void sortBySampleID(){
            GenericSorting.quickSort(0, sampleID.length, this, this);
        }

        public int getSampleIndex(String sampleID){
            String[] sampleIDArray = this.getSampleID();
            int sampleIndex = Arrays.binarySearch(sampleIDArray, sampleID);
            return sampleIndex;
        }

        public String getP1(int sampleIndex){
            return this.getP1s()[sampleIndex];
        }

        public String getP2(int sampleIndex){
            return this.getP2s()[sampleIndex];
        }

        public String getNamID(int sampleIndex){
            return this.namIDs[sampleIndex];
        }

        public static Nam2Parensts getInstance(String nam2ParentsFile){
            List<String> sampleIDs= new ArrayList<>();
            List<String> namIDs = new ArrayList<>();
            List<String> p1s = new ArrayList<>();
            List<String> p2s = new ArrayList<>();
            try (BufferedReader br = IOTool.getReader(nam2ParentsFile)) {
                String line;
                List<String> temp;
                br.readLine();
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    sampleIDs.add(temp.get(0));
                    namIDs.add(temp.get(1));
                    p1s.add(temp.get(2));
                    p2s.add(temp.get(3));
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            String[] sampleIDArray = new String[sampleIDs.size()];
            String[] namIDArray = new String[namIDs.size()];
            String[] p1sArray = new String[p1s.size()];
            String[] p2sArray = new String[p2s.size()];
            for (int i = 0; i < sampleIDArray.length; i++) {
                sampleIDArray[i] = sampleIDs.get(i);
                namIDArray[i]=namIDs.get(i);
                p1sArray[i]=p1s.get(i);
                p2sArray[i]=p2s.get(i);
            }
            return new Nam2Parensts(sampleIDArray, namIDArray, p1sArray, p2sArray);
        }


    }

    public static class Parent2IBS{

        String parents;
        double ibs;

        public Parent2IBS(String parents, double ibs){
            this.parents=parents;
            this.ibs=ibs;
        }

        public double getIbs() {
            return ibs;
        }

        public String getParents() {
            return parents;
        }
    }
}
