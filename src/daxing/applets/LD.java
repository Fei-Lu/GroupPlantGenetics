package daxing.applets;

import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.CombinatoricsUtils;
import utils.IOFileFormat;
import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.IntStream;


/**
 * plink --vcf chrAll.tauschii.vcf  --chr-set 21 --make-bed --out chrAll.tauschii
 * plink --bfile chrAll.tauschii --r2 --matrix
 * note: plink.ld outFile not sorted, and include inter-Chromosome R2
 * using plink.ld and bim File as input, getting result file for ggplot
 */
public class LD {

    List<TDoubleArrayList> triangle;
    List<SNP> header;

    public LD(File ldMatrixFile, File bimFile){
        List<TDoubleArrayList> triangle=new ArrayList<>();
        List<SNP> snpList=new ArrayList<>();
        int i=0;
        try (BufferedReader br1 = IOUtils.getTextReader(ldMatrixFile.getAbsolutePath());
             BufferedReader br2=IOUtils.getTextReader(bimFile.getAbsolutePath())) {
            String line;
            String[] temp;
            TDoubleArrayList t;
            SNP snp;
            while ((line=br1.readLine())!=null){
                temp= StringUtils.split(line, "\t");
                t=new TDoubleArrayList();
                for (i = 0; i < temp.length; i++) {
                    if (temp[i].equals("nan")){
                        t.add(-1d);
                    }else {
                        t.add(Double.parseDouble(temp[i]));
                    }

                }
                triangle.add(t);
            }
           while ((line=br2.readLine())!=null){
               temp= StringUtils.split(line, "\t");
               snp=new SNP(Integer.parseInt(temp[0]), Integer.parseInt(temp[3]));
               snpList.add(snp);
           }
           this.triangle=triangle;
           this.header=snpList;
        }catch (Exception e){
            System.out.println(i);
            e.printStackTrace();
        }
    }

    /**
     *
     * @param snpIndex1
     * @param snpIndex2
     * @return double or -1d if nan exist
     */
    public double getR2(int snpIndex1, int snpIndex2){
        if (snpIndex2>snpIndex1){
            return this.triangle.get(snpIndex2).get(snpIndex1);
        }
        return this.triangle.get(snpIndex1).get(snpIndex2);
    }

    public SNP getSNP(int snpIndex){
        return this.header.get(snpIndex);
    }

    public boolean isOnSameChromosome(int snpIndex1, int snpIndex2){
        SNP snp1=this.getSNP(snpIndex1);
        SNP snp2=this.getSNP(snpIndex2);
        return snp1.chrID==snp2.chrID ? true: false;
    }

    /**
     *
     * @param snpIndex1
     * @param snpIndex2
     * @return  positive value or -1 if not on same chromosome
     */
    public int caculateDistanceBetweenSNPs(int snpIndex1, int snpIndex2){
        if (!this.isOnSameChromosome(snpIndex1, snpIndex2)) return -1;
        SNP snp1=this.getSNP(snpIndex1);
        SNP snp2=this.getSNP(snpIndex2);
        return Math.abs(snp1.getPos()-snp2.getPos());
    }

    /**
     * 过滤不在同一染色体上的r2值, r2为-1(NaN)的值, distance大于distancesThresh的r2值
     * @param physicalDistanceR2File
     * @param distancesThresh
     */
    public void  writePhysicalDistanceR2ForLDDecay(String physicalDistanceR2File, int distancesThresh){
        try (BufferedWriter bw = IOUtils.getTextWriter(physicalDistanceR2File)) {
            bw.write("Distance\tr2\n");
//            bw.write("SNP1\tSNP2\tDistance\tr2\n");
            Iterator<int[]> iterator = CombinatoricsUtils.combinationsIterator(header.size(), 2);
            int[] combinationIndex;
            StringBuilder sb;
            int count=0;
            int cnt=0;
            int cnnt=0;
            while (iterator.hasNext()) {
                combinationIndex = iterator.next();
                if (!isOnSameChromosome(combinationIndex[0], combinationIndex[1])) {
                    count++;
                    continue;
                }
                if (this.getR2(combinationIndex[0], combinationIndex[1])<0) {
                    cnt++;
                    continue;
                }
                if (this.caculateDistanceBetweenSNPs(combinationIndex[0], combinationIndex[1])>distancesThresh){
                    cnnt++;
                    continue;
                }
                sb=new StringBuilder();
//                sb.append(this.getSNP(combinationIndex[0]).toString()).append("\t");
//                sb.append(this.getSNP(combinationIndex[1]).toString()).append("\t");
                sb.append(this.caculateDistanceBetweenSNPs(combinationIndex[0], combinationIndex[1])).append("\t");
                sb.append(this.getR2(combinationIndex[0], combinationIndex[1]));
                bw.write(sb.toString());
                bw.newLine();
            }
            System.out.println("On different chromosome: "+count);
            System.out.println("r2 < 0: "+cnt);
            System.out.println("great than "+distancesThresh+": "+cnnt);
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * 过滤r2为-1(NaN)的值, 且抽取rate
     * @param physicalDistanceR2File
     */
    public void  writePhysicalDistanceR2ForLDDecay(String physicalDistanceR2File, double rate){
        try (BufferedWriter bw = IOUtils.getTextWriter(physicalDistanceR2File)) {
            bw.write("Distance\tr2\n");
            List<TDoubleArrayList> triangle=this.triangle;
            double cnt=0;
            int total=0;
            StringBuilder sb;
            for (int i = 0; i < triangle.size(); i++) {
                for (int j = 0; j < triangle.get(i).size(); j++) {
                    if (i==j) continue;
                    total++;
                    if (this.getR2(i, j)<0){
                        cnt++;
                        continue;
                    }
                    if (Math.random()>rate) continue;
                    sb=new StringBuilder();
                    sb.append(this.caculateDistanceBetweenSNPs(i, j)).append("\t");
                    sb.append(this.getR2(i, j));
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            System.out.println(new File(physicalDistanceR2File).getName()+" r2 < 0: "+ NumberTool.format(cnt/total,
                    5)+"("+(int)cnt+", "+total+")"+ " among r2 great than 0, only "+rate+" were retained");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * Distance	r2
     * 4	1.0
     * 5	0.00302003
     * 6	0.835544
     * 7	1.0
     * 7	0.01
     * 8	1.0
     * 8	0.0357143
     * @param distanceThresh
     * @param r2Thresh
     * @param inputDistanceR2File
     * @param outFile
     */
    public static void getDistanceLessThan(int distanceThresh, double r2Thresh, String inputDistanceR2File,
                                           String outFile){
        try (BufferedReader br = IOUtils.getTextReader(inputDistanceR2File);
             BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            String line;
            line=br.readLine();
            bw.write(line);
            bw.newLine();
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                if(Double.parseDouble(temp.get(1))<r2Thresh) continue;
                if (Integer.parseInt(temp.get(0))<distanceThresh){
                    bw.write(line);
                    bw.newLine();
                }else {
                    break;
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Distance	r2
     * 1	0.146269
     * 1	0.102565
     * 1	0.0314594
     * 1	0.0635227
     * @param inputDistanceR2File
     * @param windowSize
     * @param stepSize
     * @param outFile
     */
    public static void slidingWindow(String inputDistanceR2File, int windowSize, int stepSize, String outFile){
        try (BufferedReader br = IOUtils.getTextReader(inputDistanceR2File);
             BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            br.readLine();
            bw.write("WindowFromTo\tWindowMiddleDistance\tCountInWindow\tWindowAverageR2\tSDofWindowR2");
            bw.newLine();
            double[] distancesArray;
            double[] r2Array;
            String line;
            TDoubleArrayList distanceList=new TDoubleArrayList();
            TDoubleArrayList r2List=new TDoubleArrayList();
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                distanceList.add(Integer.parseInt(temp.get(0)));
                r2List.add(Double.parseDouble(temp.get(1)));
            }
            distancesArray=distanceList.toArray();
            r2Array=r2List.toArray();
            int maxDistace=(int)distancesArray[distancesArray.length-1];
            int num=(maxDistace-windowSize)/stepSize+2;
            double[] boundaryS=IntStream.iterate(0, n->n+stepSize).limit(num).mapToDouble(n->n-0.1).toArray();
            double[] boundaryL=IntStream.iterate(windowSize, n->n+stepSize).limit(num).mapToDouble(n->n-0.1).toArray();
            int[] indexS=new int[boundaryS.length];
            int[] indexL=new int[boundaryL.length];
            for (int i = 0; i < indexL.length; i++) {
                indexS[i]= Integer.MIN_VALUE;
                indexL[i]=Integer.MIN_VALUE;
            }
            int index=Integer.MAX_VALUE;
            for (int i = 0; i < boundaryS.length; i++) {
                index=Arrays.binarySearch(distancesArray, boundaryS[i]);
                indexS[i]=-index-1;
                index=Arrays.binarySearch(distancesArray, boundaryL[i]);
                indexL[i]=-index-1;
            }
            int countInBoundary=-1;
            StringBuilder sb;
            List<Double> doubleList;
            int count=0;
            DescriptiveStatistics stats;
            for (int j = 0; j < indexS.length; j++) {
                countInBoundary=indexL[j]-indexS[j];
                if (countInBoundary<1){
                    count++;
                }
                sb=new StringBuilder();
                stats=new DescriptiveStatistics();
                for (int i = indexS[j]; i < indexL[j]; i++) {
                    stats.addValue(r2Array[i]);
                }
                sb.append((int)(boundaryS[j]+0.1)).append("_").append((int)(boundaryL[j]+0.1)).append("\t");
                sb.append((double) windowSize/2+boundaryS[j]+0.1).append("\t");
                sb.append(countInBoundary).append("\t");
                if (countInBoundary<1) {
                    sb.append("NaN").append("\t").append("NaN").append("\n");
                }else {
                    sb.append(stats.getMean()).append("\t").append(stats.getStandardDeviation()).append("\n");
                }
                bw.write(sb.toString());
            }
            if(count==0){
                System.out.println(count+" window count were 0");
            }else {
                System.out.println(count+" windows count were 0, and its r2 and SD will be setting NaN");
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * inputDir is plink result dir, using .ld and .bim file
     * @param ld_bimDir
     * @param distanceR2OutDir
     * @param numThreads
     * @param rate
     */
    public static void getDistaceR2(String ld_bimDir, String distanceR2OutDir, int numThreads, double rate){
        File[] files=new File(ld_bimDir).listFiles();
        File[] ldFiles=IOUtils.listFilesEndsWith(files, "ld");
        File[] bimFiles=IOUtils.listFilesEndsWith(files, "bim");
        Predicate<File> p=File::isHidden;
        File[] lds=Arrays.stream(ldFiles).filter(p.negate()).sorted().toArray(File[]::new);
        File[] bims=Arrays.stream(bimFiles).filter(p.negate()).sorted().toArray(File[]::new);
        String[] outNames= Arrays.stream(lds).map(File::getName).map(str->str.replaceAll("ld$", "distance_r2.txt"))
                .toArray(String[]::new);
        files=null;
        ldFiles=null;
        bimFiles=null;
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(lds.length, numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream()
                    .forEach(index-> {
                        LD ld=new LD(lds[index], bims[index]);
                        ld.writePhysicalDistanceR2ForLDDecay(new File(distanceR2OutDir, outNames[index]).getAbsolutePath(), rate);
                    });
        }
        File[] distanceR2Files=IOUtils.listRecursiveFiles(new File(distanceR2OutDir));
        BufferedReader br;
        BufferedWriter bw=IOUtils.getTextWriter(new File(distanceR2OutDir, "chrAll.distaceR2.txt").getAbsolutePath());
        try {
            String line;
            bw.write("Distance\tr2\n");
            for (int i = 0; i < distanceR2Files.length; i++) {
                br=IOUtils.getTextReader(distanceR2Files[i].getAbsolutePath());
                br.readLine();
                while ((line=br.readLine())!=null){
                    bw.write(line);
                    bw.newLine();
                }
                br.close();
//                distanceR2Files[i].delete();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        RowTableTool rowTableTool= new RowTableTool(new File(distanceR2OutDir, "chrAll.distaceR2.txt").getAbsolutePath());
        Comparator<List<String>> c=Comparator.comparing(l->Integer.parseInt(l.get(0)));
        rowTableTool.sortBy(c);
        rowTableTool.writeTextTable(new File(distanceR2OutDir, "chrAll.distaceR2.txt").getAbsolutePath(), IOFileFormat.Text);
    }

    /**
     * 多个文件合并
     * @param slidingWindowDistanceR2Dir slidingWindow dir
     * @param group
     * @param outFile
     */
    public static void mergeDistanceR2ForGGplot(String slidingWindowDistanceR2Dir, String[] group, String outFile){
        File[] files=new File(slidingWindowDistanceR2Dir).listFiles();
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files).filter(p.negate().and(File::isFile)).sorted().toArray(File[]::new);
        BufferedReader br;
        BufferedWriter bw=IOUtils.getTextWriter(outFile);
        String line;
        StringBuilder sb;
        String header;
        boolean first=true;
        try {
            for (int i = 0; i < f1.length; i++) {
                br=IOUtils.getTextReader(f1[i].getAbsolutePath());
                header=br.readLine();
                if (first){
                    bw.write(header+"\tGroup\n");
                    first=false;
                }
                while ((line=br.readLine())!=null){
                    sb=new StringBuilder();
                    sb.append(line).append("\t").append(group[i]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void getSubgenomeLD(String subgenomeDir, String subgenome, double rate, String outDir){
        File[] files=Arrays.stream(new File(subgenomeDir).listFiles()).sorted().toArray(File[]::new);
        Predicate<File> h=File::isHidden;
        Predicate<File> s=f->f.getName().substring(4,5).equals(WheatLineage.valueOf(subgenome).name());
        Predicate<File> r=f->Math.random()<rate;
        File[] f=Arrays.stream(files).filter(h.negate().and(s)).filter(r).toArray(File[]::new);
        String[] outNames=Arrays.stream(f).map(File::getName).toArray(String[]::new);
        try {
            for (int i = 0; i < outNames.length; i++) {
                Files.copy(f[i].toPath(), Paths.get(outDir, outNames[i]));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void mergeToSubgenome(String input1Dir, String inout2Dir, String mergedFile){
        File[] files1=new File(input1Dir).listFiles();
        File[] files2=new File(inout2Dir).listFiles();
        File[] mergedFiles= (File[]) ArrayUtils.addAll(files1, files2);
        BufferedReader br;
        BufferedWriter bw;
        try {
            bw=IOUtils.getTextWriter(mergedFile);
            bw.write("Distance\tR2\n");
            String line;
            for (int i = 0; i < mergedFiles.length; i++) {
                br=IOUtils.getTextReader(mergedFiles[i].getAbsolutePath());
                br.readLine();
                while ((line=br.readLine())!=null){
                    bw.write(line);
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        RowTableTool<String> rowTable=new RowTableTool<>(mergedFile);
        Comparator<List<String>> comparator= Comparator.comparing(l->Integer.parseInt(l.get(0)));
        rowTable.sortBy(comparator);
        rowTable.writeTextTable(mergedFile, IOFileFormat.Text);
    }

    private class SNP {

        int chrID;
        int pos;

        SNP(int chrID, int pos){
            this.chrID=chrID;
            this.pos=pos;
        }

        public int getChrID() {
            return chrID;
        }

        public int getPos() {
            return pos;
        }

        public String getChr(int chrID){
            int[] integers= IntStream.range(1,8).toArray();
            String[] abd={"A","B","D"};
            Map<Integer, String> chrIDmap=new HashMap<>();
            int count=0;
            for (int i = 0; i < integers.length; i++) {
                for (int j = 0; j < abd.length; j++) {
                    count++;
                    chrIDmap.put(count, integers[i]+abd[j]);
                }
            }
            return chrIDmap.get(chrID);
        }

        @Override
        public String toString() {
            StringBuilder sb=new StringBuilder();
            sb.append(this.getChr(chrID)).append("-").append(pos);
            return sb.toString();
        }
    }

//    public static void main(String[] args) {
//        LD ld=new LD(new File(args[0]), new File(args[1]));
//        ld.writePhysicalDistanceR2ForLDDecay(args[3], Integer.parseInt(args[4]));
//        Comparator<List<String>> c=Comparator.comparing(l->Integer.parseInt(l.get(0)));
//        RowTableTool<String> rowTable=new RowTableTool<>(args[3]);
//        rowTable.sortBy(c);
//        rowTable.writeTextTable(args[5], IOFileFormat.Text);
//        LD.getDistanceLessThan(Integer.parseInt(args[11]), Double.parseDouble(args[12]), args[5], args[6]);
//        LD.slidingWindow(args[6], Integer.parseInt(args[13]), Integer.parseInt(args[14]), args[7]);
//    }
}