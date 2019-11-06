package daxing.applets;

import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.math3.util.CombinatoricsUtils;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;
import java.util.stream.IntStream;

public class LD {

    List<TDoubleArrayList> matrix;
    List<SNP> header;

    public LD(String ldMatrixFile, String bimFile){
        List<TDoubleArrayList> matrix=new ArrayList<>();
        List<SNP> snpList=new ArrayList<>();
        int i=0;
        try (BufferedReader br1 = IOUtils.getTextReader(ldMatrixFile);
             BufferedReader br2=IOUtils.getTextReader(bimFile)) {
            String line;
            String[] temp;
            TDoubleArrayList t;
            SNP snp;
            while ((line=br1.readLine())!=null){
                temp= StringUtils.split(line, " ");
                t=new TDoubleArrayList();
                for (i = 0; i < temp.length; i++) {
                    if (temp[i].equals("nan")){
                        t.add(-1d);
                    }else {
                        t.add(Double.parseDouble(temp[i]));
                    }

                }
                matrix.add(t);
            }
           while ((line=br2.readLine())!=null){
               temp= StringUtils.split(line, "\t");
               snp=new SNP(Integer.parseInt(temp[0]), Integer.parseInt(temp[3]));
               snpList.add(snp);
           }
           this.matrix=matrix;
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
        return this.matrix.get(snpIndex1).get(snpIndex2);
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

    public static void slidingWindow(String inputDistanceR2File, int windowSize, int stepSize, String outFile){
        try (BufferedReader br = IOUtils.getTextReader(inputDistanceR2File);
             BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            br.readLine();
            bw.write("WindowFromTo\tWindowMiddleDistance\tCountInWindow\tWindowAverageR2");
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
            double sumOfR2=0;
            int count=0;
            int cnt=0;
            for (int j = 0; j < indexS.length; j++) {
                countInBoundary=indexL[j]-indexS[j];
                if (countInBoundary<1){
                    count++;
                }
                sb=new StringBuilder();
                for (int i = indexS[j]; i < indexL[j]; i++) {
                    sumOfR2+=r2Array[i];
                }
                sb.append((int)(boundaryS[j]+0.1)).append("_").append((int)(boundaryL[j]+0.1)).append("\t");
                sb.append((double) windowSize/2+boundaryS[j]+0.1).append("\t");
                sb.append(countInBoundary).append("\t");
                sb.append(sumOfR2/countInBoundary).append("\n");
                bw.write(sb.toString());
                sumOfR2=0;
            }
            System.out.println(cnt+" r2 great than 100000000 ");
            if(count==0){
                System.out.println(count+" window count were 0");
            }else {
                System.out.println(count+" windows count were 0, and its r2 will be setting NaN");
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
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
}
