package daxing.applets;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
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
            String header=br.readLine();
            bw.write(header);
            bw.newLine();
            int[] distancesArray;
            double[] r2Array;
            String line;
            TIntArrayList distanceList=new TIntArrayList();
            TDoubleArrayList r2List=new TDoubleArrayList();
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                distanceList.add(Integer.parseInt(temp.get(0)));
                r2List.add(Double.parseDouble(temp.get(1)));
            }
            distancesArray=distanceList.toArray();
            r2Array=r2List.toArray();
            int maxDistace=distancesArray[distancesArray.length-1];
            int num=(maxDistace-windowSize)/stepSize+2;
            int[] boundaryS=IntStream.iterate(0, n->n+stepSize).limit(num).toArray();
            int[] boundaryL=IntStream.iterate(windowSize, n->n+stepSize).limit(num).toArray();
            TIntArrayList indexS=new TIntArrayList();
            TDoubleArrayList newR2List=new TDoubleArrayList();

            double sum=0d;
            int k=0;
            int count=0;
            for (int i = 0; i < boundaryL.length; i++) {
                for (int j = k; j < r2Array.length; j++) {
                    int distance=distancesArray[j];
                    int index=boundaryL[i];
                    if (distancesArray[j]<boundaryL[i]){
                        sum+=r2Array[j];
                        count++;
                    }else {
                        bw.write(sum+"\t"+count);
                        bw.newLine();
                        newR2List.add(sum/count);
                        sum=r2Array[j];
                        count=1;
                        k=j+1;
                        break;
                    }
                }
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
