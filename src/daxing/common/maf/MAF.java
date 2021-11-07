package daxing.common.maf;

import daxing.common.factors.WheatLineage;
import daxing.common.utiles.IOTool;
import org.apache.commons.lang.math.NumberUtils;
import org.apache.commons.lang3.StringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * multiple alignment format
 */
public class MAF {

    List<MAFAlignment> mafAlignmentList;
    WheatLineage wheatLineage;

    public MAF(File mafFile, WheatLineage wheatLineage){
        this.wheatLineage=wheatLineage;
        List<MAFAlignment> mafAlignmentList= new ArrayList<>();
        MAFAlignment mafAlignment;
        List<String> linesPerAlignment= new ArrayList<>();
        String line;
        try (BufferedReader br = IOTool.getReader(mafFile)) {
            while ((line=br.readLine()).startsWith("##")){}
            linesPerAlignment.add(line);
            while ((line=br.readLine())!=null){
                if (org.apache.commons.lang.StringUtils.isBlank(line)){
                    mafAlignment = new MAFAlignment(linesPerAlignment);
                    mafAlignmentList.add(mafAlignment);
                    linesPerAlignment= new ArrayList<>();
                }else {
                    linesPerAlignment.add(line);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        this.mafAlignmentList=mafAlignmentList;
    }

    public List<MAFAlignment> getMafAlignmentList() {
        return mafAlignmentList;
    }

    public void extractSpecies(String outFile, EnumSet<Species> speciesEnumSet){
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            List<Species> speciesList = new ArrayList<>(speciesEnumSet);
            Collections.sort(speciesList);
            StringBuilder sb = new StringBuilder();
            int refChr, refPos_1_Based;
            Record refRecord;
            sb.append("Chr\tPos\t");
            for (int i = 0; i < speciesList.size(); i++) {
                sb.append(speciesList.get(i).getStr()).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            byte baseAscII;
            for (MAFAlignment mafAlignment : this.getMafAlignmentList()){
                if (mafAlignment.ifExistUnChr(speciesEnumSet)) continue;
                for (int i = 0; i < mafAlignment.alignmentSize; i++) {
                    if (mafAlignment.anyDash(speciesEnumSet, i)) continue;
                    refRecord = mafAlignment.speciesMAFRecordEnumMap.get(Species.traes);
                    refChr = refRecord.refChr;
                    refPos_1_Based= refRecord.get1basedRefPosWithoutDash(i);
                    sb.setLength(0);
                    sb.append(refChr).append(wheatLineage).append("\t");
                    sb.append(refPos_1_Based).append("\t");
                    for (Species species : speciesList){
                        baseAscII=mafAlignment.speciesMAFRecordEnumMap.get(species).baseAscIIArray[i];
                        sb.append(Sequence.getBaseFromAscII(baseAscII)).append("\t");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static class Record {
        int refChr;
        int start;
        int size;
        Strand strand;
        int srcSize;
        byte[] baseAscIIArray;

        Record(int refChr, int start, int size, Strand strand, int srcSize, byte[] ascIIArray){
            this.refChr=refChr;
            this.start=start;
            this.size=size;
            this.strand=strand;
            this.srcSize=srcSize;
            this.baseAscIIArray =ascIIArray;
        }

        public boolean isDash(int indexInBaseAscIIArray){
            return baseAscIIArray[indexInBaseAscIIArray]==45;
        }

        public int getUCSCStartPosWithoutDash(int indexInBaseAscIIArray){
            int ucscPos= this.start;
            for (int i = 0; i < indexInBaseAscIIArray; i++) {
                if (isDash(i)) continue;
                ucscPos++;
            }
            return ucscPos;
        }

        public int get1basedRefPosWithoutDash(int indexInBaseAscIIArray){
            int ucscPos=this.getUCSCStartPosWithoutDash(indexInBaseAscIIArray);
            Coordinate coordinate = Coordinate.getInstance();
            coordinate.setFrom(strand, ucscPos, srcSize);
            return coordinate.getOneStart();
        }
    }

    public static class MAFAlignment{
        double score;
        int alignmentSize;
        EnumMap<Species, Record> speciesMAFRecordEnumMap= new EnumMap<>(Species.class);

        MAFAlignment(List<String> linesPerAlignment){
            String[] temp, tem;
            score=Double.parseDouble(StringUtils.split(linesPerAlignment.get(0), "=")[1]);
            int refChr, start, size, srcSize, alignmentSize=-1;
            Strand strand;
            byte[] ascIIArray;
            Species species;
            Record record;
            for (int i = 1; i < linesPerAlignment.size(); i++) {
                temp= StringUtils.split(linesPerAlignment.get(i), " ");
                tem=StringUtils.split(temp[1], ".");
                // exist MT chr
                species= i==1 ? Species.valueOf(tem[0].substring(0,5)) : Species.valueOf(tem[0]);
                refChr= tem[1].length() < 4 ? -1 : NumberUtils.isNumber(tem[1].substring(3)) ?
                    Integer.parseInt(tem[1].substring(3)) : -1;
                start=Integer.parseInt(temp[2]);
                size = Integer.parseInt(temp[3]);
                strand=Strand.getInstanceFrom(temp[4]);
                srcSize=Integer.parseInt(temp[5]);
                alignmentSize=temp[6].length();
                ascIIArray=temp[6].getBytes();
                record= new Record(refChr, start, size, strand, srcSize, ascIIArray);
                speciesMAFRecordEnumMap.put(species, record);
            }
            this.alignmentSize=alignmentSize;
        }

        /**
         *  指定的species 是否存在叶绿体、线粒体、unknow Chr
         */
        public boolean ifExistUnChr(EnumSet<Species> speciesEnumSet){
            for (Species species: speciesEnumSet){
                if (speciesMAFRecordEnumMap.get(species).refChr < 0) return true;
            }
            return false;
        }

        /**
         * 指定的species 和 index 是否存在 dash(-)
         */
        public boolean anyDash(EnumSet<Species> speciesEnumSet, int indexInBaseAscIIArray){
            for (Species species : speciesEnumSet){
                if(speciesMAFRecordEnumMap.get(species).isDash(indexInBaseAscIIArray)) return true;
            }
            return false;
        }

    }


}
