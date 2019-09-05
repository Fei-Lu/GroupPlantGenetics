package daxing.ancestralAllele;

import daxing.common.ChrConvertionRule;
import format.position.ChrPos;
import gnu.trove.list.array.TIntArrayList;
import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * MAF类按0-44划分，但mafRecord下的position是0-based coordinates
 */
public class MAF {

    private List<MAFrecord>[] mafRecords;
    private int numThreads=45;
    private int taxonIndexForOrder;   //triticum_aestivum在MAFrecord记录中的位置，上为0，下为1
    private ChrConvertionRule chrConvertionRule;
    private static int count=0;

    /**
     *
     * @param taxonIndexForOrder must be in format like "hordeum_vulgare" or "triticum_aestivum"
     * @param chrConvertionRule
     * @param mafInputFileDir
     */
    public MAF(int taxonIndexForOrder, ChrConvertionRule chrConvertionRule, Path mafInputFileDir){
        this.taxonIndexForOrder = taxonIndexForOrder;
        this.chrConvertionRule=chrConvertionRule;
        this.initialize(mafInputFileDir);
        count++;
    }

    private void initialize(Path mafInputFileDir){
        File[] files=IOUtils.listRecursiveFiles(new File(mafInputFileDir.toString()));
        Predicate<File> p=File::isHidden;
        File[] fileArray= Arrays.stream(files).filter(p.negate()).filter(e->e.getName().endsWith("maf")).toArray(File[]::new);
        List<MAFrecord>[] mafrecord=new ArrayList[45];
        for (int i = 0; i < mafrecord.length; i++) {
            mafrecord[i]=new ArrayList<>((fileArray.length)*200/42);
        }
        ConcurrentHashMap<MAFrecord, File> maFrecordsMap=new ConcurrentHashMap<>();
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(fileArray.length, this.getNumThreads());
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream()
                    .forEach(index-> {
                        try(BufferedReader br= IOUtils.getNIOTextReader(fileArray[index].toString())){
                            br.readLine();
                            br.readLine();
                            long id=Long.MIN_VALUE;
                            long score=Long.MIN_VALUE;
                            String[] taxons;
                            String[] chr;
                            int[] startPos;
                            int[] seqLen;
                            boolean[] ifMinus;
                            int[] chrLen;
                            SeqByte[] seq;
                            String line;
                            List<String> lineList1, lineList2, lineList3, lineList4;
                            while ((line=br.readLine())!=null){
                                taxons=new String[2];
                                chr=new String[2];
                                startPos=new int[2];
                                seqLen=new int[2];
                                ifMinus=new boolean[2];
                                chrLen=new int[2];
                                seq=new SeqByte[2];
                                lineList1= PStringUtils.fastSplit(line, " ");
                                id=Long.parseLong(lineList1.get(2));
                                lineList2=PStringUtils.fastSplit(br.readLine(), "=");
                                score=Long.parseLong(lineList2.get(1));
                                lineList3=PStringUtils.fastSplit(br.readLine(), " ");
                                Predicate<String> predicate=s->s.equals("");
                                lineList3=lineList3.stream().filter(predicate.negate()).collect(Collectors.toList());
                                taxons[0]=PStringUtils.fastSplit(lineList3.get(1), ".").get(0);
                                chr[0]=PStringUtils.fastSplit(lineList3.get(1), ".").get(1);
                                startPos[0] =Integer.parseInt(lineList3.get(2));
                                seqLen[0] = Integer.parseInt(lineList3.get(3));
                                if (lineList3.get(4).equals("+")){
                                    ifMinus[0]=false;
                                }else {
                                    ifMinus[0]=true;
                                }
                                chrLen[0]=Integer.parseInt(lineList3.get(5));
                                seq[0]=new SeqByte(lineList3.get(6));
                                lineList4=PStringUtils.fastSplit(br.readLine(), " ");
                                lineList4=lineList4.stream().filter(predicate.negate()).collect(Collectors.toList());
                                taxons[1]=PStringUtils.fastSplit(lineList4.get(1), ".").get(0);
                                chr[1]=PStringUtils.fastSplit(lineList4.get(1), ".").get(1);
                                startPos[1] =Integer.parseInt(lineList4.get(2));
                                seqLen[1] = Integer.parseInt(lineList4.get(3));
                                if (lineList4.get(4).equals("+")){
                                    ifMinus[1]=false;
                                }else {
                                    ifMinus[1]=true;
                                }
                                chrLen[1]=Integer.parseInt(lineList4.get(5));
                                seq[1]=new SeqByte(lineList4.get(6));
                                maFrecordsMap.put(new MAFrecord(id, score, taxons, chr, startPos, seqLen, ifMinus,
                                        chrLen, seq), fileArray[index]);
                                br.readLine();
                            }
                        }catch (Exception e){
                            e.printStackTrace();
                        }
                    });
        }
        List<MAFrecord> mafRecordList=new ArrayList<>(maFrecordsMap.keySet());
        int chrID;
        String chr;
        int startPos;
        for(MAFrecord e: mafRecordList){
            chr=e.getChr(this.getTaxonIndexForOrder());
            startPos=e.getStartPos(this.getTaxonIndexForOrder());
            chrID=this.getChrConvertionRule().getChrIDFromOriChrName(chr, startPos);
            mafrecord[chrID].add(e);
        }
        Comparator<MAFrecord> comparator=Comparator.comparing(m->m.getStartPos(this.getTaxonIndexForOrder()));
        Arrays.stream(mafrecord).forEach(e-> Collections.sort(e, comparator));
        this.mafRecords=mafrecord;
    }

    public List<MAFrecord>[] getMafRecords() {
        return mafRecords;
    }

    public int getNumThreads() {
        return numThreads;
    }

    public int getTaxonIndexForOrder(){
        return taxonIndexForOrder;
    }

    public ChrConvertionRule getChrConvertionRule(){
        return chrConvertionRule;
    }

    public void setNumOfThreads(int numThreads){
        this.numThreads=numThreads;
    }

    public Map<ChrPos, String[]> getAllele(AllelesInfor allelesInfor){
        ConcurrentHashMap<ChrPos, String[]> chrPosOutgroupAlleleMap=new ConcurrentHashMap<>();
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetSize(this.getMafRecords().length, numThreads);
        for (int i = 0; i < indices.length; i++) {
            int[] subLibIndices = new int[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            Arrays.stream(subLibIndices).parallel().forEach(e->{
                List<MAFrecord> maFrecordList=this.getMafRecords()[e];
                List<ChrPos> chrPosList=allelesInfor.getChrPoss((short) e);
                List<String> refAlleleList=allelesInfor.getRefAllele(((short)e));
                ChrPos chrPos;
                String refAllele;
                String outGroupAllele;
                String[] refOutgroupAllele;
                int index;
                for (int j = 0; j < chrPosList.size(); j++) {
                    chrPos=chrPosList.get(j);
                    refAllele=refAlleleList.get(j);
                    index=this.binarySearch(chrPos);
                    if (index<0) continue;
                    outGroupAllele=maFrecordList.get(index).getOutgroupAllele(this.getTaxonIndexForOrder(), chrPos, this.chrConvertionRule);
                    refOutgroupAllele=new String[2];
                    refOutgroupAllele[0]=refAllele;
                    refOutgroupAllele[1]=outGroupAllele;
                    chrPosOutgroupAlleleMap.put(chrPos, refOutgroupAllele);
                }
            });
        }
        return chrPosOutgroupAlleleMap;
    }

    public void getAllele(AllelesInfor allelesInfor, File outDir){
        ConcurrentHashMap<ChrPos, String[]> chrPosOutgroupAlleleMap=new ConcurrentHashMap<>();
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetSize(this.getMafRecords().length, numThreads);
        for (int i = 0; i < indices.length; i++) {
            int[] subLibIndices = new int[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            Arrays.stream(subLibIndices).parallel().forEach(e->{
                List<MAFrecord> maFrecordList=this.getMafRecords()[e];
                List<ChrPos> chrPosList=allelesInfor.getChrPoss((short) e);
                ChrPos chrPos;
                String refAllele;
                String outGroupAllele;
                String[] refOutgroupAllele;
                int index;
                for (int j = 0; j < chrPosList.size(); j++) {
                    chrPos=chrPosList.get(j);
                    index=this.binarySearch(chrPos);
                    if (index<0) continue;
                    refAllele=maFrecordList.get(index).getRefAllele(this.getTaxonIndexForOrder(), chrPos, this.chrConvertionRule);
                    outGroupAllele=maFrecordList.get(index).getOutgroupAllele(this.getTaxonIndexForOrder(), chrPos, this.chrConvertionRule);
                    refOutgroupAllele=new String[2];
                    refOutgroupAllele[0]=refAllele;
                    refOutgroupAllele[1]=outGroupAllele;
                    chrPosOutgroupAlleleMap.put(chrPos, refOutgroupAllele);
                }
            });
        }
        List<Integer> l=new ArrayList<>();
        l.add(0);
        l.add(1);
        l.remove(this.getTaxonIndexForOrder());
        String[] taxons=this.getMafRecords()[1].get(1).getTaxon();
        String taxon=taxons[l.get(0)];
        try(BufferedWriter bw=IOUtils.getTextWriter(new File(outDir, taxon+".txt").toString())){
            StringBuilder sb=new StringBuilder();
            sb.append("CHR").append("\t").append("POS").append("\t")
                    .append("refAllele").append("\t").append(taxon).append("\n");
            bw.write(sb.toString());
            ChrPos key;
            String[] value;
            for (Map.Entry<ChrPos, String[]> entry: chrPosOutgroupAlleleMap.entrySet()){
                key=entry.getKey();
                value=entry.getValue();
                sb=new StringBuilder();
                sb.append(key.getChromosome()).append("\t").append(key.getPosition()).append("\t")
                        .append(value[0]).append("\t").append(value[1]).append("\n");
                bw.write(sb.toString());
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public void getAllele(String outDir){
        Map<Integer, BufferedWriter> integerBufferedWriterMap=new HashMap<>();
        BufferedWriter bw;
        List<Integer> l=new ArrayList<>();
        l.add(0);
        l.add(1);
        l.remove(this.getTaxonIndexForOrder());
        String[] taxons=this.getMafRecords()[1].get(1).getTaxon();
        String outgroupTaxon=taxons[l.get(0)];
        for (int i = 0; i < this.getMafRecords().length; i++) {
            bw=IOUtils.getTextWriter(new File(outDir, "triticum_aestivumChr"+PStringUtils.getNDigitNumber(3, i)+"_v_"+outgroupTaxon+".txt").getAbsolutePath());
            integerBufferedWriterMap.put(i, bw);
        }
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetSize(this.getMafRecords().length, numThreads);
        for (int i = 0; i < indices.length; i++) {
            int[] subLibIndices = new int[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            Arrays.stream(subLibIndices).parallel().forEach(e->{
                List<MAFrecord> maFrecordList=this.getMafRecords()[e];
                StringBuilder sb=new StringBuilder();
                try(BufferedWriter bufferedWriter=integerBufferedWriterMap.get(e)) {
                    sb.append("CHR").append("\t").append("POS").append("\t").append("ref").append("\t").append(outgroupTaxon).append("\n");
                    bufferedWriter.write(sb.toString());
                    Map<ChrPos, String[]> refOutAlleleMap;
                    List<ChrPos> list;
                    for (int j = 0; j < maFrecordList.size(); j++) {
                        refOutAlleleMap=maFrecordList.get(j).getRefCoordinateOfOutGroup(this.getTaxonIndexForOrder(), this.getChrConvertionRule());
                        list=new ArrayList<>(refOutAlleleMap.keySet());
                        Collections.sort(list);
                        for (int k = 0; k < list.size(); k++) {
                            sb=new StringBuilder();
                            sb.append(list.get(k).getChromosome()).append("\t").append(list.get(k).getPosition()).append("\t")
                                    .append(refOutAlleleMap.get(list.get(k))[0]).append("\t").append(refOutAlleleMap.get(list.get(k))[1]);
                            bufferedWriter.write(sb.toString());
                            bufferedWriter.newLine();
                        }
                    }
                    bufferedWriter.flush();
                    System.out.println("triticum_aestivumChr"+PStringUtils.getNDigitNumber(3, e)+"_v_"+outgroupTaxon+".txt is completed");
                }catch (Exception exception){
                    exception.printStackTrace();
                }

            });
        }
    }

    public List<int[]> getStartEnd(short chr){
        List<int[]> list=new ArrayList<>();
        List<MAFrecord> maFrecordList=this.getMafRecords()[chr];
        int[] startEnd;
        for (int i = 0; i < maFrecordList.size(); i++) {
            startEnd=maFrecordList.get(i).startEnd(this.getTaxonIndexForOrder());
            Arrays.sort(startEnd);
            list.add(startEnd);
        }
        return list;
    }

    /**
     * 判断给定的ChrPos是否位于比对内部，如果不在比对内部，则返回-1
     * @param chrPos 来源于vcf文件的ChrPos，是1-based coordinates
     * @return 给定ChrPos在染色体上对应的index
     */
    public int binarySearch(ChrPos chrPos){
        short chr=chrPos.getChromosome();
        int refPos=this.getChrConvertionRule().getRefPositionFromVCF(chrPos);
        List<int[]> startEnd=this.getStartEnd(chr);
        TIntArrayList start=new TIntArrayList();
        TIntArrayList end=new TIntArrayList();
        for (int i = 0; i < startEnd.size(); i++) {
            start.add(startEnd.get(i)[0]);
            end.add(startEnd.get(i)[1]);
        }
        int res=start.binarySearch(refPos);
        int index;
        if (res<-1){
            index=-res-2;
            if (end.get(index)>=refPos){
                return index;
            }else {
                return -1;
            }
        }else {
            return res;
        }
    }

    public static void mergeTwoFiles(String inputOutgroup1File, String inputOutgroup2File, String outFileDir){
        try {
            List<List<String>> l1= Files.newBufferedReader(Paths.get(inputOutgroup1File)).lines().skip(1).parallel()
                    .map(PStringUtils::fastSplit).collect(Collectors.toList());
            List<List<String>> l2=Files.newBufferedReader(Paths.get(inputOutgroup2File)).lines().skip(1).parallel()
                    .map(PStringUtils::fastSplit).collect(Collectors.toList());
            Map<ChrPos, String[]> map1=new HashMap<>();
            Map<ChrPos, String[]> map2=new HashMap<>();
            Map<ChrPos, String[]> map=new HashMap<>();
            short chr;
            int pos;
            String[] refAle;
            for (int i = 0; i < l1.size(); i++) {
                chr=Short.parseShort(l1.get(i).get(0));
                pos=Integer.parseInt(l1.get(i).get(1));
                refAle=new String[2];
                refAle[0]=l1.get(i).get(2);
                refAle[1]=l1.get(i).get(3);
                map1.put(new ChrPos(chr, pos), refAle);
            }
            for (int i = 0; i < l2.size(); i++) {
                chr=Short.parseShort(l2.get(i).get(0));
                pos=Integer.parseInt(l2.get(i).get(1));
                refAle=new String[2];
                refAle[0]=l2.get(i).get(2);
                refAle[1]=l2.get(i).get(3);
                map2.put(new ChrPos(chr, pos), refAle);
            }
            ChrPos key1, key2;
            String[] value1, value2;
            String[] refOut1_2;
            for (Map.Entry<ChrPos, String[]> entry: map1.entrySet()){
                refOut1_2=new String[3];
                key1=entry.getKey();
                chr=key1.getChromosome();
                pos=key1.getPosition();
                value1=entry.getValue();
                if (map2.containsKey(key1)){
                    refOut1_2[0]=value1[0];
                    refOut1_2[1]=value1[1];
                    refOut1_2[2]=map2.get(key1)[1];
                    map.put(new ChrPos(chr, pos), refOut1_2);
                }
                else {
                    refOut1_2[0]=value1[0];
                    refOut1_2[1]=value1[1];
                    refOut1_2[2]="-";  //"-"表示对应的outgroup没有allele
                    map.put(new ChrPos(chr, pos), refOut1_2);
                }
            }
            Set<ChrPos> set2=map2.keySet();
            set2.removeAll(map1.keySet());
            for (Map.Entry<ChrPos, String[]> entry: map2.entrySet()){
                key2=entry.getKey();
                value2=entry.getValue();
                chr=key2.getChromosome();
                pos=key2.getPosition();
                refOut1_2=new String[3];
                refOut1_2[0]=value2[0];
                refOut1_2[1]="-";    //"-"表示对应的outgroup没有allele
                refOut1_2[2]=value2[1];
                map.put(new ChrPos(chr, pos), refOut1_2);
            }
            List<ChrPos> list=new ArrayList<>(map.keySet());
            Collections.sort(list);
            String[] tasons=new String[2];
            List<String> taxonName1List=PStringUtils.fastSplit(new File(inputOutgroup1File).getName(), "_v_");
            List<String> taxonName2List=PStringUtils.fastSplit(new File(inputOutgroup2File).getName(), "_v_");
            tasons[0]=taxonName1List.get(1).replaceAll(".txt$", "");
            tasons[1]=taxonName2List.get(1).replaceAll(".txt$", "");
            BufferedWriter bw=IOUtils.getTextWriter(new File(outFileDir,taxonName1List.get(0)+"_ancestralAllele_"+tasons[0]+"_v_"+tasons[1]+".txt").getAbsolutePath());
            StringBuilder sb;
            sb=new StringBuilder();
            sb.append("CHR").append("\t").append("POS").append("\t").append("refAllele").append("\t")
                    .append(tasons[0]).append("\t").append(tasons[1]).append("\n");
            bw.write(sb.toString());
            for (int i = 0; i < list.size(); i++) {
                chr=list.get(i).getChromosome();
                pos=list.get(i).getPosition();
                refOut1_2=map.get(list.get(i));
                sb=new StringBuilder();
                sb.append(chr).append("\t").append(pos).append("\t").append(refOut1_2[0].toUpperCase())
                        .append("\t").append(refOut1_2[1].toUpperCase()).append("\t").append(refOut1_2[2].toUpperCase()).append("\n");
                bw.write(sb.toString());
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void merge(String inputDir, String mergeOutDir){
        File[] files=IOUtils.listRecursiveFiles(new File(inputDir));
        Predicate<File> p=File::isHidden;
        File[] f=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        for (int i = 0; i < f.length; i=i+2) {
            MAF.mergeTwoFiles(f[i].getAbsolutePath(), f[i+1].getAbsolutePath(), mergeOutDir);
        }
    }

    /**
     *
     * @param inputDir
     * @param outDir
     */
    public static void getAncestralAllele(String inputDir, String outDir){
        File[] input=new File(inputDir).listFiles();
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).sorted().toArray(File[]::new);
        Map<Integer, BufferedReader> chrBufferReaderMap=new HashMap<>();
        Map<Integer, BufferedWriter> chrBufferWriterMap=new HashMap<>();
        IntStream.range(0,45).forEach(e->{
            chrBufferReaderMap.put(e, IOUtils.getNIOTextReader(files[e].getAbsolutePath()));
            chrBufferWriterMap.put(e, IOUtils.getNIOTextWriter(new File(outDir, files[e].getName()).getAbsolutePath()));
        });
        BufferedReader br;
        BufferedWriter bw;
        List<String> lines;
        String line;
        String header;
        List<String> lineList;
        try{
            for (int i = 0; i < files.length; i++) {
                br=chrBufferReaderMap.get(i);
                header=br.readLine();
                lines=new ArrayList<>();
                while ((line=br.readLine())!=null){
                    lineList=PStringUtils.fastSplit(line);
                    if (lineList.get(3).equals("N")) continue;
                    if (lineList.get(3).equals("-")) continue;
                    if (lineList.get(4).equals("N")) continue;
                    if (lineList.get(4).equals("-")) continue;
                    if (!lineList.get(3).equals(lineList.get(4))) continue;
                    lines.add(line);
                }
                br.close();
                bw=chrBufferWriterMap.get(i);
                bw.write(header);
                bw.newLine();
                for (int j = 0; j < lines.size(); j++) {
                    bw.write(lines.get(j));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     *
     * @param inputDir
     * @param outDir
     */
    public static void getAncestralAlleleParallel(String inputDir, String outDir){
        File[] input=new File(inputDir).listFiles();
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).sorted().toArray(File[]::new);
        Map<Integer, BufferedReader> chrBufferReaderMap=new HashMap<>();
        Map<Integer, BufferedWriter> chrBufferWriterMap=new HashMap<>();
        IntStream.range(0,45).forEach(e->{
            chrBufferReaderMap.put(e, IOUtils.getNIOTextReader(files[e].getAbsolutePath()));
            chrBufferWriterMap.put(e, IOUtils.getNIOTextWriter(new File(outDir, files[e].getName()).getAbsolutePath()));
        });
        IntStream.range(0, 45).parallel().forEach(index->{
            BufferedReader br;
            BufferedWriter bw;
            List<String> lines;
            String line;
            String header;
            List<String> lineList;
            try{
                br=chrBufferReaderMap.get(index);
                header=br.readLine();
                lines=new ArrayList<>();
                while ((line=br.readLine())!=null){
                    lineList=PStringUtils.fastSplit(line);
                    if (lineList.get(3).equals("N")) continue;
                    if (lineList.get(3).equals("-")) continue;
                    if (lineList.get(4).equals("N")) continue;
                    if (lineList.get(4).equals("-")) continue;
                    if (!lineList.get(3).equals(lineList.get(4))) continue;
                    lines.add(line);
                }
                br.close();
                bw=chrBufferWriterMap.get(index);
                bw.write(header);
                bw.newLine();
                for (int j = 0; j < lines.size(); j++) {
                    bw.write(lines.get(j));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }catch (Exception e){
                e.printStackTrace();
            }
        });
    }


}
