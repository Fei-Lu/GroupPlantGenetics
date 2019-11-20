package daxing.ancestralAllele;

import daxing.common.ChrConvertionRule;
import daxing.common.WheatLineage;
import format.position.ChrPos;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TShortHashSet;
import utils.Benchmark;
import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Path;
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
     * @param taxonIndexForOrder triticum_aestivum在MAFrecord记录中的位置，上为0，下为1
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
        List<String> d_lineage= WheatLineage.valueOf("D").getChr();
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
            integerList.stream().parallel()
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
                                br.readLine();
                                if (d_lineage.contains(chr[taxonIndexForOrder])) continue;
                                maFrecordsMap.put(new MAFrecord(id, score, taxons, chr, startPos, seqLen, ifMinus,
                                        chrLen, seq), fileArray[index]);
                            }
                        }catch (Exception e){
                            e.printStackTrace();
                        }
                    });
        }
        List<MAFrecord> mafRecordList=new ArrayList<>(maFrecordsMap.keySet());
        int chrID;
        String chr;
        int startPos1_based;
        for(MAFrecord e: mafRecordList){
            chr=e.getChr(this.getTaxonIndexForOrder());
            startPos1_based=e.getStartPos1_based(this.getTaxonIndexForOrder());
            chrID=this.getChrConvertionRule().getVCFChrFromRefChrPos(chr, startPos1_based);
            mafrecord[chrID].add(e);
        }
        Comparator<MAFrecord> comparator=Comparator.comparing(m->m.getStartPos1_based(this.getTaxonIndexForOrder()));
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
        String taxon=this.getAnotherTaxonName();
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
        long start=System.nanoTime();
        ConcurrentHashMap<Integer, BufferedWriter> integerBufferedWriterMap=new ConcurrentHashMap<>();
        BufferedWriter bw;
        String outgroupTaxon=this.getAnotherTaxonName();
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
                    System.out.println("triticum_aestivumChr"+PStringUtils.getNDigitNumber(3, e)+"_v_"+outgroupTaxon+".txt is completed in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
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
        int refPos=this.getChrConvertionRule().getRefPosFromVCFChrPos(chrPos);
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

    public String getAnotherTaxonName(){
        List<Integer> l=new ArrayList<>();
        l.add(0);
        l.add(1);
        l.remove(this.getTaxonIndexForOrder());
        String[] taxons=this.getMafRecords()[1].get(1).getTaxon();
        String outgroupTaxon=taxons[l.get(0)];
        return outgroupTaxon;
    }

    public static void mergeTwoFiles(File inputOutgroup1File, File inputOutgroup2File, File outFile){
        try {
            long start=System.nanoTime();
            BufferedReader br1=IOUtils.getTextReader(inputOutgroup1File.getAbsolutePath());
            BufferedReader br2=IOUtils.getTextReader(inputOutgroup2File.getAbsolutePath());
            BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath());
            br1.readLine();
            br2.readLine();
            TShortHashSet chr=new TShortHashSet();
            String line;
            List<String> lineList;
            TIntArrayList posList1=new TIntArrayList();
            TIntArrayList posList2=new TIntArrayList();
            List<String> refList1=new ArrayList<>();
            List<String> altList1=new ArrayList<>();
            List<String> altList2=new ArrayList<>();
            while ((line=br1.readLine())!=null){
                lineList = PStringUtils.fastSplit(line);
                chr.add(Short.parseShort(lineList.get(0)));
                posList1.add(Integer.parseInt(lineList.get(1)));
                refList1.add(lineList.get(2).toUpperCase());
                altList1.add(lineList.get(3).toUpperCase());
            }
            br1.close();
            while ((line=br2.readLine())!=null){
                lineList = PStringUtils.fastSplit(line);
                chr.add(Short.parseShort(lineList.get(0)));
                posList2.add(Integer.parseInt(lineList.get(1)));
                altList2.add(lineList.get(3).toUpperCase());
            }
            br2.close();
            if (chr.size()>1){
                System.out.println("Wrong input files!"+"\t"+inputOutgroup1File.getName()+"\t"+inputOutgroup2File.getName()+" Program quits.");
//                System.exit(0);
                return;
            }
            TIntHashSet mergedPosSet=new TIntHashSet(posList1);
            mergedPosSet.addAll(posList2);
            int[] mergedPos=mergedPosSet.toArray();
            Arrays.sort(mergedPos);
            int index1=Integer.MIN_VALUE;
            int index2=Integer.MIN_VALUE;
            StringBuilder sb=new StringBuilder();
            sb.append("CHR").append("\t").append("POS").append("\t").append("refAllele").append("\t")
                    .append("outgroup1").append("\t").append("outgroup2").append("\n");
            bw.write(sb.toString());
            for (int i = 0; i < mergedPos.length; i++) {
                index1=posList1.binarySearch(mergedPos[i]);
                index2=posList2.binarySearch(mergedPos[i]);
                if (index1 < 0 || index2 < 0) continue;
                sb=new StringBuilder();
                sb.append(chr.iterator().next()).append("\t").append(mergedPos[i]).append("\t").append(refList1.get(index1))
                        .append("\t").append(altList1.get(index1)).append("\t").append(altList2.get(index2));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(outFile.getName()+" is completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void merge(String inputDir, String mergeOutDir){
        File[] files=IOUtils.listRecursiveFiles(new File(inputDir));
        Predicate<File> p=File::isHidden;
        File[] f=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        Comparator<File> fileComparator=Comparator.comparing(File::getName);
        Arrays.sort(f, fileComparator);
        String[] outNames=Arrays.stream(f).map(File::getName).map(str->str.substring(0, 26))
                .map(str->str.replaceAll("_v_$",".txt")).distinct().toArray(String[]::new);
        int[] aa=IntStream.iterate(0, n->n+2).limit(45).toArray();
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetSize(aa.length, 7);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(index->{
                MAF.mergeTwoFiles(f[index], f[index+1], new File(mergeOutDir, outNames[index/2]));
            });
        }
    }

    public static void sort(String inputDir, String outDir){
        File[] input=IOUtils.listRecursiveFiles(new File(inputDir));
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).toArray(File[]::new);
        String[] outName=Arrays.stream(files).map(File::getName).toArray(String[]::new);
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetSize(files.length, 7);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(index->{
                MAF.sort(files[index], new File(outDir, outName[index]));
            });
        }
    }

    public static void sort(File inputFile, File outFile){
        long start=System.nanoTime();
        try(BufferedReader br=IOUtils.getTextReader(inputFile.getAbsolutePath());
            BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath())){
            String line, header;
            header=br.readLine();
            List<String> lines=new ArrayList<>();
            while ((line=br.readLine())!=null){
                lines.add(line);
            }
            Comparator<String> posComparator=Comparator.comparing(str->Integer.parseInt(PStringUtils.fastSplit(str).get(1)));
            Collections.sort(lines, posComparator);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < lines.size(); i++) {
                bw.write(lines.get(i));
                bw.newLine();
            }
            bw.flush();
            System.out.println(outFile.getName()+" is completed in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
        }catch (Exception e){
            e.printStackTrace();
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
        long start=System.nanoTime();
        File[] input=new File(inputDir).listFiles();
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).sorted().toArray(File[]::new);
        Map<Integer, BufferedReader> chrBufferReaderMap=new HashMap<>();
        Map<Integer, BufferedWriter> chrBufferWriterMap=new HashMap<>();
        IntStream.range(0, files.length).forEach(e->{
            chrBufferReaderMap.put(e, IOUtils.getNIOTextReader(files[e].getAbsolutePath()));
            chrBufferWriterMap.put(e, IOUtils.getNIOTextWriter(new File(outDir, files[e].getName()).getAbsolutePath()));
        });
        IntStream.range(0, files.length).parallel().forEach(index->{
            BufferedReader br;
            BufferedWriter bw;
            List<String> lines;
            String line;
            String header="CHR\tPOS\tAncestralAllele";
            List<String> lineList;
            try{
                br=chrBufferReaderMap.get(index);
                br.readLine();
                lines=new ArrayList<>();
                StringBuilder sb=new StringBuilder();
                while ((line=br.readLine())!=null){
                    lineList=PStringUtils.fastSplit(line);
                    if (!lineList.get(3).equals(lineList.get(4))) continue;
                    sb.append(lineList.get(0)).append("\t").append(lineList.get(1)).append("\t")
                            .append("\t").append(lineList.get(3));
                    lines.add(sb.toString());
                    sb=new StringBuilder();
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
            System.out.println(files[index].getName()+" is completed in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
        });
    }


}
