package daxing.ancestralAllele;

import daxing.common.ChrConvertionRule;
import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;
import java.io.BufferedReader;
import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Predicate;
import java.util.stream.Collectors;


public class MAF {

    private List<MAFrecord>[] mafRecords;
    private int numThreads=12;
    private String taxonForOrder;
    private ChrConvertionRule chrConvertionRule;

    /**
     *
     * @param taxonForOrder must be in format like "hordeum_vulgare" or "triticum_aestivum"
     * @param chrConvertionRule
     * @param mafInputFileDir
     */
    public MAF(String taxonForOrder, ChrConvertionRule chrConvertionRule, Path mafInputFileDir){
        this.taxonForOrder=taxonForOrder;
        this.chrConvertionRule=chrConvertionRule;
        this.initialize(mafInputFileDir);
    }

    private void initialize(Path mafInputFileDir){
        File[] files=IOUtils.listRecursiveFiles(new File(mafInputFileDir.toString()));
        Predicate<File> p=e->e.getName().contains(".DS_Store");
        File[] fileArray= Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
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
            integerList.stream()
                    .forEach(index-> {
                        try(BufferedReader br= IOUtils.getNIOTextReader(fileArray[index].toString())){
                            br.readLine();
                            br.readLine();
                            long id=Long.MIN_VALUE;
                            int score=Integer.MIN_VALUE;
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
                                score=Integer.parseInt(lineList2.get(1));
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
        String[] taxons=mafRecordList.get(0).getTaxon();
        int indexOfTaxon=Arrays.binarySearch(taxons, this.getTaxonForOrder());
        int chrID;
        String chr;
        int startPos;
        for(MAFrecord e: mafRecordList){
            chr=e.getChr(indexOfTaxon);
            startPos=e.getStartPos(indexOfTaxon);
            chrID=this.getChrConvertionRule().getChrIDFromOriChrName(chr, startPos);
            mafrecord[chrID].add(e);
        }
        Comparator<MAFrecord> comparator=Comparator.comparing(m->m.getStartPos(indexOfTaxon));
        Arrays.stream(mafrecord).forEach(e-> Collections.sort(e, comparator));
        this.mafRecords=mafrecord;
    }

    public List<MAFrecord>[] getMafRecords() {
        return mafRecords;
    }

    public int getNumThreads() {
        return numThreads;
    }

    public String getTaxonForOrder(){
        return taxonForOrder;
    }

    public ChrConvertionRule getChrConvertionRule(){
        return chrConvertionRule;
    }

    public void setNumOfThreads(int numThreads){
        this.numThreads=numThreads;
    }


}
