package daxing.applets;

import daxing.common.DateTime;
import daxing.common.PGF;
import daxing.common.RowTableTool;
import daxing.common.VCF;
import pgl.infra.dna.FastaByte;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.IntStream;

public class BadMutation {

    public static void getFaSubs(String geneHC_File, String pgf_File, String genomeFa_Dir, String chrGeneFaOutDir,
                                 String exon_VCFDir, String geneFaOutDir, String subsititonDir){
        long start=System.nanoTime();
        System.out.println(DateTime.getDateTimeOfNow()+" start");
        RowTableTool<String> hcGeneTable=new RowTableTool<>(geneHC_File);
        Set<String> hcGeneSet =new HashSet<>(hcGeneTable.getColumn(0));
        Predicate<PGF.Gene> hcGenes= gene -> hcGeneSet.contains(gene.getGeneName());
        PGF pgf=new PGF(pgf_File);
        pgf.removeIf(hcGenes.negate());
        pgf.writeCDSSequencePerChr(genomeFa_Dir, chrGeneFaOutDir);
        writeSubstitionFile(pgf,chrGeneFaOutDir, exon_VCFDir, geneFaOutDir, subsititonDir);
        System.out.println("completed in "+ Benchmark.getTimeSpanHours(start)+" hours");
        System.out.println(DateTime.getDateTimeOfNow()+" end");
    }

    private static void writeSubstitionFile(PGF pgf, String chrGeneFaDir, String exon_VCFDir, String geneFaOutDir,
                                            String subsititonDir){
        long start=System.nanoTime();
        Predicate<File> predicate=File::isHidden;
        File[] chGeneFiles= IOUtils.listRecursiveFiles(new File(chrGeneFaDir));
        File[] f=Arrays.stream(chGeneFiles).filter(predicate.negate()).toArray(File[]::new);
        List<String>[] vcfPosOnAllChr= VCF.getColumnList(exon_VCFDir, 1);
        PGF[] chrPgf=pgf.getGeneOnAllChr();
        IntStream.range(0, f.length).parallel().forEach(e->
                writeSubstitionFile(chrPgf[e], f[e], vcfPosOnAllChr[e], new File(geneFaOutDir), new File(subsititonDir)));
        System.out.println("all substitution file were write to "+subsititonDir+" in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
    }

    private static void writeSubstitionFile(PGF chrPgf, File chrGeneFaFile, List<String> vcfPosOnChr, File geneFaOutDir,
                                            File subsititonOutDir){
        int[] vcfPosOnChrArray=vcfPosOnChr.stream().mapToInt(Integer::parseInt).toArray();
        FastaByte chrFa=new FastaByte(chrGeneFaFile.getAbsolutePath());
        String[] transcriptName=chrFa.getNames();
        TIntArrayList geneCDSPosList;
        int[] geneCDSPosArray;
        if (chrPgf.getGeneNumber()!=chrFa.getSeqNumber()){
            System.out.println("error, check pgf and chrGeneFaFile");
            System.exit(1);
        }
        BufferedWriter bwSubs;
        StringBuilder sb;
        int chr=chrPgf.getChrs().get(0);
        try {
            for (int i = 0; i < chrPgf.getGeneNumber(); i++) {
                int longestTrancriptIndex=chrPgf.getLongestTranscriptIndex(i);
                geneCDSPosList=chrPgf.getCDSPosOnGene(i, longestTrancriptIndex);
                geneCDSPosArray=geneCDSPosList.toArray();
                geneCDSPosList.retainAll(vcfPosOnChrArray);
                if (geneCDSPosList.size()==0) continue;
//                if (chrPgf.getTranscriptStrand(i, longestTrancriptIndex)==0){
//                    System.out.println(chrPgf.getTranscriptName(i, longestTrancriptIndex));
//                }
                chrFa.writeFasta(new File(geneFaOutDir, transcriptName[i]+".fasta").getAbsolutePath(), i, IOFileFormat.Text);
                geneCDSPosList.sort();
                bwSubs=IOUtils.getTextWriter(new File(subsititonOutDir, transcriptName[i]+".subs").getAbsolutePath());
                sb=new StringBuilder();
                for (int j = 0; j < geneCDSPosList.size(); j++) {
                    if (chrPgf.getTranscriptStrand(i, longestTrancriptIndex)==0){
                        int index= Arrays.binarySearch(geneCDSPosArray, geneCDSPosList.get(j));
                        index=((geneCDSPosArray.length)/3+1)-(index/3+1);
                        sb.append(index).append("\t").append(chr).append("_").append(geneCDSPosList.get(j));
                    }else {
                        int index=Arrays.binarySearch(geneCDSPosArray, geneCDSPosList.get(j));
                        sb.append(index/3+1).append("\t").append(chr).append("_").append(geneCDSPosList.get(j));
                    }
                    bwSubs.write(sb.toString());
                    bwSubs.newLine();
                    sb=new StringBuilder();
                }
                bwSubs.flush();
                bwSubs.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

//    public static void main(String[] args) {
//        BadMutation.getFaSubs("", "", "", "", "", "", "");
//    }

}
