package xuebo.analysis.wheat_epi;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.set.hash.HashIntSet;
import com.koloboke.collect.set.hash.HashIntSets;
import com.koloboke.collect.set.hash.HashLongSet;
import com.koloboke.collect.set.hash.HashLongSets;
import java.io.DataOutputStream;
import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.FastaByte;

import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;

public class ReferenceKmerLib {
    int kmerLength = -1;
    HashIntSet intSet = null;
    HashLongSet longSet = null;
    
    public ReferenceKmerLib(int kmerLength, String inputGenomeFileS) {
        this.kmerLength = kmerLength;
        this.createKmerSet(inputGenomeFileS);
    }

    void createKmerSet (String referenceGenomeFileS) {
//        referenceGenomeFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/maize_chr12.fa";
//        inputGenomeFileS = "/Users/feilu/Documents/database/maize/reference/AGPv4/maizeAGPv4.fa";
//        inputGenomeFileS = "/Users/feilu/Documents/database/maize/reference/download/Zea_mays.AGPv4.dna.chromosome.10.fa.gz";
        FastaByte f = new FastaByte(referenceGenomeFileS);
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIBaseByteMap();
        System.out.println("Building kmer library from reference...");
        System.out.println("KmerLength = "+String.valueOf(kmerLength)+ " bp");
        long start = System.nanoTime();
        if (kmerLength == 16) {
            this.intSet = this.getIntKmerSet(f, ascIIByteMap);
        }
        else if (kmerLength == 32) {
            this.longSet = this.getLongKmerSet(f, ascIIByteMap);
        }
        System.out.println(Benchmark.getTimeSpanSeconds(start)+" seconds used to build Kmer library");
    }
    
    void writeBinaryFile (String libFileS) {
//        libFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/kmerLib/kmerLib.bin";
        try {
            DataOutputStream dos  = IOUtils.getBinaryWriter(libFileS);
            dos.writeInt(kmerLength);
            if (kmerLength == 16) {
                int[] kmerD = intSet.toIntArray();
                dos.writeInt(kmerD.length);
                for (int i = 0; i < kmerD.length; i++) {
                    dos.writeInt(kmerD[i]);
                }
            }
            else if (kmerLength == 32) {
                long[] kmerD = longSet.toLongArray();
                dos.writeInt(kmerD.length);
                for (int i = 0; i < kmerD.length; i++) {
                    dos.writeLong(kmerD[i]);
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Kmer library is written to " + libFileS);
    }
    
    private HashLongSet getLongKmerSet (FastaByte f, HashByteByteMap ascIIByteMap) {
        int genomeSize = (int)f.getTotalSeqLength();
        HashLongSet kmerSet = HashLongSets.newMutableSet(genomeSize);
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            int mark = 0;
            boolean flag = false;
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                flag = false;
                for (int j = mark; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) {
                    mark = i + 1;
                    continue;
                }
                else {
                    mark = i + kmerLength;
                }
                long kmerL = BaseEncoder.getLongSeqFromSubByteArray(bArray, i, i + kmerLength);
                if (!kmerSet.contains(kmerL)) kmerSet.add(kmerL);
                int pos = i+1;
                if (pos%50000000 == 0) {
                    System.out.println("Chromosome: "+f.getName(k)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos) + ". Kmer set size: " + String.valueOf(kmerSet.size()));
                }
            }
        }
        return kmerSet;
    }
    
    private HashIntSet getIntKmerSet (FastaByte f, HashByteByteMap ascIIByteMap) {
        int genomeSize = (int)f.getTotalSeqLength();
        HashIntSet kmerSet = HashIntSets.newMutableSet(genomeSize);
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            int mark = 0;
            boolean flag = false;
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                flag = false;
                for (int j = mark; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) {
                    mark = i + 1;
                    continue;
                }
                else {
                    mark = i + kmerLength;
                }
                int kmerL = BaseEncoder.getIntSeqFromSubByteArray(bArray, i, i + kmerLength);
                if (!kmerSet.contains(kmerL)) kmerSet.add(kmerL);
                int pos = i+1;
                if (pos%50000000 == 0) {
                    System.out.println("Chromosome: "+f.getName(k)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos) + ". Kmer set size: " + String.valueOf(kmerSet.size()));
                }
            }
        }
        return kmerSet;
    }
}

