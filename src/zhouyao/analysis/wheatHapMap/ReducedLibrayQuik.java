///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//package zhouyao.analysis.wheatHapMap;
//
//import java.io.File;
//import java.util.Arrays;
//import java.util.List;
//import java.util.concurrent.atomic.AtomicInteger;
//
///**
// *
// * @author yaozhou
// */
//public class ReducedLibrayQuik {
//    ReducedLibrayQuik(String inFile,String cutter1,String cutter2, String outFile) {
//
//        this.getSubString(inFile,cutter1,cutter2,outFile);
////        this.outputFiles(outfileS1, outfileS2);
//    }
//    String infileDirS = "/Users/yaozhou/Documents/RefGenome/test/ref" ;
//    String abdFileS = "/Users/yaozhou/Documents/RefGenome/test/ABD/abd_iwgscV1.fa.gz";
//    String aFileS = "/Users/yaozhou/Documents/RefGenome/test/A/a_iwgscV1.fa.gz";
//    String bFileS = "/Users/yaozhou/Documents/RefGenome/test/B/b_iwgscV1.fa.gz";
//    String abFileS = "/Users/yaozhou/Documents/RefGenome/test/AB/ab_iwgscV1.fa.gz";
//    String dFileS = "/Users/yaozhou/Documents/RefGenome/test/D/d_iwgscV1.fa.gz";
//    public void getSubString(String inFile,String cutter1,String cutter2, String outFile){
//        File[] fArray = new File(infileDirS).listFiles();
//        fArray = IOUtils.listFilesEndsWith(fArray, ".gz");
//        List<File> fList = Arrays.asList(fArray);
//        FastaBit[] faArray = new FastaBit[fList.size()];
//        AtomicInteger atoCnt = new AtomicInteger();
//        fList.stream().forEach(f -> {
//            try {
//                FastaBit fb = new FastaBit(f.getAbsolutePath());
//                faArray[atoCnt.getAndIncrement()] = fb;
//            }
//            catch (Exception e) {
//                e.printStackTrace();
//            }
//        });
//        FastaBit afb = new FastaBit (faArray);
//        afb.sortByName();
//        afb.writeFasta(abdFileS, IOFileFormat.TextGzip);
//        boolean[] ifOut = new boolean[afb.getSeqNumber()];
//        for (int i = 0; i < afb.getSeqNumber(); i++) {
//            if (afb.getName(i).split("_")[0].endsWith("A")) ifOut[i] = true;
//            if (afb.getName(i).contains("Mit")) ifOut[i] = true;
//            if (afb.getName(i).contains("Chl")) ifOut[i] = true;
//            if (afb.getName(i).contains("Un")) ifOut[i] = true;
//        }
//        afb.writeFasta(aFileS, ifOut, IOFileFormat.TextGzip);
//        ifOut = new boolean[afb.getSeqNumber()];
//        for (int i = 0; i < afb.getSeqNumber(); i++) {
//            if (afb.getName(i).split("_")[0].endsWith("A") || afb.getName(i).split("_")[0].endsWith("B")) ifOut[i] = true;
//            if (afb.getName(i).contains("Mit")) ifOut[i] = true;
//            if (afb.getName(i).contains("Chl")) ifOut[i] = true;
//            if (afb.getName(i).contains("Un")) ifOut[i] = true;
//        }
//        afb.writeFasta(abFileS, ifOut, IOFileFormat.TextGzip);
//        ifOut = new boolean[afb.getSeqNumber()];
//        for (int i = 0; i < afb.getSeqNumber(); i++) {
//            if (afb.getName(i).split("_")[0].endsWith("D")) ifOut[i] = true;
//            if (afb.getName(i).contains("Mit")) ifOut[i] = true;
//            if (afb.getName(i).contains("Chl")) ifOut[i] = true;
//            if (afb.getName(i).contains("Un")) ifOut[i] = true;
//        }
//        afb.writeFasta(dFileS, ifOut, IOFileFormat.TextGzip);
//        ifOut = new boolean[afb.getSeqNumber()];
//        for (int i = 0; i < afb.getSeqNumber(); i++) {
//            if (afb.getName(i).split("_")[0].endsWith("B")) ifOut[i] = true;
//            if (afb.getName(i).contains("Mit")) ifOut[i] = true;
//            if (afb.getName(i).contains("Chl")) ifOut[i] = true;
//            if (afb.getName(i).contains("Un")) ifOut[i] = true;
//        }
//        afb.writeFasta(bFileS, ifOut, IOFileFormat.TextGzip);
//    }
//}
