package daxing.filterSNP;

import format.position.ChrPos;
import utils.Benchmark;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.IntStream;

public class FilterSNPGo {

    /**
     * 根据从VCF得到的DepthDB, 以top density 0.75为标准过滤VCF，输出chr pos
     * @param depthsDBDir chr pos depth sd
     * @param outPutDir chr pos
     * @param topDensity 0.75
     */
    public static void getTopDensity(String depthsDBDir, String outPutDir, double topDensity){
        File[] input= IOUtils.listRecursiveFiles(new File(depthsDBDir));
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).sorted().toArray(File[]::new);
        String[] outFileName= Arrays.stream(files).map(File::getName).map(str->str.replaceAll("depth.txt.gz$", "chrpos.txt.gz")).sorted().toArray(String[]::new);
        IntStream.range(0, files.length).parallel().forEach(e->{
            Cells.getTopCellDensity(files[e], new File(outPutDir, outFileName[e]), topDensity);
        });
    }

    public static void getFilteredVCF(String vcfDir, String chrPosDir, String filteredVCFOutDir){
        File[] vcfFilesInput=IOUtils.listRecursiveFiles(new File(vcfDir));
        File[] chrPosFilesInput=IOUtils.listRecursiveFiles(new File(chrPosDir));
        Predicate<File> p=File::isHidden;
        File[] vcfFiles=Arrays.stream(vcfFilesInput).filter(p.negate()).sorted().toArray(File[]::new);
        File[] chrPosFiles=Arrays.stream(chrPosFilesInput).filter(p.negate()).sorted().toArray(File[]::new);
        if (vcfFiles.length!=chrPosFiles.length){
            System.out.println("please check vcfDir and chrPosDir, they must have equal number files");
            System.exit(1);
        }
        String[] filteredVCFOutName=Arrays.stream(chrPosFiles).map(File::getName)
                .map(str->str.replaceAll("chrpos.txt.gz$", "filtered0.75.vcf.gz")).sorted().toArray(String[]::new);
        IntStream.range(0, vcfFiles.length).parallel().forEach(e->{
            FilterSNPGo.getFilteredVCFfile(vcfFiles[e], chrPosFiles[e], new File(filteredVCFOutDir, filteredVCFOutName[e]));
        });
    }

    private static void getFilteredVCFfile(File vcfFile, File chrPosFile, File filteredVCFOutFile){
        try(BufferedReader br1=IOUtils.getTextReader(vcfFile.getAbsolutePath());
            BufferedReader br2=IOUtils.getTextGzipReader(chrPosFile.getAbsolutePath());
            BufferedWriter bw=IOUtils.getTextGzipWriter(filteredVCFOutFile.getAbsolutePath())){
            long start=System.nanoTime();
            ChrPos[] chrPosArray=br2.lines().skip(1).map(PStringUtils::fastSplit)
                    .map(e->new ChrPos(Short.parseShort(e.get(0)), Integer.parseInt(e.get(1)))).sorted().toArray(ChrPos[]::new);
            String line;
            List<String> lineList;
            short chr;
            int pos;
            ChrPos chrPos;
            int index;
            StringBuilder sb=new StringBuilder();
            while ((line=br1.readLine()).startsWith("##")){
                sb.append(line).append("\n");
            }
            sb.append(line).append("\n");
            bw.write(sb.toString());
            while ((line=br1.readLine())!=null){
                lineList=PStringUtils.fastSplit(line.substring(0, 20));
                chr= Short.parseShort(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                index=Arrays.binarySearch(chrPosArray, new ChrPos(chr, pos));
                if(index<0) continue;
                bw.write(line);
                bw.newLine();
            }
            System.out.println(filteredVCFOutFile.getName()+" is complicated in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

//    public static void main(String[] args) {
//        FilterSNPGo.getTopDensity(args[0], args[1], Double.parseDouble(args[2]));
//        FilterSNPGo.getFilteredVCF(args[0], args[1], args[2]);
//    }
}
