package daxing.individualIntrogression;

import daxing.common.factors.LoadType;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.DateTime;
import daxing.common.utiles.IOTool;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import daxing.load.ancestralSite.SNPAnnotation;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.IntStream;

/**
 * 以群体fd window为区间, 计算introgressed个体和nonintrogressed个体的 total del count
 */
public class PopulationIndividualFdStart {

    public static void start(String exonSNPAnnoDir, String exonVCFDir, String taxa_InfoDBFile,
                             SNPAnnotation.MethodCallDeleterious methodCallDeleterious, String popFdFile,
                             String individualFdDir,
                             String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        String[] subDir= {"001_siteGridPerChr","002_siteGridAllChr", "003_introgressionDonorBurden"};
        File[] subFiles = new File[subDir.length];
        for (int i = 0; i < subFiles.length; i++) {
            subFiles[i] = new File(outDir, subDir[i]);
            subFiles[i].mkdir();
        }
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        String[] outNames = exonVCFFiles.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".siteGridLoad" +
                ".txt.gz")).toArray(String[]::new);
        IntStream.range(0, exonVCFFiles.size()).parallel().forEach(e-> getSiteGridLoad(exonAnnoFiles.get(e), exonVCFFiles.get(e),
                taxa_InfoDBFile, e+1, methodCallDeleterious, new File(subFiles[0], outNames[e])));
        Path siteGridPath = Paths.get(new File(subFiles[1], "siteGrid.txt.gz").getAbsolutePath()).toAbsolutePath();
        merge(subFiles[0], siteGridPath.toFile());
        Path introgressionDonorBurdenPath = Paths.get(new File(subFiles[2], "introgressionDonorBurden.txt.gz").getAbsolutePath());
        PopulationIndividualFd.writeWindowSize(popFdFile, individualFdDir, taxa_InfoDBFile,siteGridPath.toString(),
                introgressionDonorBurdenPath.toString());
        System.out.println(DateTime.getDateTimeOfNow());
    }

    private static void getSiteGridLoad(File exonSNPAnnoFile, File exonVCFFile,
                              String taxa_InfoDBFile, int chr, SNPAnnotation.MethodCallDeleterious methodCallDeleterious, File outFile){
        try (BufferedReader br = IOTool.getReader(exonVCFFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){}
            temp= PStringUtils.fastSplit(line);
            List<String> taxonNames=temp.subList(9, temp.size());
            ChrSNPAnnoDB transcriptDB=new ChrSNPAnnoDB(exonSNPAnnoFile);
            IndividualChrPosLoad[] taxonLoads=new IndividualChrPosLoad[taxonNames.size()];
            for (int i = 0; i < taxonLoads.length; i++) {
                taxonLoads[i]=new IndividualChrPosLoad(taxonNames.get(i), chr);
            }
            int pos, depth;
            String genotype;
            List<String> genotypeList, depthList;
            boolean isRefAlleleAncestral;
            boolean isSyn, isNonsyn, isDeleterious;
            byte genotypeByte;
            LoadType loadType;
            boolean ifHeter, ifHomozygousDerived;
            Optional<LoadType> optionalLoadType;
            TIntArrayList posList = new TIntArrayList();
            List<LoadType> loadTypeList = new ArrayList<>();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                pos=Integer.parseInt(temp.get(1));
                if (!transcriptDB.contain(chr, pos)) continue;
                if (!transcriptDB.hasAncestral(chr, pos)) continue;
                isSyn=transcriptDB.isSyn(chr, pos);
                isNonsyn=transcriptDB.isNonsyn(chr, pos);
                if (!(isSyn || isNonsyn)) continue;
                isDeleterious=transcriptDB.isDeleterious(chr, pos, methodCallDeleterious);
                isRefAlleleAncestral=transcriptDB.isRefAlleleAncestral(chr, pos);
                optionalLoadType=LoadType.getInstanceFromString(transcriptDB.getVariantType(chr, pos));
                loadType=optionalLoadType.orElseThrow(IllegalArgumentException::new);
                loadType= isDeleterious ? LoadType.Del : loadType;
                posList.add(pos);
                loadTypeList.add(loadType);
                for (int i = 0; i < taxonNames.size(); i++) {
                    genotypeList=PStringUtils.fastSplit(temp.get(i+9), ":");
                    genotype=genotypeList.get(0);
                    if (genotype.equals("./.")) continue;
//                    if (genotype.equals("0/1")) continue;
                    depthList=PStringUtils.fastSplit(genotypeList.get(1),",");
                    depth=Integer.parseInt(depthList.get(0))+Integer.parseInt(depthList.get(1));
                    if (depth < 2) continue;
                    genotypeByte= IndividualChrPosLoad.caculateGenotype(genotype, isRefAlleleAncestral);
                    ifHeter = genotypeByte==2 ? true : false;
                    ifHomozygousDerived = genotypeByte==1 ? true : false;
                    taxonLoads[i].addChrPos(chr, pos, loadType, ifHeter, ifHomozygousDerived);
                }
            }
            RowTableTool<String> table = new RowTableTool<>(taxa_InfoDBFile);
            Predicate<List<String>> lrclP = l->l.get(10).equals("Landrace") || l.get(10).equals("Cultivar");
            table.removeIf(lrclP.negate());
            List<String> taxaLRCL= table.getColumn(0);
            Map<String, String> taxaIntrogressIDMap=table.getHashMap(0, 35);
            String[] taxaLRCL_ID = taxaLRCL.stream().map(taxaIntrogressIDMap::get).toArray(String[]::new);
            IndividualChrPosLoad[] taxaLoad = Arrays.stream(taxonLoads).parallel().filter(grid -> taxaLRCL.contains(grid.getTaxonName())).toArray(IndividualChrPosLoad[]::new);
            assert taxaLoad.length==taxaLRCL.size() : "check taxon name";
            int siteIndex, taxonIndex;
            StringBuilder sb = new StringBuilder();
            sb.append("Chr\tPos\tLoadType\t").append(String.join("\t", taxaLRCL_ID));
            bw.write(sb.toString());
            bw.newLine();
            Arrays.sort(taxaLoad);
            for (int i = 0; i < posList.size(); i++) {
                sb.setLength(0);
                sb.append(chr).append("\t").append(posList.get(i)).append("\t");
                sb.append(loadTypeList.get(i)).append("\t");
                for (String taxon : taxaLRCL) {
                    taxonIndex = Arrays.binarySearch(taxaLoad, new IndividualChrPosLoad(taxon, -1));
                    assert taxonIndex >=0 : "check taxon name";
                    siteIndex = taxaLoad[taxonIndex].getIndex(posList.get(i));
                    if (siteIndex < 0){
                        sb.append("-1,-1").append("\t");
                    }else {
                        sb.append(taxaLoad[taxonIndex].getIfHeter(siteIndex)).append(",");
                        sb.append(taxaLoad[taxonIndex].getIfHomozygousDerived(siteIndex)).append("\t");
                    }
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void merge(File siteGridPerChrDir, File outFile){
        List<File> siteGridFiles = IOTool.getFileListInDirEndsWith(siteGridPerChrDir.getAbsolutePath(), ".txt.gz");
        try {
            BufferedReader br=IOTool.getReader(siteGridFiles.get(0));
            BufferedWriter bw= IOTool.getWriter(outFile);
            String line=br.readLine();
            bw.write(line);
            bw.newLine();
            while ((line=br.readLine())!=null){
                bw.write(line);
                bw.newLine();
            }
            br.close();
            for (int i = 1; i < siteGridFiles.size(); i++) {
                br = IOTool.getReader(siteGridFiles.get(i));
                br.readLine();
                while ((line=br.readLine())!=null){
                    bw.write(line);
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
