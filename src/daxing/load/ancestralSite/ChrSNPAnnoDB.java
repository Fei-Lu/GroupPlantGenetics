package daxing.load.ancestralSite;

import daxing.common.chrrange.ChrRange;
import daxing.common.chrrange.ChrRanges;
import daxing.common.utiles.IOTool;
import daxing.common.table.RowTableTool;
import pgl.infra.pos.ChrPos;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.*;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ChrSNPAnnoDB {

    private List<SNPAnnotation> snpAnnotationList;

    public ChrSNPAnnoDB(File exonSNPAnnoFile){
        this.snpAnnotationList =this.getGeneSNPAnno(exonSNPAnnoFile);
        Collections.sort(snpAnnotationList);
    }

    private List<SNPAnnotation> getGeneSNPAnno(File exonSNPAnnoFile){
        List<SNPAnnotation> geneAnno=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(exonSNPAnnoFile)) {
            br.readLine();
            String line;
            while ((line=br.readLine())!=null){
                geneAnno.add(getSNPAnnotation2(line));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return geneAnno;
    }

    /**
     * for vmap2 643
     * @param line
     * @return
     */
    private SNPAnnotation getSNPAnnotation(String line){
        List<String> temp= PStringUtils.fastSplit(line);
        short chr=Short.parseShort(temp.get(1));
        int pos=Integer.parseInt(temp.get(2));
        char refBase=temp.get(3).charAt(0);
        char altBase=temp.get(4).charAt(0);
        char majorBase=temp.get(5).charAt(0);
        double[] aaf=new double[2];
        double maf=Double.parseDouble(temp.get(7));
        aaf[0]=Double.parseDouble(temp.get(8));
        aaf[1]=Double.parseDouble(temp.get(9));
        String geneName=temp.get(10).substring(0,18);
        SNPAnnotation.Region region=SNPAnnotation.Region.valueOf(temp.get(11));
        String variant_type=temp.get(12);
        String ancestral=temp.get(15);
        String derived_SIFT=temp.get(16);
        String daf=temp.get(17);
        String[] dafs=new String[2];
        dafs[0]=temp.get(18);
        dafs[1]=temp.get(19);
        String gerp=temp.get(20);
        String recombinationRate=null;
        return new SNPAnnotation(chr, pos, refBase, altBase, geneName, majorBase, ancestral, maf
                , aaf, daf, dafs, region, variant_type, derived_SIFT, gerp, recombinationRate,
                null, null, null, null, null,
                null, null, null,null,null,null,
                null, null, null, null);
    }

    /**
     * for vmap2 1061
     * @param line
     * @return
     */
    private SNPAnnotation getSNPAnnotation2(String line){
        List<String> temp= PStringUtils.fastSplit(line);
        short chr=Short.parseShort(temp.get(1));
        int pos=Integer.parseInt(temp.get(2));
        char refBase=temp.get(3).charAt(0);
        char altBase=temp.get(4).charAt(0);
        char majorBase=temp.get(5).charAt(0);
        double maf=Double.parseDouble(temp.get(7));
//        aaf[0]=Double.parseDouble(temp.get(8));
//        aaf[1]=Double.parseDouble(temp.get(9));
        String geneName=temp.get(8).substring(0,18);
        String ancestral=temp.get(9);
        String daf=temp.get(10);
        String gerp=temp.get(11);
        SNPAnnotation.Region region=SNPAnnotation.Region.valueOf(temp.get(12));
        String variant_type=temp.get(13);
        String derived_SIFT=temp.get(16);
        String effect_vep=temp.get(17);
        String impact_vep=temp.get(18);
        String effect_snpEff=temp.get(19);
        String impact_snpEff=temp.get(20);
        String gerp16way=temp.get(22);
        String aaPos=temp.get(23);
        String listS2=temp.get(25);
        String phylop=temp.get(26);
        String aaSize=temp.get(27);
        String phylopRefMask=temp.get(28);
        String alleleAgeJ=temp.get(29);
        String alleleAgeM=temp.get(30);
        String alleleAgeR=temp.get(31);

        String recombinationRate=null;
        double[] aaf = null;
        String[] dafs = null;

        String derived_PolyPhen2_prediction=temp.get(38);
        String derived_PolyPhen2_class=temp.get(39);
        String derived_PolyPhen2_prob=temp.get(40);

        return new SNPAnnotation(chr, pos, refBase, altBase, geneName, majorBase, ancestral, maf
                , aaf, daf, dafs, region, variant_type, derived_SIFT, gerp, recombinationRate,
                effect_snpEff, impact_snpEff, effect_vep, impact_vep,gerp16way,aaPos, aaSize, listS2, phylopRefMask,
                alleleAgeJ,alleleAgeM, alleleAgeR, derived_PolyPhen2_prediction, derived_PolyPhen2_class, derived_PolyPhen2_prob);
    }

    public int getSize(){
        return this.snpAnnotationList.size();
    }

    public List<SNPAnnotation> getSnpAnnotationList() {
        return snpAnnotationList;
    }

    /**
     *
     * @param chr chr
     * @param pos pos
     * @return negative value, not in DB
     */
    public int binarySearch(int chr, int pos){
        return Collections.binarySearch(this.snpAnnotationList, new ChrPos((short) chr, pos));
    }

    public SNPAnnotation getSNPAnnotation(int index){
        return this.snpAnnotationList.get(index);
    }

    public String getGeneName(int chr, int pos){
        int snpIndex= this.binarySearch(chr, pos);
        if (snpIndex < 0){
            System.out.println("error");
            System.exit(1);
        }
        return this.getSNPAnnotation(snpIndex).getSNPInfo();
    }

    public String[] getGeneName(){
        List<SNPAnnotation> snpAnnotations=this.snpAnnotationList;
        Set<String> geneNames=new HashSet<>();
        String geneName;
        for (SNPAnnotation snpAnnotation : snpAnnotations) {
            geneName = snpAnnotation.getSNPInfo();
            geneNames.add(geneName);
        }
        return geneNames.stream().sorted().toArray(String[]::new);
    }

    public SNPAnnotation getSNP(int chr, int pos){
        int index=this.binarySearch(chr, pos);
        return this.snpAnnotationList.get(index);
    }

    public boolean contain(int chr, int pos){
        int index=this.binarySearch(chr, pos);
        return index >= 0;
    }

    /**
     *
     * @param chr chr
     * @param pos pos
     * @return variantType
     */
    public String getVariantType(int chr, int pos){
        int index=this.binarySearch(chr, pos);
        if (index < 0) return "INTERGENIC";
        return this.getSNP(chr, pos).variant_type;
    }

    public boolean hasAncestral(int chr, int pos){
        return this.getSNP(chr, pos).hasAncestral();
    }

    public boolean isRefAlleleAncestral(int chr, int pos){
        return this.getSNP(chr, pos).isRefAlleleAncestral();
    }

    public boolean isSyn(int chr, int pos){
        return this.getSNP(chr, pos).isSyn();
    }

    public boolean isNonsyn(int chr, int pos){
        return this.getSNP(chr, pos).isNonSyn();
    }

    /**
     * for SIFT && GERP
     * @param chr
     * @param pos
     * @return
     */
    public boolean isDeleterious(int chr, int pos){
        return this.getSNP(chr, pos).isDeleterious();
    }

    public boolean isDeleterious(int chr, int pos, SNPAnnotation.MethodCallDeleterious methodCallDeleterious){
        return this.getSNP(chr, pos).isDeleterious(methodCallDeleterious);
    }

    public int getDelNum(){
        List<SNPAnnotation> snpAnnotations=this.snpAnnotationList;
        Predicate<SNPAnnotation> p = SNPAnnotation::isDeleterious;
        return (int) snpAnnotations.stream().filter(p).count();
    }

    public int getNum(Predicate<SNPAnnotation> p){
        return (int) this.snpAnnotationList.stream().filter(p).count();
    }

    public ChrSNPAnnoDB filter(Predicate<SNPAnnotation> p){
        this.snpAnnotationList= this.snpAnnotationList.stream().filter(p).collect(Collectors.toList());
        return this;
    }

    public void write(String outFile){
        StringBuilder sb=new StringBuilder();
        SNPAnnotation snpAnnotation;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Chr\tPos");
            bw.newLine();
            for (SNPAnnotation annotation : snpAnnotationList) {
                snpAnnotation = annotation;
                sb.setLength(0);
                sb.append(snpAnnotation.getChromosome()).append("\t").append(snpAnnotation.getPos());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * split DB by chr
     * @param snpAnnpDBFile
     * @param outDir
     */
    public static void splitSNPAnnoDB(String snpAnnpDBFile, String outDir){
        Set<String> set= RowTableTool.getColumnSet(snpAnnpDBFile, 1);
        List<String> chrList= new ArrayList<>(set);
        Collections.sort(chrList);
        BufferedWriter[] bws = new BufferedWriter[chrList.size()];
        String file;
        for (int i = 0; i < chrList.size(); i++) {
            file = "chr"+PStringUtils.getNDigitNumber(3, Integer.parseInt(chrList.get(i)))+"_"+new File(snpAnnpDBFile).getName();
            bws[i] = IOTool.getWriter(new File(outDir, file));
        }
        try (BufferedReader br = IOTool.getReader(snpAnnpDBFile)) {
            String line, header;
            List<String> temp;
            header = br.readLine();
            for (int i = 0; i < bws.length; i++) {
                bws[i].write(header);
                bws[i].newLine();
            }
            while ((line=br.readLine())!=null){
                temp =PStringUtils.fastSplit(line);
                int index = Collections.binarySearch(chrList, temp.get(1));
                bws[index].write(line);
                bws[index].newLine();
            }
            for (int i = 0; i < bws.length; i++) {
                bws[i].flush();
                bws[i].close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void filterRange(String rangeFile, String snpAnnpDBDir, String outDir){
        ChrRanges chrRanges=getChrRange(rangeFile);
        chrRanges.sortBy(ChrRanges.SortType.POSITION);
        List<File> snpAnnDBFiles=IOTool.getFileListInDirEndsWith(snpAnnpDBDir, ".gz");
        String[] outNames=snpAnnDBFiles.stream().map(File::getName).toArray(String[]::new);
        IntStream.range(0, snpAnnDBFiles.size()).forEach(e->{
            try (BufferedReader br = IOTool.getReader(snpAnnDBFiles.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                List<String> temp;
                String line, refChr;
                line=br.readLine();
                bw.write(line);
                bw.newLine();
                int chrID, pos, refPos;
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chrID = Integer.parseInt(temp.get(1));
                    pos = Integer.parseInt(temp.get(2));
                    refChr= RefV1Utils.getChromosome(chrID, pos);
                    refPos = RefV1Utils.getPosOnChromosome(chrID, pos);
                    if (chrRanges.contain(refChr, refPos)) continue;
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    private static ChrRanges getChrRange(String rangeFile){
        List<ChrRange> rangeList=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(rangeFile)) {
            String line, refChr;
            int refStart, refEnd;
            ChrRange chrRange;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                refChr=temp.get(0);
                refStart=Integer.parseInt(temp.get(1));
                refEnd=Integer.parseInt(temp.get(2))+1;
                chrRange = new ChrRange(refChr, refStart, refEnd);
                rangeList.add(chrRange);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new ChrRanges(rangeList);
    }
}
