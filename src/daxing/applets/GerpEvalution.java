package daxing.applets;

import daxing.common.ChrConvertionRule;
import daxing.common.NumberTool;
import daxing.common.PGF;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;

/**
 * chr001 chr002 chr003, ...
 */
public class GerpEvalution {

    private TIntArrayList[] pos;
    private TDoubleArrayList[] rsScore;

    /**
     * Chr     Pos     GerpNeutralRate GerpScore
     * 1       1       0       0
     * 1       2       0       0
     * 1       3       0       0
     * 1       4       0       0
     * 1       5       0       0
     * @param gerpScoresDir
     */
    public GerpEvalution(String gerpScoresDir){
        this.initilize(gerpScoresDir);
    }

    private void initilize(String gerpScoresDir){
        File[] files= IOUtils.listRecursiveFiles(new File(gerpScoresDir));
        Predicate<File> p=File::isHidden;
        File[] f= Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        int[] chrs= Arrays.stream(f).map(File::getName).map(str->Integer.parseInt(str.substring(3,6)))
                                    .mapToInt(Integer::intValue).toArray();
        pos=new TIntArrayList[f.length];
        rsScore=new TDoubleArrayList[f.length];
        for (int i = 0; i < f.length; i++) {
            this.readGerpFile(f[i], chrs[i]);
        }
    }

    private void readGerpFile(File gerpFile, int chr){
        TIntArrayList posList=new TIntArrayList();
        TDoubleArrayList gerpScoreList=new TDoubleArrayList();
        try (BufferedReader br = IOUtils.getTextReader(gerpFile.getAbsolutePath())) {
            String line;
            List<String> lineList;
            int chrLine=-1;
            int posLine=-1;
            double gerpScore=-1d;
            br.readLine();
            while ((line=br.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                chrLine=Integer.parseInt(lineList.get(0));
                posLine=Integer.parseInt(lineList.get(1));
                gerpScore=Double.parseDouble(lineList.get(3));
                if (chrLine!=chr){
                    System.out.println("file error, check your "+gerpFile.getName());
                    System.exit(1);
                }
//                if (gerpScore==0) continue;
                posList.add(posLine);
                gerpScoreList.add(gerpScore);
            }
            pos[chr-1]=posList;
            rsScore[chr-1]=gerpScoreList;
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public TDoubleArrayList[] getRsScore() {
        return rsScore;
    }

    public TIntArrayList[] getPos() {
        return pos;
    }

    public int getTotalSites(){
        int count=0;
        for (int i = 0; i < this.getRsScore().length; i++) {
            count +=this.getRsScore()[i].size();
        }
        return count;
    }

    public int evaluateGerpScore(){
        int greatThan0=0;
        int equal0=0;
        int lessThan0=0;
        double total=this.getTotalSites();
        for (int i = 0; i < rsScore.length; i++) {
            for (int j = 0; j < rsScore[i].size(); j++) {
                if (rsScore[i].get(j)>0){
                    greatThan0++;
                }else if (rsScore[i].get(j)==0){
                    equal0++;
                }else {
                    lessThan0++;
                }

            }
        }
        System.out.println("Total sites: "+(int)total+"\t"+"Constrained regions: "+greatThan0+"\n"+ "Less than 0: "+lessThan0+
                "\n"+"equal 0: "+equal0);
        return greatThan0;
    }

    public double ratio(double start, double end){
        int count=0;
        double total=this.getTotalSites();
        for (int i = 0; i < this.getRsScore().length; i++) {
            for (int j = 0; j < this.getRsScore()[i].size(); j++) {
                if (this.getRsScore()[i].get(j)<start || this.getRsScore()[i].get(j)>=end) continue;
                count++;
            }
        }
        System.out.println("Total sites: "+NumberTool.parse((int)total)+"\t"+"["+start+", "+end+") "+NumberTool.parse(count)+
                "\t"+ NumberTool.format( 100*count/total,
                4)+"%");
        return count;
    }

    /**
     * Chr	Pos	GerpScore
     * 1A	6997	0.861
     * 1A	19328	0.202
     * 1A	22122	0.861
     * ...   ....   ....
     * @param outFile
     * @param chrConvertionRule
     */
    public void write(File outFile, ChrConvertionRule chrConvertionRule){
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            bw.write("Chr\tPos\tGerpScore");
            bw.newLine();
            StringBuilder sb;
            String chr=null;
            int refPos=-1;
            for (int i = 0; i < pos.length; i++) {
                for (int j = 0; j < pos[i].size(); j++) {
                    chr=chrConvertionRule.getRefChrFromVCFChr(i+1);
                    refPos=chrConvertionRule.getRefPosFromVCFChrPos(i+1, pos[i].get(j));
                    sb=new StringBuilder();
                    sb.append(chr).append("\t").append(refPos).append("\t").append(rsScore[i].get(j));
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * if IsGene is intergenic, BioType will be "."
     * IsGene: intergenic, gene
     * BioType:  intron, ., 5'UTR, 3'UTR, CDS
     * Chr	Pos	GerpScore	IsGene	BioType
     * 1	6997	0.861	intergenic	.
     * 1	19328	0.202	intergenic	.
     * 1	22122	0.861	intergenic	.
     * 1	28355	0.239	intergenic	.
     * @param outFile
     * @param pgf
     */
    public void writeAnnotation(String outFile, PGF pgf){
        int chr=-1;
        int pos=-1;
        double gerpScore=Double.MIN_VALUE;
        int geneIndex=Integer.MIN_VALUE;
        int longestTranscriptIndex;
        int tempIndex;
        StringBuilder sb;
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            bw.write("Chr\tPos\tGerpScore\tIsGene\tBioType");
            bw.newLine();
            for (int i = 0; i < this.getPos().length; i++) {
                for (int j = 0; j < this.getPos()[i].size(); j++) {
                    chr=i+1;
                    pos=this.getPos()[i].get(j);
                    gerpScore=this.getRsScore()[i].get(j);
                    sb=new StringBuilder();
                    sb.append(chr).append("\t").append(pos).append("\t").append(gerpScore).append("\t");
                    geneIndex= pgf.getGeneIndex(chr, pos);
                    if (geneIndex<0){
                        sb.append("intergenic").append("\t").append(".");
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }else {
                        sb.append("gene").append("\t");
                        longestTranscriptIndex= pgf.getLongestTranscriptIndex(geneIndex);
                        tempIndex= pgf.get5UTRIndex(geneIndex, longestTranscriptIndex, chr, pos);
                        if (tempIndex > -1){
                            sb.append("5'UTR");
                            bw.write(sb.toString());
                            bw.newLine();
                            continue;
                        }
                        tempIndex=pgf.get3UTRIndex(geneIndex, longestTranscriptIndex, chr, pos);
                        if (tempIndex > -1){
                            sb.append("3'UTR");
                            bw.write(sb.toString());
                            bw.newLine();
                            continue;
                        }
                        tempIndex=pgf.getCDSIndex(geneIndex, longestTranscriptIndex, chr, pos);
                        if (tempIndex > -1){
                            sb.append("CDS");
                            bw.write(sb.toString());
                            bw.newLine();
                            continue;
                        }
                        sb.append("intron");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }

    }

}
