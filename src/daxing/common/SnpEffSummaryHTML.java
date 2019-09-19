package daxing.common;

import gnu.trove.list.array.TIntArrayList;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * SnpEff version: SnpEff 4.3t (build 2017-11-24 10:18), by Pablo Cingolani
 */
public class SnpEffSummaryHTML {

    private int numberOfVariants;
    private int numberOfEffects;
    private int chrLen;
    private int[] numberOfEffectsByImpact;   //HIGH, LOW, MODERATE, MODIFIER
    private int[] numberOfEffectsByFunctionalClass; //MISSENSE, NONSENSE, SILENT
    private int[] numberOfEffectsByType;  //3_prime_UTR_variant, 5_prime_UTR_premature_start_codon_gain_variant, ...
    private int[] numberOfEffectsByRegion;  //DOWNSTREAM, EXON, INTERGENIC, INTRON, ...

    public SnpEffSummaryHTML(File inputFileOfSummaryHTML){
        this.parseSummaryHTML(inputFileOfSummaryHTML);
    }

    private void parseSummaryHTML(File inputFileOfSummaryHTML){
        try (BufferedReader br = IOUtils.getTextReader(inputFileOfSummaryHTML.getAbsolutePath())) {
            String line;
            boolean first=false;
            TIntArrayList numberOfEffectsByImpactList=new TIntArrayList();
            TIntArrayList numberOfEffectsByFunctionalClassList=new TIntArrayList();
            TIntArrayList numberOfEffectsByTypeList=new TIntArrayList();
            TIntArrayList numberOfEffectsByRegionList=new TIntArrayList();
            while ((line = br.readLine()) != null) {
                if (line.contains("Number of variants (before filter)")){
                    line=br.readLine();
                    numberOfVariants=Integer.parseInt(line.replaceAll("\\D+", ""));
                    first=true;
                }
                if (line.contains("Number of effects") && first){
                    line=br.readLine();
                    numberOfEffects=Integer.parseInt(line.replaceAll("\\D+", ""));
                    first=false;
                    br.readLine();br.readLine();br.readLine();
                    br.readLine();br.readLine();br.readLine();
                    br.readLine();
                    line=br.readLine();
                    chrLen=Integer.parseInt(line.replaceAll("\\D+", ""));
                }
                if (line.contains("HIGH")){
                    numberOfEffectsByImpactList.add(this.extractKeyWord_1(br));
                    numberOfEffectsByImpactList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByImpactList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByImpactList.add(this.extractKeyWord_5(br));
                }
                if (line.contains("MISSENSE")){
                    numberOfEffectsByFunctionalClassList.add(this.extractKeyWord_1(br));
                    numberOfEffectsByFunctionalClassList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByFunctionalClassList.add(this.extractKeyWord_5(br));
                }
                if (line.contains("3_prime_UTR_variant")) {
                    numberOfEffectsByTypeList.add(this.extractKeyWord_1(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                }
                if (line.contains("DOWNSTREAM")){
                    numberOfEffectsByRegionList.add(this.extractKeyWord_1(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                }
            }
            numberOfEffectsByImpact=numberOfEffectsByImpactList.toArray();
            numberOfEffectsByFunctionalClass=numberOfEffectsByFunctionalClassList.toArray();
            numberOfEffectsByType=numberOfEffectsByTypeList.toArray();
            numberOfEffectsByRegion=numberOfEffectsByRegionList.toArray();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private int extractKeyWord_1(BufferedReader br){
        int num=Integer.MIN_VALUE;
        try {
            String line;
            String temp=null;
            br.readLine();
            line=br.readLine();
            Pattern pattern=Pattern.compile("([1-9]+[,]?)+\\d*");
            Matcher matcher=pattern.matcher(line);
            while (matcher.find()){
                temp=matcher.group(0);
            }
            List<String> l= PStringUtils.fastSplit(temp, ",");
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < l.size(); i++) {
                sb.append(l.get(i));
            }
            num=Integer.parseInt(sb.toString());
        } catch (IOException e) {
            e.printStackTrace();
        }
        return num;
    }

    private int extractKeyWord_5(BufferedReader br){
        int num=Integer.MIN_VALUE;
        try {
            String line;
            String temp=null;
            br.readLine();
            br.readLine();
            br.readLine();
            br.readLine();
            br.readLine();
            line=br.readLine();
            Pattern pattern=Pattern.compile("([1-9]+[,]?)+\\d*");
            Matcher matcher=pattern.matcher(line);
            while (matcher.find()){
                temp=matcher.group(0);
            }
            List<String> l= PStringUtils.fastSplit(temp, ",");
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < l.size(); i++) {
                sb.append(l.get(i));
            }
            num=Integer.parseInt(sb.toString());
        } catch (IOException e) {
            e.printStackTrace();
        }
        return num;
    }

    public int getNumberOfVariants() {
        return numberOfVariants;
    }

    public int getNumberOfEffects() {
        return numberOfEffects;
    }

    public int getChrLen() {
        return chrLen;
    }

    public int[] getNumberOfEffectsByImpact() {
        return numberOfEffectsByImpact;
    }

    public int[] getNumberOfEffectsByFunctionalClass() {
        return numberOfEffectsByFunctionalClass;
    }

    public int[] getNumberOfEffectsByType() {
        return numberOfEffectsByType;
    }

    public int[] getNumberOfEffectsByRegion() {
        return numberOfEffectsByRegion;
    }

    public double getVariantsRate(){
        return (double) chrLen/numberOfVariants;
    }

    public double getMissenseSilentRatio(){
        return (double) numberOfEffectsByFunctionalClass[0]/numberOfEffectsByFunctionalClass[2];
    }

    public double[] getImpactPercent(){
        double sum= Arrays.stream(this.numberOfEffectsByImpact).sum();
        return Arrays.stream(this.numberOfEffectsByImpact).mapToDouble(e->e/sum).toArray();
    }

    public double[] getFunctionPercent(){
        double sum= Arrays.stream(numberOfEffectsByFunctionalClass).sum();
        return Arrays.stream(this.numberOfEffectsByFunctionalClass).mapToDouble(e->e/sum).toArray();
    }

    public double[] getEffectsTypePercent(){
        double sum= Arrays.stream(numberOfEffectsByType).sum();
        return Arrays.stream(this.numberOfEffectsByType).mapToDouble(e->e/sum).toArray();
    }

    public double[] getEffectsRegionPercent(){
        double sum= Arrays.stream(numberOfEffectsByRegion).sum();
        return Arrays.stream(this.numberOfEffectsByRegion).mapToDouble(e->e/sum).toArray();
    }


}
