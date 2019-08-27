package daxing.ancestralAllele;

import format.position.ChrPos;
import utils.IOUtils;
import utils.PStringUtils;
import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class AllelesInfor {

    private List<ChrPos>[] chrPoss;
    private List<String>[] refAlleles;
    private List<List<String>>[] altAlleles;

    public AllelesInfor(Path vcfInputFile){
        this.initialize(vcfInputFile);
    }

    private void initialize(Path vcfInputFile){
        chrPoss=new List[45];
        refAlleles=new List[45];
        altAlleles=new List[45];
        for (int i = 0; i < chrPoss.length; i++) {
            chrPoss[i]=new ArrayList<>();
            refAlleles[i]=new ArrayList<>();
            altAlleles[i]=new ArrayList<>();
        }
        try(BufferedReader br= IOUtils.getNIOTextReader(vcfInputFile.toString())){
            String line;
            List<String> lineList;
            short chr;
            int pos;
            List<String> altAllele;
            while ((line=br.readLine())!=null){
                if (line.startsWith("#")) continue;
                lineList= PStringUtils.fastSplit(line);
                chr=Short.parseShort(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                altAllele=PStringUtils.fastSplit(lineList.get(4), ",");
                chrPoss[chr].add(new ChrPos(chr, pos));
                refAlleles[chr].add(lineList.get(3));
                altAlleles[chr].add(altAllele);
            }

        }catch (Exception e){
            e.printStackTrace();
        }

    }

    public List<ChrPos>[] getChrPoss() {
        return chrPoss;
    }

    public List<String>[] getRefAlleles() {
        return refAlleles;
    }

    public List<List<String>>[] getAltAlleles() {
        return altAlleles;
    }

    public List<ChrPos> getChrPoss(short chr){
        return this.getChrPoss()[chr];
    }

    public List<String> getRefAllele(short chr){
        return this.getRefAlleles()[chr];
    }

    public List<List<String>> getAltAllele(short chr){
        return this.getAltAlleles()[chr];
    }

}
