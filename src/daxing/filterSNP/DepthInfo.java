package daxing.filterSNP;

import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class DepthInfo {

    private List<Dot> dotList;

    public DepthInfo(String inputFile){
        this.initialize(inputFile);
    }

    private void initialize(String inputFile){
        List<Dot> dotList=new ArrayList<>();
        String line=null;
        try(BufferedReader br= IOUtils.getTextGzipReader(inputFile)){
            List<String> lineList;
            short chr;
            int pos;
            double depth;
            double sd;
            br.readLine();
            while ((line=br.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                chr=Short.parseShort(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                depth=Double.parseDouble(lineList.get(2));
                sd=Double.parseDouble(lineList.get(3));
                dotList.add(new Dot(chr, pos, depth, sd));
            }
            Collections.sort(dotList);
            this.dotList=dotList;
        }catch (Exception e){
            e.printStackTrace();
            System.out.println(inputFile+"\t"+line);
        }
    }

    public List<Dot> getDotList() {
        return dotList;
    }
}
