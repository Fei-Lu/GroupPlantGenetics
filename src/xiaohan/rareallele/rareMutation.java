package xiaohan.rareallele;

import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @ author: yxh
 * @ created: 2021-07-13 : 2:21 PM
 */
public class rareMutation {

    public rareMutation(String[] args){
        this.getAnnotation(args);
    }

    public void getAnnotation(String[] args){
        HashSet<String> positionSet = new HashSet<>();
        HashMap<String,String> positionMap = new HashMap<>();
        HashMap<String,String> positionMap1 = new HashMap<>();
        String annotationFile = args[0];
        BufferedReader br = IOUtils.getTextReader(annotationFile);
        String temp = null;
        String[] temps = null;
        String position = null;
        String annotation = null;
        String annotation1 = null;
        try{
            while((temp = br.readLine())!=null){
                temps = temp.split("\t");
                if(temp.startsWith("ID"))continue;
                position = temps[0].split("-")[0]+"_"+temps[0].split("-")[1];
                annotation = temps[11];
                annotation1 = temps[12];
                positionSet.add(position);
                positionMap.put(position,annotation);
                positionMap1.put(position,annotation1);
            }
            System.out.println("Finished mapping  ----------------------------------------------------------------------");
            for (int i = 0; i < 42; i++) {
                int chr = i+1;
                String infileDir = args[1];
                String outfileDir = args[2];
                File out = new File(new File(outfileDir,"/chr"+chr).getAbsolutePath());
                out.mkdir();
                File[] fs = new File(infileDir,"/chr"+chr).listFiles();
                fs = IOUtils.listFilesEndsWith(fs, "_SNP.bed");
                for (int j = 0; j < fs.length; j++) {
                    System.out.println("Reading file : "+ fs[j].getName() + "-------------------------------------------");
                    BufferedReader br1 = IOUtils.getTextReader(fs[j].getAbsolutePath());
                    BufferedWriter bw = IOUtils.getTextWriter(new File(out,fs[j].getName()).getAbsolutePath());
                    while((temp = br1.readLine())!=null){
                        temps = temp.split("\t");
                        position = temps[0].replace("chr","")+"_"+temps[2];
                        if(positionMap.get(position)!=null){
                            bw.write(temp+"\t"+positionMap.get(position)+"\t"+positionMap1.get(position)+"\n");
                        }else {
                            bw.write(temp+"\tintergenic\tNon_conserved\n");
                        }
                    }
                    br1.close();
                    bw.flush();
                    bw.close();
                }
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void main(String[] args){
        new rareMutation(args);
    }
}
