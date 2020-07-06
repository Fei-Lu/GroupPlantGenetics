package yafei;
//import utils.IOUtils;
import java.io.*;
import java.util.*;
import pgl.infra.utils.IOUtils;
/**
 * @author Yafei Guo
 * @create 2020-07-01 1:25 PM
 */
public class MergeVD {
    public static StringBuilder getSign(Integer num, String sign){
        StringBuilder sb = new StringBuilder();
        String s = sign + "/" + sign;
        for (int i = 0; i < num-1; i++) {
            sb.append(s+"\t");
        }
        sb.append(s);
        return sb;
    }
    public static File[] sortPathFile(File Path, String estring){
        File[] files = Path.listFiles(new FileFilter() {
            @Override
            public boolean accept(File file) {
                if(file.getName().endsWith(estring) && file.isFile()){
                    return true;
                }
                return false;   //否则过滤掉
            }
        });
        List sortfiles = Arrays.asList(files);
        Collections.sort(sortfiles, new Comparator<File>(){
            public int compare(File o1, File o2) {
                return o1.getName().compareTo(o2.getName());
            }
        });
        return files;
    }
    public static Map getDownID(String Inputfile) throws IOException {
        BufferedReader reader = IOUtils.getTextGzipReader(Inputfile);
        Map map = new HashMap();
        String line;
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("\t",2));
            map.put(list.get(0),list.get(1));
        }
        reader.close();
        return map;
    }
    public static List getID(String Inputfile, String num ) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(new File(Inputfile)));
        List listId = new ArrayList<>();
        String line;
        while ((line = reader.readLine()) != null) {
            List list = Arrays.asList(line.split("-"));
            String i = String.valueOf(list.get(0));
            if(num.equals(list.get(0))){
                listId.add(line);
            }
        }
        reader.close();
        return listId;
    }
    public static void writeOut(String inputVmap, String inputId, String output ) throws IOException {
        List listV;
        listV = Arrays.asList(sortPathFile(new File(inputVmap), "all.vcf.gz"));
        BufferedWriter writer = new BufferedWriter(new FileWriter(output));
        int count = 0;
        for (int i = 0; i <42; i++) {
            BufferedReader reader;
            reader = IOUtils.getTextGzipReader(String.valueOf(listV.get(i)));
            System.out.println(listV.get(i));
            List listid = getID(inputId, String.valueOf(i+1));
            String line;
            while ((line = reader.readLine()) != null) {
                List<String> list = Arrays.asList(line.split("\t"));
                if(list.size() > 2){
                    if(listid.contains(list.get(2))){
                        writer.write(list.get(2)+"\n");
                        count++;
                    }
                }
            }
            reader.close();
        }
        writer.write(count+"\n");
        writer.close();
    }
    public static void main(String[] args) throws IOException {
//        System.out.println(getDot(10));
        List listV;
        listV = Arrays.asList(sortPathFile(new File("/data1/home/yafei/Project3/Vmap1.1"), "all.vcf.gz"));
        List listD = Arrays.asList(sortPathFile(new File("/data1/home/yafei/Project3/Download"), "txt.gz"));
        List<String> out = new ArrayList();
        String a = "chr";
        String b = ".vcf";
        for (int i = 1; i < 43; i++) {
            String name = a + i + b;
            out.add(name);
        }
        for (int i = 0; i < 42; i++) {
            Map ID;
            ID = getDownID(String.valueOf(listD.get(i)));
            BufferedReader reader;
            reader = IOUtils.getTextGzipReader(String.valueOf(listV.get(i)));
            BufferedWriter writer = new BufferedWriter(new FileWriter("/data1/home/yafei/Project3/linageVcf/"+ out.get(i)));
            for (int j = 0; j < 22 ; j++) {
                writer.write(reader.readLine()+"\n");
            }
            String line;
            while ((line = reader.readLine()) != null) {
                List<String> list = Arrays.asList(line.split("\t"));
                if(ID.keySet().contains(list.get(2))){
                    writer.write(line + "\t" + ID.get(list.get(2))+"\n");
                }else{
                    writer.write(line + "\t" + getSign(855,".")+"\n");
                }
            }
            writer.close();
            reader.close();
        }
    }
}
