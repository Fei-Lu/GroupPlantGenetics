package yafei;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Yafei Guo
 * @create 2020-07-04 9:59 AM
 * 1. 找到三个文件分离的位点，输出到Sep_ID文件中
 * 2. 找Vmap1.1与三文件并集的交集，输出到out.txt文件中
 * 3. 找到三个文件的交集
 */
public class Venn {
    public static void getSepSite(String inputFile,String outputFile) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(new File(inputFile)));
        String line = reader.readLine();
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));

        while((line = reader.readLine()) != null) {
            List list = Arrays.asList(line.split("\t",10));
            List list2 = Arrays.asList(String.valueOf(list.get(9)).split("\t"));
            Set<String> set = new HashSet<String>();
            for (int i = 0; i < list2.size(); i++) {
                List list3 = Arrays.asList(String.valueOf(list2.get(i)).split(":"));
                set.add((String) list3.get(0));
            }
            if(set.contains("./.")){
                if(set.size() == 2){
                    writer.write((String) list.get(2)+"\n");
                }
            }else{
                if(set.size() == 1){
                    writer.write((String) list.get(2)+"\n");
                }
            }

        }
        reader.close();
        writer.close();
    }
    public static void main(String[] args) throws IOException {
        getSepSite("/data1/home/yafei/Project3/Download/new_EC","/data1/home/yafei/Project3/Download/EC_noSep.txt");

    }
}
