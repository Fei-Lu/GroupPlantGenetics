/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yafei;

/**
 *
 * @author guoyafei
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
/**
 * @author guoyafei
 */
public class GBSpos {
    private static Set<String> getGBSPos(String Path1) throws FileNotFoundException, IOException {
        FileInputStream inputStream = new FileInputStream(Path1);
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        List<List<String>> result = new ArrayList(); 
        String line = null; 
        Set<String> numbers = new HashSet<String>(); 
        while ((line = reader.readLine()) != null) {                
            List<String> list = Arrays.asList(line.split("\t"));  
            StringBuilder sb = new StringBuilder();
            sb.append(list.get(0)+"_"+list.get(1));
            numbers.add(sb.toString());
        }
        reader.close();
        return numbers;    
    }  
    public static void getLibraryPos(String Path1, String Path2) throws FileNotFoundException, IOException {
        File writename = new File("/data2/yafei/Hapscanner/posAllele_Filterd_java.txt"); // 相对路径，如果没有则要建立一个新的output.txt文件
	writename.createNewFile(); // 创建新文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));   
        FileInputStream inputStream = new FileInputStream(Path2);
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        String line = null;
        
        Set<String> numbers = new HashSet<String>();
        numbers = getGBSPos(Path1);
        int m = 0;
        while ((line = reader.readLine()) != null) {                
            List<String> list = Arrays.asList(line.split("\t"));
            StringBuilder sb = new StringBuilder();
            sb.append(list.get(0)+"_"+list.get(1));
            if(numbers.contains(sb.toString())){
                out.write(line+"\n");
//                System.out.println(line);
            }
        }
        out.flush();
        out.close();
        reader.close();
    }
}
