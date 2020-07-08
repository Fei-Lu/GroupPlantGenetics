/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yafei;

//import utils.IOUtils;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import pgl.infra.utils.IOUtils;
/**
 *
 * @author guoyafei
 */
public class PileupCount {
    
//    private static Set<String> getVCFPos(String Path2) throws FileNotFoundException, IOException {
//        FileInputStream inputStream = new FileInputStream(Path2);
//        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
//        List<List<String>> result = new ArrayList(); 
//        String line = null; 
//        Set<String> numbers = new HashSet<String>(); 
//        while ((line = reader.readLine()) != null) {                
//            List<String> list = Arrays.asList(line.split(" "));  
//            StringBuilder sb = new StringBuilder();
//            sb.append(list.get(0)+"_"+list.get(1));
//            numbers.add(sb.toString());
//        }
//        reader.close();
//        return numbers;    
//    }  
    public static void readsCount(String Path1, String Path2) throws FileNotFoundException, IOException {
        File writename = new File("/data2/yafei/Hapscanner/output.txt"); // 相对路径，如果没有则要建立一个新的output.txt文件
	writename.createNewFile(); // 创建新文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));   
        
        List<File> fList = IOUtils.getFileListInDirEndsWith(Path1, "bam");
        
//        Set<String> numbers = new HashSet<String>();
//        numbers = getVCFPos(Path2);
        HashMap<String,Integer> hsTest= new HashMap<String,Integer>();
        
        for (int i = 0; i < fList.size(); i++) {
         
            FileInputStream inputStream = new FileInputStream(fList.get(i).getAbsolutePath());
            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
            String line = null;
   
            while ((line = reader.readLine()) != null) {                
                List<String> list = Arrays.asList(line.split(" "));
                
                StringBuilder sb = new StringBuilder();
                
                sb.append(list.get(0)+"_"+list.get(1));
                
//                if(flag){
                    
//                    hsTest.put(sb.toString(), list.get(2));
//                    
//                    hsTest.put("young", 23);
//                    
//                    hsTest.put(null, 0);
//                    
//                    System.out.println("the size of hsTest is " + hsTest.size());
//
//                    
//                    out.write(line+"\n");
////                System.out.println(line);
//                    System.out.println(sb.toString());
//                }
//                sb.append(fList.get(i).getAbsolutePath());
//            }
//        }
//
//        out.flush();
//        out.close();
//        reader.close();
//        
        }
        
    }
    
}}
