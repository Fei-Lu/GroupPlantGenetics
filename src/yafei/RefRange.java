/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yafei;
import daxing.common.StringTool;
import gnu.trove.list.array.TIntArrayList;
//import utils.IOUtils;

import java.io.*;

import static java.lang.Math.abs;
import pgl.infra.utils.IOUtils;
/**
 *
 * @author guoyafei
 * 
 */
public class RefRange extends StringTool {
        /**
        * 以行为单位读取文件        
        * @param path
        * @throws java.io.IOException        
        */
        public static void readFromFiles(String path) throws IOException {           
            File writename = new File("/data2/yafei/output"); // 相对路径，如果没有则要建立一个新的output.txt文件
	    writename.createNewFile(); // 创建新文件
	    BufferedWriter out = new BufferedWriter(new FileWriter(writename));
            File file = new File(path);
            String[] fileList = file.list(); 
            
            for(int i=0;i < fileList.length;i++){
                BufferedReader reader = null;
                try {
                    String filename = "/data1/publicData/wheat/reference/v1.0/byChr/"+ (String) fileList[i];
                    reader = IOUtils.getTextGzipReader(filename);
                    String chr = reader.readLine();
                    StringBuilder sb = new StringBuilder(500_000_000);
                    String str1 = reader.readLine();
                    while (str1 != null) {
//                    for(i=0;i<5000;i++){
                        sb.append(str1);
                        str1 = reader.readLine();
                    }
                    String seq = sb.toString();
                    
//                    System.out.println(seq.length());                    
                    TIntArrayList index1 = getIndexOfSubStr(seq, "CCGG");
//                    int g = getCountOfSubStr(seq, "GCC");
                    TIntArrayList index2 = getIndexOfSubStr(seq, "GGATCC");
//                    int y = getCountOfSubStr(seq, "GATCC");                    
//                    System.out.println(index1.get(300000)+" "+g+" "+index2.get(3000)+" "+y);  
                    int[] temp1 = new int[2];
                    int[] temp2 = new int[2];
                    for (int m = 0; m < index1.size(); m++) {
                        for (int n = 0; n < index2.size(); n++){
//                            System.out.println(index1.get(m)+" "+index2.get(n));
                            int abs = abs(index1.get(m)-index2.get(n));
//                            System.out.println(index1.get(m)+" "+abs);
                            if (abs > 1000) {
                                if(index1.get(m)-index2.get(n) > 0){
                                    continue;
                                }else{
                                    break;
                                }
//                        }
                            }else{
                                int min = (index1.get(m) - index2.get(n) < 0) ? index1.get(m) : index2.get(n);
                                int max = (index1.get(m) - index2.get(n) < 0) ? index2.get(n) : index1.get(m);
                                if(temp1[0] != 0){
                                    temp2[0] = min;
                                    temp2[1] = max;
                                }else{
                                    temp1[0] = min;
                                    temp1[1] = max;
                                }
                                if(temp2[0] != 0 & temp2[0] <= temp1[1]){
                                    temp1[1] = temp2[1];
                                }else{
                                    if(temp2[0] == 0){
                                        continue;
                                    }else{
                                    min = temp1[0];
                                    max = temp1[1];
                                    sb=new StringBuilder();
                                    sb.append(chr).append("\t").append(min).append("\t").append(max).append("\n");
//                                    System.out.print(chr + " " + min + " " + max + "\r\n");
                                    temp1[0] = temp2[0];
                                    temp1[1] = temp2[1];
                                    out.write(sb.toString());
                                    }
                                }
                            }
                        }
                    } 
                    out.flush();
                }catch (IOException e) {
                } finally {
                    if (reader != null) {
                        try {
                            reader.close();
                        } catch (IOException e1) {
                        }
                    }
                }
            }
            out.close();    
        }
}