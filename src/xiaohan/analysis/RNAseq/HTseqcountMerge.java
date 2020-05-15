/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;

import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author yxh
 */
public class HTseqcountMerge {
    public HTseqcountMerge() {
        this.CountMerge();
//        this.Merge();
    }
//  
    public void CountMerge() {//顺序是HTSeq里面的顺序
        String DirS = "/Users/yxh/Documents/eQTL/009ERCC test/SiPAS-Truseq_1M/SiPAS1M/geneCount";
//        String inputDirS = "/Users/yxh/Documents/RNA-seq/006test/HTseqcount/coleoptile";
//        String outputDirS = "/Users/yxh/Documents/RNA-seq/006test/DEseq2";
//        String inputDirS = "/data1/home/xiaohan/rnaseq/root/HTseqcount";
//        String outputDirS = "/data1/home/xiaohan/rnaseq/root/DESeq2";
//        String inputDirS = "/Users/yxh/Documents/RNA-seq/test/seq/coleoptile";
        List<String> nameList = new ArrayList<>();
        String subFqDirS = new File(DirS).getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        Arrays.sort(fs);
        final List<File> fList = Arrays.asList(fs);
        int[][] count = new int[92][fList.size()];
        
        List<Integer> indexList = new ArrayList<>();
        for (int i = 0; i < fList.size(); i++) {
            indexList.add(i);
        }
        indexList.stream().forEach(index -> {
            File f = fList.get(index);
            String temp = null;
            String[] tem = null;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                int cnt = 0;int counter = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> tList = PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    if (tem[0].startsWith("ERCC")) {
                        nameList.add(tem[0]);
                        int index1 = counter;
                        count[index1][index] = Integer.parseInt(tem[1]);
                        counter++;
                    }
                    cnt++;
                    if (cnt%10000 == 0) System.out.println(cnt);
                }
                br.close();
            } catch (Exception ex) {
                System.out.println(tem[0] + "\t1234");
                ex.printStackTrace();

            }
        });
        
                String outputFileS = new File(DirS, "countmergecol.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            sb.append("Gene" + "\t");
            for (int i = 0; i < fList.size(); i++) {
                sb.append(fList.get(i).getName().replace("Count.txt", "") + "\t");
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < count.length; i++) {
                sb = new StringBuilder();
                for (int j = 0; j < fList.size(); j++) {
                    if (j == 0) {
                        sb.append(nameList.get(i) + "\t");
                    }
                    sb.append(count[i][j] + "\t");
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new HTseqcountMerge();
    }

        
          
//    public void Merge(){
//        String inputDirS = "/Users/yxh/Documents/RNA-seq/test/seq/root";
//        String outputDirS = "/Users/yxh/Documents/RNA-seq/test/seq/root/count";
//        String OriginFile = "/Users/yxh/Documents/RNA-seq/test/seq/root";
//        RowTable<String> t = new RowTable<>(OriginFile);
//        String[] namelist = new String[107891];
//        namelist[0] = "Gene";
//        for(int i = 0 ;i<t.getRowNumber();i++){
//          namelist[i+1] = t.getCell(i, 0).toString(); 
//        }
//        List[] l = new List[96];
//        for(int i = 0 ; i < l.length ; i ++){
//            l[i] = new ArrayList();
//        }
//        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, ".count");
//        List<File> fList = Arrays.asList(fs);
//        HashSet<String> nameSet = new HashSet<String>();
//        for (int i = 0; i < fs.length; i++) {
//            nameSet.add(fs[i].getName().replace(".count", ""));
//        }
//        nameSet.stream().forEach((String p) -> {  
//            HashMap countStrain = new HashMap();
//            RowTable<String> tt = new RowTable<>(new File(inputDirS,p+".count").getAbsolutePath());
//            for(int i = 0;i < tt.getRowNumber();i++){
//               countStrain.put(t.getCell(i, 1),t.getCell(i, 0)); 
//            }
//            StringBuilder sb = new StringBuilder();
//            sb.append(nameSet.toArray().toString());
//            BufferedWriter [] bw = new BufferedWriter[fList.size()];
//            for(int j = 0;j<fList.size();j++){
//            bw [j]= IOUtils.getTextWriter(new File(outputDirS,"countResult").getAbsolutePath());
//            }
//            
//        });
//    }
//    
//    public static File[] listFilesEndsWith (File[] fAll, String endStr) {
//        ArrayList<File> al = new ArrayList();
//        for (int i = 0; i < fAll.length; i++) {
//            if (fAll[i].getName().endsWith(endStr)){ 
//                al.add(fAll[i]);
//            }
//        }
//        return al.toArray(new File[al.size()]);
//    }

//        fList.parallelStream().forEach(f -> {
//            String temp = null;
//            String[] tem = null;
//            try {
//                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
//                while ((temp = br.readLine()) != null) {
//                    List<String> tList = PStringUtils.fastSplit(temp);
//                    tem = tList.toArray(new String[tList.size()]);
//                    //tem =  temp.split(temp);
////                    nameList.add(tem[0]);
//                    if (tem[0].startsWith("TraesCS")) {
//                        if (!nameList.contains(tem[0])) {
//                            nameList.add(tem[0]);
//                        }
//                        int index = nameList.indexOf(tem[0]);
//                        //int index2 = fList.indexOf();
//                        count[index][fList.indexOf(f.getName())] = Integer.parseInt(tem[1]);
//                    }
//                }
//            } catch (Exception ex) {
//                System.out.println(tem[0] + "\t1234");
//                ex.printStackTrace();
//
//            }
//        });


}
