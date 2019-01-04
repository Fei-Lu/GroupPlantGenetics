package daxing.md5;

import utils.IOUtils;
import utils.PArrayUtils;
import javax.xml.bind.DatatypeConverter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.security.MessageDigest;
import java.util.*;

/**
 *
 * @author xudaxing
 */
public class MD5 {
    public static String getMD5FromString(String str) {
        String myMD5=null;
        try {
            MessageDigest md = MessageDigest.getInstance("md5");
            md.update(Byte.parseByte(str));
            byte[] digest = md.digest();
            myMD5 = DatatypeConverter.printHexBinary(digest).toUpperCase();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return myMD5;
    }

    public static boolean checkMD5ForString(String str,String hash){
        String strHash=MD5.getMD5FromString(str);
        return strHash.equals(hash);
    }

    public static String getMD5FromFile(File file){
        String myMD5=null;
        try{
            MessageDigest md=MessageDigest.getInstance("md5");
            md.update(Files.readAllBytes(file.toPath()));
            byte[] digest=md.digest();
            myMD5=DatatypeConverter.printHexBinary(digest).toUpperCase();
        }
        catch (Exception e){
            e.printStackTrace();
        }
        return myMD5;
    }

    public static boolean checkMD5ForFile(File file, String hash){
        String fileHash=MD5.getMD5FromFile(file);
        return fileHash.equals(hash);
    }

    public static void getMD5FromDir(File inputDirS, File outfileMD5){
        int numThreads = 32;
        if(Runtime.getRuntime().availableProcessors()<32) {
            numThreads=Runtime.getRuntime().availableProcessors();
        }
        File[] file = IOUtils.listRecursiveFiles(inputDirS);
        if(file.length<numThreads){
            numThreads=file.length;
        }
        Map<String,String> md5ValuePathMap=new Hashtable<>();
        BufferedWriter bw=null;
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetNumber(file.length, numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            Arrays.stream(subLibIndices)
                    .filter(index-> (!(file[index].getName().contains(".DS_Store"))))
                    .forEach(index-> {
                        String md5Value=MD5.getMD5FromFile(file[index]);
                        String fName = file[index].getAbsolutePath();
                        md5ValuePathMap.put(md5Value, fName);
                    });
        }
        try{
            bw=IOUtils.getTextWriter(outfileMD5.getAbsolutePath());
            for(Map.Entry<String,String> entry:md5ValuePathMap.entrySet()){
                //System.out.println(entry.getKey()+"  "+entry.getValue());
                bw.write(entry.getKey()+"  "+entry.getValue());
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    public static void checkMD5ForDir(File hashFileS){
        int numThreads = 32;
        if(Runtime.getRuntime().availableProcessors()<32) {
            numThreads=Runtime.getRuntime().availableProcessors();
        }
        List<String> md5ValuePath=new ArrayList<>();
        String line;
        try(BufferedReader br=IOUtils.getTextReader(hashFileS.getAbsolutePath())){
            while((line=br.readLine())!=null){
                md5ValuePath.add(line);
            }
        }
        catch(Exception e){
            e.printStackTrace();
        }
        if(md5ValuePath.size()<numThreads){
            numThreads=md5ValuePath.size();
        }
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetNumber(md5ValuePath.size(), numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            Arrays.stream(subLibIndices).forEach(index->{
                String valuePath=md5ValuePath.get(index);
                String[] md5Value=valuePath.split("  ");
                String value=md5Value[0];
                String path=md5Value[1];
                boolean f=MD5.checkMD5ForFile(new File(path), value);
                if(!f){
                    System.out.println("False  "+path+": "+value);
                }
                else{
                    System.out.println("True  "+path+": "+value);
                }
            });
        }
    }

}

