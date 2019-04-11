package daxing.informal;

import utils.IOUtils;
import utils.PStringUtils;
import java.io.BufferedReader;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class LibraryOfGRT {

    public static void splitLibrary(String barcodeFile, String libraryFastqMapFile){
        //
    }

    public static String getFlowCellLaneIndex(String r1FastqFile, String r2FastqFile){
        Set<String> s=new HashSet<>();
        String temp1;
        String temp2;
        StringBuilder sb1;
        StringBuilder sb2;
        try(BufferedReader brR1= IOUtils.getTextReader(r1FastqFile);
            BufferedReader brR2= IOUtils.getTextReader(r2FastqFile)){
            while ((temp1=brR1.readLine())!=null){
                temp2=brR2.readLine();
                sb1=new StringBuilder();
                sb2=new StringBuilder();
                List<String> l1=PStringUtils.fastSplit(temp1,":");
                List<String> l2=PStringUtils.fastSplit(temp2,":");
                sb1.append(l1.get(2)).append("_").append(l1.get(3)).append("_").append(l1.get(9));
                sb2.append(l2.get(2)).append("_").append(l2.get(3)).append("_").append(l2.get(9));
                s.add(sb1.toString());
                s.add(sb2.toString());
                brR1.readLine();brR1.readLine();brR1.readLine();
                brR2.readLine();brR2.readLine();brR2.readLine();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        if(s.size()==1) {
            return s.iterator().next();
        }else if(s.size()>1){
            System.out.println(r1FastqFile+" and "+r2FastqFile+" has more than one flowcellLaneIndex, please check it" );
            System.exit(1);
        }
        return null;
    }

//    public static void getModifiedBarcodeFile(String inputBarcodeFile, String inputLibraryFastqMap,
//                                              String outputModifiedBarcode, String outputModifiedLibraryFastqMap){
//        String temp1;
//        List<String> flowcellLaneIndex=new ArrayList<>();
//        try(BufferedReader brMap=IOUtils.getTextGzipReader(inputLibraryFastqMap)){
//            brMap.readLine();
//            while ((temp1=brMap.readLine())!=null){
//                List<String> l=PStringUtils.fastSplit(temp1);
//                flowcellLaneIndex.add(LibraryOfGRT.getFlowCellLaneIndex(l.get(1),l.get(2)));
//            }
//        }catch (Exception e){
//            e.printStackTrace();
//        }
//        ColumnTable<String> cBarcode=new ColumnTable<>(inputBarcodeFile);
//        ColumnTable<String> cMap=new ColumnTable<>(inputLibraryFastqMap);
//        cBarcode.removeColumn(0);
//    }
}
