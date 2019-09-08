package daxing.filterSNP;

import utils.IOUtils;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Predicate;

public class FilterSNPGo {

    public static void getTopDensity(String depthsDBDir, String outPutDir, double topDensity){
        File[] input= IOUtils.listRecursiveFiles(new File(depthsDBDir));
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).toArray(File[]::new);
        String[] outFileName= Arrays.stream(files).map(File::getName).map(str->str.replaceAll("depth.txt.gz$", "chrpos.txt.gz")).toArray(String[]::new);
        Map<File, String> map=new HashMap<>();
        for (int i = 0; i < files.length; i++) {
            map.put(files[i], outFileName[i]);
        }
        Arrays.stream(files).parallel().forEach(f -> {
            Cells.getTopCellDensity(f.getAbsolutePath(), new File(outPutDir, map.get(f)).getAbsolutePath(), topDensity);
        });
    }

    public static void main(String[] args) {
        FilterSNPGo.getTopDensity(args[0], args[1], Double.parseDouble(args[2]));
    }
}
