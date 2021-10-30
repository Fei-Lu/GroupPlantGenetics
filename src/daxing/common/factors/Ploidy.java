package daxing.common.factors;

import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author xudaxing
 *
 */
public enum Ploidy {

    HEXAPLOID("hexaploid", 3, "AABBDD","ABD"),
    TETRAPLOID("tetraploid", 2, "AABB","AB"),
    DIPLOID("diploid", 1, "DD","D");

    private final String ploidy;
    private final int subNum;
    private final String subChar;
    private final String subSingleChr;

    Ploidy(String ploidy, int subgenomowNum, String subChar, String subSingleChar) {
        this.ploidy=ploidy.toLowerCase();
        this.subNum=subgenomowNum;
        this.subChar=subChar;
        this.subSingleChr = subSingleChar;
    }

    public String getPloidy() {
        return ploidy;
    }

    public int getSubgenomewNum() {
        return subNum;
    }

    public String getSubChar() {
        return subChar;
    }

    public String getSubSingleChr() {
        return subSingleChr;
    }

    public String[] getSubgenomeArray(){
        String subSingleChr= this.subSingleChr;
        char[] temp= subSingleChr.toCharArray();
        String[] res=new String[temp.length];
        for (int i = 0; i < res.length; i++) {
            res[i]=String.valueOf(temp[i]);
        }
        return res;
    }

    private static final Map<String,Ploidy> ploidyToEnumMap=
            Stream.of(values()).collect(Collectors.toMap(Ploidy::getPloidy, e->e));

    public Optional<Ploidy> getInstanceFromPloidy(String ploidy){
        return Optional.ofNullable(ploidyToEnumMap.get(ploidy));
    }

    private static final Map<String,Ploidy> subToEnumMap=Stream.of(values()).collect(Collectors.toMap(Ploidy::getSubChar, e->e));

    public static Optional<Ploidy> getInstanceFromSubChar(String subgenomeChar){
        return Optional.ofNullable(subToEnumMap.get(subgenomeChar));
    }

    private static final Map<String,Ploidy> subSingleToEnumMap=Stream.of(values()).collect(Collectors.toMap(Ploidy::getSubSingleChr, e->e));

    public static Optional<Ploidy> getInstanceFromSubSingleChar(String subgenomeSingleChar){
        return Optional.ofNullable(subSingleToEnumMap.get(subgenomeSingleChar));
    }

    private static Map<Integer, Ploidy> subgenomeNumToEnumMap= Stream.of(values()).collect(Collectors.toMap(Ploidy::getSubgenomewNum, e->e));

    public static Optional<Ploidy> getInstanceFromSubNum(int subgenomeNum){
        return Optional.ofNullable(subgenomeNumToEnumMap.get(subgenomeNum));
    }

}
