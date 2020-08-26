package daxing.common;

/**
 * @author xudaxing
 *
 */
public enum Ploidy {

    hexaploid("hexaploid", 3, "AABBDD","ABD"),
    tetraploid("tetraploid", 2, "AABB","AB"),
    diploid("diploid", 1, "DD","D");

    public static final Ploidy HEXAPLOID= Ploidy.hexaploid;
    public static final Ploidy TETRAPLOID= Ploidy.tetraploid;
    public static final Ploidy DIPLOID= Ploidy.diploid;

    String ploidy;
    int subNum;
    String subChar;
    String subSingleChr;

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

    public static Ploidy newInstanceFromPloidy(String ploidy){
        switch (ploidy.toLowerCase()) {
            case "hexaploid":
                return HEXAPLOID;
            case "tetraploid":
                return TETRAPLOID;
            case "diploid":
                return DIPLOID;
            default:
                System.out.println("please check your parameter: "+ploidy);
        }
        return null;
    }

    public static Ploidy newInstanceFromSubChar(String subgenomeChar){
        switch (subgenomeChar.toUpperCase()){
            case "AABBDD":
                return HEXAPLOID;
            case "AABB":
                return TETRAPLOID;
            case "DD":
                return DIPLOID;
            default:
                System.out.println("please check your parameter: "+subgenomeChar);
        }
        return null;
    }

    public static Ploidy newInstanceFromSubSingleChar(String subgenomeSingleChar){
        switch (subgenomeSingleChar.toUpperCase()){
            case "ABD":
                return HEXAPLOID;
            case "AB":
                return TETRAPLOID;
            case "D":
                return DIPLOID;
            default:
                System.out.println("please check your parameter: "+subgenomeSingleChar);
        }
        return null;
    }

    public static Ploidy newInstanceFromSubNum(int subgenomeNum){
        switch (subgenomeNum){
            case 3:
                return HEXAPLOID;
            case 2:
                return TETRAPLOID;
            case 1:
                return DIPLOID;
            default:
                System.out.println("please check your parameter: "+subgenomeNum);
        }
        return null;
    }

}
