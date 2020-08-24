package daxing.common;

/**
 * @author xudaxing
 *
 */
public enum Ploidy {

    hexaploid("hexaploid", 3, "AABBDD"),
    tetraploid("tetraploid", 2, "AABB"),
    diploid("diploid", 1, "DD");

    public static final Ploidy HEXAPLOID= Ploidy.hexaploid;
    public static final Ploidy TETRAPLOID= Ploidy.tetraploid;
    public static final Ploidy DIPLOID= Ploidy.diploid;

    String ploidy;
    int subNum;
    String subChar;

    Ploidy(String ploidy, int subgenomowNum, String subChar) {
        this.ploidy=ploidy.toLowerCase();
        this.subNum=subgenomowNum;
        this.subChar=subChar;
    }

    public String getPloidy() {
        return ploidy;
    }

    public int getSubgenomowNum() {
        return subNum;
    }

    public String getSubChar() {
        return subChar;
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

    public static Ploidy newInstanceFromSubChar(String subgenomowChar){
        switch (subgenomowChar.toUpperCase()){
            case "AABBDD":
                return HEXAPLOID;
            case "AABB":
                return TETRAPLOID;
            case "DD":
                return DIPLOID;
            default:
                System.out.println("please check your parameter: "+subgenomowChar);
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
