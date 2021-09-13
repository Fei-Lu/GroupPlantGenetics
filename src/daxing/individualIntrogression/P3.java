package daxing.individualIntrogression;

public enum P3 {

    WILD_EMMER(0, "WE"), DOMESTICATED_EMMER(1, "DE"), FREE_THRESHING_TETRAPLOIDS(2, "FT"), STRANGULATA(3, "AT");

    public static final P3 WE= P3.WILD_EMMER;
    public static final P3 DE= P3.DOMESTICATED_EMMER;
    public static final P3 FT= P3.FREE_THRESHING_TETRAPLOIDS;
    public static final P3 AE= P3.STRANGULATA;

    int index;
    String abbreviation;

    P3(int index, String p3Abbreviation) {
        this.index=index;
        this.abbreviation=p3Abbreviation;
    }

    public int getIndex() {
        return index;
    }

    public String getAbbreviation() {
        return abbreviation;
    }

    public static final P3 newInstanceFrom(int index){
        switch (index) {
            case 0:
                return WE;
            case 1:
                return DE;
            case 2:
                return FT;
            case 3:
                return AE;
            default:
                System.out.println("please check your parameter: "+index);
        }
        return null;
    }

    public static final P3 newInstanceFrom(String p3Abbreviation){
        switch (p3Abbreviation) {
            case "WE":
                return WE;
            case "DE":
                return DE;
            case "FT":
                return FT;
            case "AT":
                return AE;
            default:
                System.out.println("please check your parameter: "+p3Abbreviation);
        }
        return null;
    }

}
