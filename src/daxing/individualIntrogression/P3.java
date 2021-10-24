package daxing.individualIntrogression;

import java.util.Arrays;

public enum P3 {

    WILD_EMMER("WE"), DOMESTICATED_EMMER("DE"), FREE_THRESHING_TETRAPLOIDS("FTT"), STRANGULATA("AT");

    String abbreviation;

    P3(String p3Abbreviation) {
        this.abbreviation=p3Abbreviation;
    }

    public int getIndex() {
        return Arrays.binarySearch(P3.values(), this);
    }

    public String getAbbreviation() {
        return abbreviation;
    }

    public static P3 newInstanceFrom(int index){
        return P3.values()[index];
    }

    public static final P3 newInstanceFrom(String p3Abbreviation){
        switch (p3Abbreviation) {
            case "WE":
                return P3.WILD_EMMER;
            case "DE":
                return P3.DOMESTICATED_EMMER;
            case "FTT":
                return P3.FREE_THRESHING_TETRAPLOIDS;
            case "AT":
                return P3.STRANGULATA;
            default:
                System.out.println("please check your parameter: "+p3Abbreviation);
        }
        return null;
    }

}
