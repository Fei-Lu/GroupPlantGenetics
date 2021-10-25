package daxing.individualIntrogression;

import java.util.Arrays;

public enum P2 {

    CULTIVAR("CL"), LANDRACE("LR");

    String abbreviation;

    P2(String abbreviation) {
        this.abbreviation=abbreviation;
    }

    public int getIndex() {
        return Arrays.binarySearch(P2.values(),this);
    }

    public String getAbbreviation() {
        return abbreviation;
    }

    public static P2 newInstanceFrom(int index){
        return P2.values()[index];
    }

    public static P2 newInstanceFrom(String abbreviation){
        switch (abbreviation) {
            case "CL":
                return P2.CULTIVAR;
            case "LR":
                return P2.LANDRACE;
            default:
                System.out.println("please check your parameter: "+abbreviation);
        }
        return null;
    }
}
