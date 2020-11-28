package daxing.individualIntrogression;

public enum P2 {

    Cultivar(0, "CL"), Landrace(1,"LR");

    public static final P2 CL= P2.Cultivar;
    public static final P2 LR= P2.Landrace;

    int index;
    String abbreviation;

    P2(int index, String abbreviation) {
        this.index=index;
        this.abbreviation=abbreviation;
    }

    public int getIndex() {
        return index;
    }

    public String getAbbreviation() {
        return abbreviation;
    }

    public static final P2 newInstanceFrom(int index){
        switch (index) {
            case 0:
                return CL;
            case 1:
                return LR;
            default:
                System.out.println("please check your parameter: "+index);
        }
        return null;
    }

    public static final P2 newInstanceFrom(String abbreviation){
        switch (abbreviation) {
            case "CL":
                return CL;
            case "LR":
                return LR;
            default:
                System.out.println("please check your parameter: "+abbreviation);
        }
        return null;
    }
}
