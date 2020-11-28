package daxing.individualIntrogression;

public enum P3 {

    Wild_emmer(0, "WE"), Domesticated_emmer(1, "DE"), Free_threshing_tetraploid(2, "FT"), Ae(3, "AT");

    public static final P3 WE= P3.Wild_emmer;
    public static final P3 DE= P3.Domesticated_emmer;
    public static final P3 FT= P3.Free_threshing_tetraploid;
    public static final P3 AE= P3.Ae;

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
