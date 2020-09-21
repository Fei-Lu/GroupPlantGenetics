package daxing.individualIntrogression;

public enum P2 {

    Cultivar(0), Landrace(1);

    public static final P2 CL= P2.Cultivar;
    public static final P2 LR= P2.Landrace;

    int index;

    P2(int index) {
        this.index=index;
    }

    public int getIndex() {
        return index;
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
}
