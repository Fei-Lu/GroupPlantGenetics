package daxing.individualIntrogression;

public enum P3 {

    Wild_emmer(0), Domesticated_emmer(1), Free_threshing_tetraploid(2), Ae(3);

    public static final P3 WE= P3.Wild_emmer;
    public static final P3 DE= P3.Domesticated_emmer;
    public static final P3 FT= P3.Free_threshing_tetraploid;
    public static final P3 AE= P3.Ae;

    int index;

    P3(int index) {
        this.index=index;
    }

    public int getIndex() {
        return index;
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
}
