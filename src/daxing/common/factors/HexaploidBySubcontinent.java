package daxing.common.factors;

public enum HexaploidBySubcontinent {

    LR_AF(0), LR_AM(1), LR_CSA(2), LR_EA(3), LR_EU(4), LR_WA(5), CL(6);

    private final int index;

    HexaploidBySubcontinent(int index){
        this.index=index;
    }

    public int getIndex(){
        return index;
    }
}
