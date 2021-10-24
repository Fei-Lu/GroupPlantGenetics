package daxing.common.factors;

import java.util.Arrays;

public enum HexaploidBySubcontinent {

    LR_AF, LR_AM, LR_CSA, LR_EA, LR_EU, LR_WA, CL;

    public int getIndex(){
        return Arrays.binarySearch(HexaploidBySubcontinent.values(), this);
    }
}
