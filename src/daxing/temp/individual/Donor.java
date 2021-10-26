package daxing.temp.individual;

import java.util.Arrays;

public enum Donor{
    WE,DE,FTT,AT,NONE;

    public int getIndex(){
        return Arrays.binarySearch(Donor.values(), this);
    }
}
