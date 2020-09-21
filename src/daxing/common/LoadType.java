package daxing.common;

public enum LoadType {

    Syn("SYNONYMOUS",0), Non("NONSYNONYMOUS",1), Del("DELETERIOUS",2);

    public static final LoadType SYNONYMOUS= LoadType.Syn;
    public static final LoadType NONSYNONYMOUS= LoadType.Non;
    public static final LoadType DELETERIOUS= LoadType.Del;


    String loadType;
    int index;
    LoadType(String loadType, int index) {
        this.loadType=loadType;
        this.index=index;
    }

    public String getLoadType() {
        return loadType;
    }

    public int getIndex() {
        return index;
    }

    public static LoadType newInstanceFrom(String loadType){
        switch (loadType.toUpperCase()) {
            case "SYNONYMOUS":
                return SYNONYMOUS;
            case "NONSYNONYMOUS":
                return NONSYNONYMOUS;
            case "DELETERIOUS":
                return DELETERIOUS;
            default:
                System.out.println("please check your parameter: "+loadType);
        }
        return null;
    }
}
