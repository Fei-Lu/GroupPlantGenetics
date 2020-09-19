package daxing.common;

public enum LoadType {

    Syn("SYNONYMOUS"), Non("NONSYNONYMOUS"), Del("DELETERIOUS");

    public static final LoadType SYNONYMOUS= LoadType.Syn;
    public static final LoadType NONSYNONYMOUS= LoadType.Non;
    public static final LoadType DELETERIOUS= LoadType.Del;


    String loadType;
    LoadType(String loadType) {
        this.loadType=loadType;
    }

    public String getLoadType() {
        return loadType;
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
