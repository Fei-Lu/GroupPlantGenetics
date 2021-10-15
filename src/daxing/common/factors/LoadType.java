package daxing.common.factors;

public enum LoadType {

    Syn("SYNONYMOUS",0), Non("NONSYNONYMOUS",1), Del("DELETERIOUS",2);

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
                return LoadType.Syn;
            case "NONSYNONYMOUS":
                return LoadType.Non;
            case "DELETERIOUS":
                return LoadType.Del;
            default:
                System.out.println("please check your parameter: "+loadType);
        }
        return null;
    }
}
