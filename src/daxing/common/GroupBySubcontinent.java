package daxing.common;

public enum GroupBySubcontinent {

    WE("Wild_emmer"),DE("Domesticated_emmer"),FT("Free_threshing_tetraploid"),AT("Ae.tauschii"),
    LR_EU("LR_EU"),
    LR_WA("LR_WA"),LR_AF("LR_Africa"),LR_CSA("LR_CSA"),LR_AM("LR_America"),
    LR_EA("LR_EA"), CL("Cultivar");

    String group;

    GroupBySubcontinent(String group) {
        this.group=group;
    }

    public String getGroup() {
        return group;
    }

    public static GroupBySubcontinent newInstanceFrom(String group){
        switch (group){
            case "Wild_emmer":
                return WE;
            case "Domesticated_emmer":
                return DE;
            case "Free_threshing_tetraploid":
                return FT;
            case "Ae.tauschii":
                return AT;
            case "LR_EU":
                return LR_EU;
            case "LR_WA":
                return LR_WA;
            case "LR_Africa":
                return LR_AF;
            case "LR_CSA":
                return LR_CSA;
            case "LR_America":
                return LR_AM;
            case "LR_EA":
                return LR_EA;
            case "Cultivar":
                return CL;
        }
        return null;
    }
}
