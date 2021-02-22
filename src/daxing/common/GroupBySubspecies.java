package daxing.common;

public enum GroupBySubspecies {

    WE("Wild_emmer"),DE("Domesticated_emmer"),FT("Free_threshing_tetraploid"),AT("Ae.tauschii"),
    LR("Landrace"), CL("Cultivar");

    String group;

    GroupBySubspecies(String group) {
        this.group=group;
    }

    public String getGroup() {
        return group;
    }

    public static GroupBySubspecies newInstanceFrom(String group){
        switch (group){
            case "Wild_emmer":
                return WE;
            case "Domesticated_emmer":
                return DE;
            case "Free_threshing_tetraploid":
                return FT;
            case "Ae.tauschii":
                return AT;
            case "Landrace":
                return LR;
            case "Cultivar":
                return CL;
        }
        return null;
    }

}
