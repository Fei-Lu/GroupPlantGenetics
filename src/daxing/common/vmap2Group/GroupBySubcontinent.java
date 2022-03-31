package daxing.common.vmap2Group;

import daxing.common.factors.WheatLineage;

import java.util.EnumSet;

public enum GroupBySubcontinent implements GroupType{

    WE("Wild_emmer",0),DE("Domesticated_emmer",1),FT("Free_threshing_tetraploid",2),AT("Ae.tauschii",3),
    LR_EU("LR_EU",4),
    LR_WA("LR_WA",5),LR_AF("LR_Africa",6),LR_CSA("LR_CSA",7),LR_AM("LR_America",8),
    LR_EA("LR_EA",9), CL("Cultivar",10);

    String group;
    int index;

    GroupBySubcontinent(String group) {
        this.group=group;
    }
    GroupBySubcontinent(String group, int index) {
        this.group=group;
        this.index = index;
    }

    @Override
    public String getGroup() {
        return group;
    }

    @Override
    public String getGroupAbbrev() {
        return this.name();
    }

    public int getIndex() {
        return index;
    }

    public static EnumSet<GroupBySubcontinent> getSubgenomeGroupByContinent(WheatLineage subgenome){
        if (subgenome==WheatLineage.A || subgenome==WheatLineage.B){
            return EnumSet.complementOf(EnumSet.of(GroupBySubcontinent.AT));
        }else if (subgenome==WheatLineage.D){
            return EnumSet.range(GroupBySubcontinent.AT, GroupBySubcontinent.CL);
        }
        return null;
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
