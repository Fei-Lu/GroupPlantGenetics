package daxing.common.vmap2Group;

import daxing.common.factors.SubgenomeCombination;
import daxing.common.factors.WheatLineage;
import java.util.EnumSet;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public enum GroupBySubcontinent implements GroupType{

    WE("Wild_emmer",0, SubgenomeCombination.AB,"WE"),
    DE("Domesticated_emmer",1, SubgenomeCombination.AB,"DE"),
    FTT("Free_threshing_tetraploid",2, SubgenomeCombination.AB, "FTT"),
    AT("Ae.tauschii",3, SubgenomeCombination.D,"AT"),
    LR_EU("LR_EU",4, SubgenomeCombination.ABD,"BW"),
    LR_WA("LR_WA",5,SubgenomeCombination.ABD,"BW"),
    LR_AF("LR_Africa",6,SubgenomeCombination.ABD,"BW"),
    LR_CSA("LR_CSA",7,SubgenomeCombination.ABD,"BW"),
    LR_AM("LR_America",8,SubgenomeCombination.ABD,"BW"),
    LR_EA("LR_EA",9,SubgenomeCombination.ABD,"BW"),
    CL("Cultivar",10,SubgenomeCombination.ABD,"BW");

    String group;
    int index;
    SubgenomeCombination subgenomeCombination;
    String group_duplicated;

    GroupBySubcontinent(String group, int index, SubgenomeCombination subgenomeCombination, String group_duplicated) {
        this.group=group;
        this.index=index;
        this.subgenomeCombination=subgenomeCombination;
        this.group_duplicated = group_duplicated;
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

    public String getGroup_duplicated() {
        return group_duplicated;
    }

    public SubgenomeCombination getSubgenomeCombination() {
        return subgenomeCombination;
    }

    public static EnumSet<GroupBySubcontinent> getSubgenomeGroupByContinent(WheatLineage subgenome){
        if (subgenome==WheatLineage.A || subgenome==WheatLineage.B){
            return EnumSet.complementOf(EnumSet.of(GroupBySubcontinent.AT));
        }else if (subgenome==WheatLineage.D){
            return EnumSet.range(GroupBySubcontinent.AT, GroupBySubcontinent.CL);
        }
        return null;
    }

    private static Map<String, GroupBySubcontinent> groupToEnumMap= Stream.of(values()).collect(Collectors.toMap(GroupBySubcontinent::getGroup, e->e));

    public static Optional<GroupBySubcontinent> newInstanceFrom(String group){
        return Optional.ofNullable(groupToEnumMap.get(group));
    }
}
