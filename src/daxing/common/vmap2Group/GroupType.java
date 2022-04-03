package daxing.common.vmap2Group;

public interface GroupType {

    String getGroupAbbrev();

    String getGroup();

    static GroupType newInstanceFrom(String groups, Group group){
        switch (group){
            case Subspecies:
                return GroupBySubspecies.newInstanceFrom(groups);
            case Subcontinent:
                return GroupBySubcontinent.newInstanceFrom(groups).get();
            case Ploidy:
                return GroupByPloidy.newInstanceFrom(groups);
        }
        return null;
    }
}
