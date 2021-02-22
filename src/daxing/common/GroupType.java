package daxing.common;

public interface GroupType {

    String getGroupAbbrev();

    String getGroup();

    static GroupType newInstanceFrom(String groupType, Group group){
        switch (group){
            case Subspecies:
                return GroupBySubspecies.newInstanceFrom(groupType);
            case Subcontinent:
                return GroupBySubcontinent.newInstanceFrom(groupType);
        }
        return null;
    }
}
