package daxing.common.vmap2Group;

public enum GroupByPloidy implements GroupType{

    DIPLOID("DD"), TETRAPLOID("AABB"), HEXAPLOID("AABBDD");

    String group;

    GroupByPloidy(String group) {
        this.group=group;
    }

    @Override
    public String getGroupAbbrev() {
        return this.name();
    }

    @Override
    public String getGroup() {
        return group;
    }

    public static GroupByPloidy newInstanceFrom(String group){
        switch (group){
            case "DD":
                return DIPLOID;
            case "AABB":
                return TETRAPLOID;
            case "AABBDD":
                return HEXAPLOID;
        }
        return null;
    }
}
