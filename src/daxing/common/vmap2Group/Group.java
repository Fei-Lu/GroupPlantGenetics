package daxing.common.vmap2Group;

import daxing.common.SubgenomeCombination;

import java.util.ArrayList;
import java.util.List;

public enum Group {
    Subspecies, Subcontinent;

    private List<String> getGroupAB() {
        List<String> groupABList=new ArrayList<>();
        switch (this){
            case Subcontinent:
                for (GroupBySubcontinent group : GroupBySubcontinent.values()){
                    if (group.equals(GroupBySubcontinent.AT)) continue;
                    groupABList.add(group.getGroup());
                }
                return groupABList;
            case Subspecies:
                for (GroupBySubspecies group : GroupBySubspecies.values()){
                    if (group.equals(GroupBySubspecies.AT)) continue;
                    groupABList.add(group.getGroup());
                }
                return groupABList;
        }
        return groupABList;
    }

    /**
     *
     * @return Wild_emmer not WE
     */
    private List<String> getGroupD(){
        List<String> groupDList=new ArrayList<>();
        switch (this){
            case Subcontinent:
                for (GroupBySubcontinent group : GroupBySubcontinent.values()){
                    if (group.equals(GroupBySubcontinent.WE)) continue;
                    if (group.equals(GroupBySubcontinent.DE)) continue;
                    if (group.equals(GroupBySubcontinent.FT)) continue;
                    groupDList.add(group.getGroup());
                }
                return groupDList;
            case Subspecies:
                for (GroupBySubspecies group : GroupBySubspecies.values()){
                    if (group.equals(GroupBySubspecies.WE)) continue;
                    if (group.equals(GroupBySubspecies.DE)) continue;
                    if (group.equals(GroupBySubspecies.FT)) continue;
                    groupDList.add(group.getGroup());
                }
                return groupDList;
        }
        return groupDList;
    }

    public List<String> getGroup(SubgenomeCombination subgenomeCombination){
        switch (subgenomeCombination){
            case AB:
                return this.getGroupAB();
            case D:
                return this.getGroupD();
        }
        return null;
    }
}
