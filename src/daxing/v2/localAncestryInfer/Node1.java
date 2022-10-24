package daxing.v2.localAncestryInfer;

import java.util.ArrayList;
import java.util.List;

public class Node1 {

    List<Node1> switchChildren;
    int position;
    int extend;
    int haplotypeIndex;

    Node1(int position, int extend, int haplotypeIndex){
        this.position=position;
        this.extend=extend;
        this.haplotypeIndex=haplotypeIndex;
        this.switchChildren = new ArrayList<>();
    }

}
