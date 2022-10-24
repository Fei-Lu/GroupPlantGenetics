package daxing.v2.localAncestryInfer;

import java.util.ArrayList;
import java.util.List;

public class Tree {

    Node1 root;

    private Node1 addRecursive(Node1 current, int position, int extend, int haplotypeIndex) {
        if (current == null) {
            return new Node1(position, extend, haplotypeIndex);
        }

        if (haplotypeIndex == current.haplotypeIndex){
            current.extend++;
            return current;
        }else {
            current.switchChildren.add(new Node1(position, extend, haplotypeIndex));
            return current;
        }
    }

    public void add(int position, int extend, int haplotypeIndex) {
        root = addRecursive(root, position, extend, haplotypeIndex);
    }
}
