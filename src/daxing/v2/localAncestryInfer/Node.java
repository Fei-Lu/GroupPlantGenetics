package daxing.v2.localAncestryInfer;

import java.util.ArrayList;
import java.util.List;

public class Node {

    private List<Node> children = new ArrayList<>();
    private Node parent= null;

    private int position;
    private int extend=1;
    private int haplotypeIndex;


    public Node(int position, int haplotypeIndex){
        this.position = position;
        this.haplotypeIndex = haplotypeIndex;
    }

    public Node(int position, int haplotypeIndex, Node parent){
        this.position = position;
        this.haplotypeIndex = haplotypeIndex;
        this.setParent(parent);
    }

    public void addChild(Node child){
        // set parent before add child
       child.setParent(this);
       this.children.add(child);
    }

    public void setParent(Node parent){
       this.parent=parent;
    }


    public void addChild(int position, int haplotypeIndex){
        Node childNode = new Node(position, haplotypeIndex);
        this.addChild(childNode);
    }

    public void addChildren(List<Node> children){
        for (Node node:children){
            node.setParent(this);
        }
        this.children.addAll(children);
    }

    public void extend(){
        this.extend1();
    }

    private void extend1(){
        this.extend = extend + 1;
    }


}
