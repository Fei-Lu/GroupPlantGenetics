package daxing.v2.localAncestryInfer;

import java.util.LinkedList;
import java.util.Queue;

public class BinaryTree {

    Node2 root;

    private Node2 addRecursive(Node2 current, int value) {
        if (current == null) {
            return new Node2(value);
        }

        if (value < current.value) {
            current.left = addRecursive(current.left, value);
        } else if (value > current.value) {
            current.right = addRecursive(current.right, value);
        } else {
            // value already exists
            return current;
        }

        return current;
    }

    public void add(int value) {
        root = addRecursive(root, value);
    }

    private boolean containsNodeRecursive(Node2 current, int value) {
        if (current == null) {
            return false;
        }
        if (value == current.value) {
            return true;
        }
        return value < current.value ? containsNodeRecursive(current.left, value) : containsNodeRecursive(current.right, value);
    }

    public boolean containsNode(int value) {
        return containsNodeRecursive(root, value);
    }

    private Node2 deleteRecursive(Node2 current, int value) {
        if (current == null) {
            return null;
        }

        if (value == current.value) {

            // a node has no children – this is the simplest case;
            // we just need to replace this node with null in its parent node
            if (current.left == null && current.right == null) {
                return null;
            }

            // a node has exactly one child –
            // in the parent node, we replace this node with its only child.
            if (current.right == null) {
                return current.left;
            }

            if (current.left == null) {
                return current.right;
            }

            // Finally, we have to handle the case where the node has two children.
            // First, we need to find the node that will replace the deleted node.
            // We'll use the smallest node of the soon to be deleted node's right sub-tree:
            // Then we assign the smallest value to the node to delete, and after that, we'll delete it from the right sub-tree:
            int smallestValue = findSmallestValue(current.right);
            current.value = smallestValue;
            current.right = deleteRecursive(current.right, smallestValue);
            return current;

        }
        if (value < current.value) {
            current.left = deleteRecursive(current.left, value);
            return current;
        }
        current.right = deleteRecursive(current.right, value);
        return current;
    }

    private int findSmallestValue(Node2 root) {
        return root.left == null ? root.value : findSmallestValue(root.left);
    }

    public void delete(int value) {
        root = deleteRecursive(root, value);
    }

    /**
     * Depth-First Search, in-order
     * Depth-first search is a type of traversal that goes deep as much as possible in every child before exploring the next sibling.
     * There are several ways to perform a depth-first search: in-order, pre-order and post-order.
     * The in-order traversal consists of first visiting the left sub-tree, then the root node, and finally the right sub-tree:
     * @param node
     */
    public void traverseInOrder(Node2 node) {
        if (node != null) {
            traverseInOrder(node.left);
            System.out.print(" " + node.value);
            traverseInOrder(node.right);
        }
    }

    /**
     * Depth-First Search, pre-order
     * Depth-first search is a type of traversal that goes deep as much as possible in every child before exploring the next sibling.
     * There are several ways to perform a depth-first search: in-order, pre-order and post-order.
     * The in-order traversal consists of first visiting the left sub-tree, then the root node, and finally the right sub-tree:
     * @param node
     */
    public void traversePreOrder(Node2 node) {
        if (node != null) {
            System.out.print(" " + node.value);
            traversePreOrder(node.left);
            traversePreOrder(node.right);
        }
    }

    /**
     * Depth-First Search, post-order
     * Depth-first search is a type of traversal that goes deep as much as possible in every child before exploring the next sibling.
     * There are several ways to perform a depth-first search: in-order, pre-order and post-order.
     * The in-order traversal consists of first visiting the left sub-tree, then the root node, and finally the right sub-tree:
     * @param node
     */
    public void traversePostOrder(Node2 node) {
        if (node != null) {
            traversePostOrder(node.left);
            traversePostOrder(node.right);
            System.out.print(" " + node.value);
        }
    }

    /**
     * Breadth-First Search
     * This is another common type of traversal that visits all the nodes of a level before going to the next level.
     * This kind of traversal is also called level-order, and visits all the levels of the tree starting from the root, and from left to right.
     * For the implementation, we'll use a Queue to hold the nodes from each level in order. We'll extract each node
     * from the list, print its values, then add its children to the queue.
     */
    public void traverseLevelOrder() {
        if (root == null) {
            return;
        }

        Queue<Node2> nodes = new LinkedList<>();
        nodes.add(root);

        while (!nodes.isEmpty()) {

            Node2 node = nodes.remove();

            System.out.print(" " + node.value);

            if (node.left != null) {
                nodes.add(node.left);
            }

            if (node.right != null) {
                nodes.add(node.right);
            }
        }
    }


}
