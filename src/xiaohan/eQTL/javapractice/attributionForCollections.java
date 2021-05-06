package xiaohan.eQTL.javapractice;

public class attributionForCollections {

    public attributionForCollections() {
        this.arraylength();
    }

    public void arraylength(){
        int[][] Array1 = new int[12][5];
        System.out.println(Array1.length);
        System.out.println(Array1[1].length);
    }

    public static void main(String[] args) {
        new attributionForCollections();
    }
}
