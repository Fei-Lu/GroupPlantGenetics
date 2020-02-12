package daxing.common;

public class GeneWindow {

    String geneName=null;
    int start=Integer.MIN_VALUE;  //inclusive
    int end=Integer.MIN_VALUE;    //exclusive

    /**
     *
     * @param geneName
     * @param start row start index
     * @param end row end index
     */
    public GeneWindow(String geneName, int start, int end){
        this.geneName=geneName;
        this.start=start;
        this.end=end;
    }

    public String getGeneName() {
        return geneName;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getMid(){
        return (start+end-1)/2;
    }

    public int getRowNum(){
        return end-start;
    }

    public String getChr(){
        return this.geneName.substring(7,9);
    }
}
