package daxing.v2.localAncestryInfer;

public class SolutionElement {

    WindowSource.Source source;
    int start; // inclusive
    int end; // exclusive

    public SolutionElement(WindowSource.Source source, int start, int end){
        this.source=source;
        this.start = start;
        this.end=end;
    }

    public void extend(){
        this.end=this.end+1;
    }
}
