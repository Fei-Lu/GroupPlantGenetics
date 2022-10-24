package daxing.v2.localAncestryInfer;

import java.util.Objects;

public class SolutionElement implements Cloneable {

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

    public SolutionElement clone(){
        SolutionElement clone = null;
        try {
            clone = (SolutionElement) super.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
        }
        return clone;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public WindowSource.Source getSource() {
        return source;
    }

    @Override
    public boolean equals(Object ob){
        if (ob == this) return true;
        if (!(ob instanceof SolutionElement)) return false;
        SolutionElement other = (SolutionElement) ob;
        boolean currencyCodeEquals = (this.source == null && other.source == null)
                || (this.source != null && this.source.equals(other.source));

        return this.start == other.start && currencyCodeEquals && this.end == other.end;
    }

    @Override
    public int hashCode() {
        return Objects.hash(source, start, end);
    }
}
