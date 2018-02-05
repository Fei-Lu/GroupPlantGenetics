/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

/**
 * Hold a range or interval on chromosome
 * Chromosome and position generally should be greater than 0
 * Note position of Integer.Min_Value is not allowed, since this value is used to fast locate range index of chromosome
 * @author Fei Lu
 */
public class Range implements Comparable<Range>{
    public int chr;
    public int start;
    /**Exclusive*/
    public int end;
    
    public Range (int chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }
    
    /**
     * Return chromosome of this range
     * @return 
     */
    public int getRangeChromosome () {
        return this.chr;
    }
    
    /**
     * Return start point of this range, inclusive
     * @return 
     */
    public int getRangeStart () {
        return this.start;
    }
    
    /**
     * Return end point of this range, exclusive
     * @return 
     */
    public int getRangeEnd () {
        return this.end;
    }
    
    public int getRangeSize () {
        return this.end - this.start;
    }
    
    public void setRangeChromosome (int chr) {
        this.chr = chr;
    }
    
    public void setRangeStart (int start) {
        this.start = start;
    }
    
    public void setRangeEnd (int end) {
        this.end = end;
    }
    
    /**
     * Return a string of all members
     * @return 
     */
    public String getOutputString () {
       StringBuilder sb = new StringBuilder();
       sb.append(chr).append("\t").append(start).append("\t").append(end);
       return sb.toString();
    }
    
    /**
     * Return if this range has overlap with another range
     * @param other
     * @return 
     */
    public boolean isOverlap (Range other) {
        if (this.chr != other.chr) return false;
        if (this.end <= other.start || other.end <= this.start) return false;
        return true;
    }
    
    /**
     * Return if this range is within another range
     * @param other
     * @return 
     */
    public boolean isWithin (Range other) {
        if (this.chr != other.chr) return false;
        if (this.start >= other.start && this.end <= other.end) return true;
        return false;
    }
    
    /**
     * Return if this range contains another range
     * @param other
     * @return 
     */
    public boolean isContain (Range other) {
        if (this.chr != other.chr) return false;
        if (this.start <= other.start && this.end >= other.end) return true;
        return false;
    }
    
    /**
     * Return an intersection range between this range and another range
     * Return null if there is not any intersection
     * @param other
     * @return 
     */
    public Range intersection (Range other) {
        if (!this.isOverlap(other)) return null;
        if (this.start < other.start) return new Range(this.chr, other.start, this.end < other.end? this.end : other.end);
        return new Range(this.chr, this.start, this.end < other.end? this.end : other.end);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Range other = (Range) obj;
        if (this.chr != other.chr) {
            return false;
        }
        if (this.start != other.start) {
            return false;
        }
        if (this.end != other.end) {
            return false;
        }
        return true;
    }
    
    public int hashCode() {
        return start;
    }
    
    @Override
    public int compareTo(Range t) {
        if (chr == t.chr) {
            if (start == t.start) {
                if (end == t.end) return 0;
                else if (end < t.end) return -1;
                return 1;
            }
            else if (start < t.start) return -1;
            return 1;
        }
        else if (chr < t.chr) return -1;
        return 1;
    }
}
