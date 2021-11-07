package daxing.common.maf;

/**
 * Coordinate Transforms for multiple alignment format
 * ucscRevStart = chromSize - oneEnd
 * ucscRevEnd   = chromSize - (oneStart - 1)
 * [oneRevStart, oneRevEnd] = [ucscRevStart+1, ucscRevEnd] = [chromSize-oneEnd+1, chromSize-oneStart+1]
 */
public class Coordinate {

    int chromSize;

    /**
     * 1-based [oneStart, oneEnd]
     */
    int oneStart;
    int oneEnd;

    /**
     * 0-based [oneStart-1, oneEnd)
     */
    int zeroStart;
    int zeroEnd;

    /**
     * 0-based half-open [ucscRevStart, ucscRevEnd)
     */
    int ucscRevStart;
    int ucscRevEnd;

    /**
     * 1-based [oneRevStart, oneRevEnd]
     */
    int oneRevStart;
    int oneRevEnd;

    protected static Coordinate holder=null;

    Coordinate(){
        this.Initialize();
    }

    private void Initialize(){
        oneStart=Integer.MIN_VALUE;
        oneEnd=Integer.MIN_VALUE;
        zeroStart=Integer.MIN_VALUE;
        zeroEnd=Integer.MIN_VALUE;
        chromSize=Integer.MIN_VALUE;
        ucscRevStart=Integer.MIN_VALUE;
        ucscRevEnd=Integer.MIN_VALUE;
        oneRevStart=Integer.MIN_VALUE;
        oneRevEnd=Integer.MIN_VALUE;
    }

    public static Coordinate getInstance(){
        if (holder==null){
            holder= new Coordinate();
            return holder;
        }else {
            holder.Clear();
            return holder;
        }
    }

    public void Clear(){
        this.Initialize();
    }

    private void setFromForward(int zeroStart, int chromSize){
        this.chromSize= chromSize;
        this.zeroStart=zeroStart;
        this.zeroEnd=zeroStart+chromSize;
        this.oneStart=zeroStart+1;
        this.oneEnd=zeroStart+chromSize;
        this.ucscRevStart=chromSize-oneEnd;
        this.ucscRevEnd=chromSize-zeroStart;
        this.oneRevStart=chromSize-oneEnd+1;
        this.oneRevEnd=chromSize-oneStart+1;
    }

    private void setFromReverse(int ucscRevStart, int chromSize){
        this.chromSize=chromSize;
        this.ucscRevStart=ucscRevStart;
        this.ucscRevEnd=ucscRevStart+chromSize;
        this.oneRevStart=ucscRevStart+1;
        this.oneRevEnd=ucscRevStart+chromSize;;
        this.oneStart=chromSize+1-ucscRevEnd;
        this.oneEnd=chromSize-ucscRevStart;
        this.zeroStart=oneStart-1;
        this.zeroEnd=oneEnd;
    }

    public void setFrom(Strand strand, int ucscStart, int chromSize){
        switch (strand){
            case Forward:
                setFromForward(ucscStart, chromSize);
                break;
            case Reverse:
                setFromReverse(ucscStart, chromSize);
                break;
        }
    }

    public int getChromSize() {
        return chromSize;
    }

    public int getZeroStart() {
        return zeroStart;
    }

    public int getZeroEnd() {
        return zeroEnd;
    }

    public int getUcscRevStart() {
        return ucscRevStart;
    }

    public int getUcscRevEnd() {
        return ucscRevEnd;
    }

    public int getOneEnd() {
        return oneEnd;
    }

    public int getOneStart() {
        return oneStart;
    }

    public int getOneRevEnd() {
        return oneRevEnd;
    }

    public int getOneRevStart() {
        return oneRevStart;
    }
}
