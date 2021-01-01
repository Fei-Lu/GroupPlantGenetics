package jijin;

public class SplitReads {
    String tags =null;
    String srName = null;
    int chrom = 0 ;
    int start=0;
    int end =0;


    //byte strand = Byte.MIN_VALUE;
    //String biotype = null;
    //String description = null;

    //int longestTranscriptIndex = -1;


    public  SplitReads (String tags, String srName, int chrom, int start, int end) {
        this.tags = tags;
        this.srName = srName;
        this.chrom = chrom;
        this.start = start;
        this.end=end;
    }

    public String gettags(){return this.tags;}

    public String getsineme(){return this.srName;}

    public int getchrom(){return this.chrom;}

    public int getstart(){
        return  this.start;
    }

    public int getend(){
        return this.end;
    }

}
