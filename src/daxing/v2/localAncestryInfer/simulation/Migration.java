package daxing.v2.localAncestryInfer.simulation;

import java.util.List;

public class Migration {

    /**
     * When source and dest are provided, this is a Asymmetric migration
     */
    String source;
    String dest;

    /**
     * Migration rate units are always “per generation” regardless of the chosen value for time_units.
     */
    double rate;


    /**
     * If the source and dest demes do not have identical start or end times,
     * then by default the migration will occur over the period of time when both demes exist simultaneously.
     *
     * To obtain greater control over when migrations occur, we can use the start_time and end_time migration
     * properties.
     */
    int start_time;
    int end_time;

    /**
     * It’s common to model migrants in both directions simultaneously, with the same rate, namely Asymmetric migration
     * The demes property can list arbitrarily many demes.
     * Symmetric migration will occur between all pairwise combinations demes
     */
    List<String> demes;

    public Migration(){}

    public Migration(String source, String dest, double rate, int start_time, int end_time, List<String> demes){
        this.source=source;
        this.dest=dest;
        this.rate=rate;
        this.start_time=start_time;
        this.end_time=end_time;
        this.demes=demes;
    }

    public int getStart_time() {
        return start_time;
    }

    public int getEnd_time() {
        return end_time;
    }

    public List<String> getDemes() {
        return demes;
    }

    public double getRate() {
        return rate;
    }

    public String getDest() {
        return dest;
    }

    public String getSource() {
        return source;
    }

    public void setEnd_time(int end_time) {
        this.end_time = end_time;
    }

    public void setStart_time(int start_time) {
        this.start_time = start_time;
    }

    public void setDemes(List<String> demes) {
        this.demes = demes;
    }

    public void setDest(String dest) {
        this.dest = dest;
    }

    public void setRate(double rate) {
        this.rate = rate;
    }

    public void setSource(String source) {
        this.source = source;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("Migration{");
        sb.append("source='").append(source).append('\'');
        sb.append(", dest='").append(dest).append('\'');
        sb.append(", rate=").append(rate);
        sb.append(", start_time=").append(start_time);
        sb.append(", end_time=").append(end_time);
        sb.append(", demes=").append(demes);
        sb.append('}');
        return sb.toString();
    }
}
