package daxing.v2.localAncestryInfer.simulation;

/**
 * A deme always has at least one epoch.
 * The population parameters for a deme are allowed to change over time, but are fixed within an epoch.
 *
 * Why don’t epochs have a start_time too ?
 * The start time of an epoch can be inferred indirectly,
 * by looking at the deme’s start_time (for epoch 0),
 * or by looking at the previous epoch’s end_time (for epoch 1, 2, etc.).
 */
public class Epoch {

    int start_size;
    int end_time;

    public Epoch(){

    }

    public Epoch(int start_size, int end_time){
        this.start_size = start_size;
        this.end_time=end_time;
    }

    public int getStart_size() {
        return start_size;
    }

    public int getEnd_time() {
        return end_time;
    }

    public void setEnd_time(int end_time) {
        this.end_time = end_time;
    }

    public void setStart_size(int start_size) {
        this.start_size = start_size;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("Epoch{");
        sb.append("start_size=").append(start_size);
        sb.append(", end_time=").append(end_time);
        sb.append('}');
        return sb.toString();
    }
}