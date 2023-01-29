package daxing.v2.localAncestryInfer.simulation;

import it.unimi.dsi.fastutil.doubles.DoubleList;

import java.util.List;

/**
 * To model migration that is limited to a very short period of time, we can define one or more pulses.
 * A pulse has a proportions list property and a sources list property (analogous to the proportions and ancestors
 * properties of a deme). Each pulse proportion defines the proportion of the dest deme that is made up of ancestry
 * from the corresponding source deme at the instant after the pulse’s time.
 *
 * The exact duration of a pulse is not defined by the Demes specification.
 * Software which implements a continuous-time model (such as a coalescent simulator) might treat a pulse as
 * occurring instantaneously. In contrast, software which implements a discrete-time model is free to treat the
 * pulse as occurring over a single time step (such as a single generation).
 */
public class Pulse {

    /**
     * When multiple pulses are specified with the same time, we recommend using multiple sources with the desired
     * final ancestry proportions
     */
    List<String> sources;
    String dest;
    DoubleList proportions;
    double time;

    public Pulse(){

    }

    public Pulse(List<String> sources, String dest, DoubleList proportions, double time){
        this.sources=sources;
        this.dest=dest;
        this.proportions=proportions;
        this.time=time;
    }

    public DoubleList getProportions() {
        return proportions;
    }

    public double getTime() {
        return time;
    }

    public String getDest() {
        return dest;
    }

    public List<String> getSources() {
        return sources;
    }

    public void setProportions(DoubleList proportions) {
        this.proportions = proportions;
    }

    public void setDest(String dest) {
        this.dest = dest;
    }

    public void setSources(List<String> sources) {
        this.sources = sources;
    }

    public void setTime(double time) {
        this.time = time;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("Pulse{");
        sb.append("sources=").append(sources);
        sb.append(", dest='").append(dest).append('\'');
        sb.append(", proportions=").append(proportions);
        sb.append(", time=").append(time);
        sb.append('}');
        return sb.toString();
    }
}
