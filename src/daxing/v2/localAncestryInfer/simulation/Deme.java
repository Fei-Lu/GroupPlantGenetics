package daxing.v2.localAncestryInfer.simulation;

import it.unimi.dsi.fastutil.doubles.DoubleList;
import java.util.List;

/**
 * A deme is a collection of individuals that share a set of population parameters at any given time.
 */
public class Deme {

    /**
     * name of this Deme
     */
    String name;

    /**
     * When a deme has multiple ancestors, these appear in the ancestors list as one might expect.
     * But for multiple ancestors we need to also specify the proportion of ancestry inherited from each ancestor.
     * This is done using the deme’s proportions list property.
     */
    List<String> ancestors;
    DoubleList proportions;

    /**
     * When multiple ancestors exists, there are two types of events to consider,
     * namely branch event and split event.
     *
     * An alternative way of modelling a population split is for the ancestral deme to remain alive after the split.
     * We will refer to this as a branch event, rather than a split event. Because the descendant deme's ancestor
     * exists until 0 generations ago, so we must explicitly provide a start_time for it.
     *
     * When involved in split event, ancestral deme are not allowed alive after the split.
     * When a deme has an ancestor, its start_time does not default to .inf.
     * In this case, the start_time for the deme is inherited from the end_time of the ancestor
     *
     * With multiple ancestors, the start_time of the descendant deme
     * does not default to the end_time of any of its ancestors.
     * So the start_time must always be specified for a deme with multiple ancestors.
     */
    int start_time;

    /**
     * We partition a deme’s interval of existence into distinct epochs.
     */
    List<Epoch> epochs;

    public Deme(){

    }

    public Deme(String name, List<String> ancestors, DoubleList proportions, int start_time, List<Epoch> epochs){
        this.name=name;
        this.ancestors=ancestors;
        this.proportions=proportions;
        this.start_time=start_time;
        this.epochs=epochs;
    }

    public DoubleList getProportions() {
        return proportions;
    }

    public int getStart_time() {
        return start_time;
    }

    public List<Epoch> getEpochs() {
        return epochs;
    }

    public List<String> getAncestors() {
        return ancestors;
    }

    public String getName() {
        return name;
    }

    public void setAncestors(List<String> ancestors) {
        this.ancestors = ancestors;
    }

    public void setEpochs(List<Epoch> epochs) {
        this.epochs = epochs;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setProportions(DoubleList proportions) {
        this.proportions = proportions;
    }

    public void setStart_time(int start_time) {
        this.start_time = start_time;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("Deme{");
        sb.append("name='").append(name).append('\'');
        sb.append(", ancestors=").append(ancestors);
        sb.append(", proportions=").append(proportions);
        sb.append(", start_time=").append(start_time);
        sb.append(", epochs=").append(epochs);
        sb.append('}');
        return sb.toString();
    }
}
