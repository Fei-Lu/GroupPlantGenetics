package daxing.v2.localAncestryInfer.demography;

public class Default {

    Epoch epoch;
    Deme deme;

    public Default(){

    }

    public Default(Epoch epoch, Deme deme){
        this.epoch=epoch;
        this.deme=deme;
    }

    public Deme getDeme() {
        return deme;
    }

    public Epoch getEpoch() {
        return epoch;
    }

    public void setDeme(Deme deme) {
        this.deme = deme;
    }

    public void setEpoch(Epoch epoch) {
        this.epoch = epoch;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("Default{");
        sb.append("epoch=").append(epoch);
        sb.append(", deme=").append(deme);
        sb.append('}');
        return sb.toString();
    }
}
