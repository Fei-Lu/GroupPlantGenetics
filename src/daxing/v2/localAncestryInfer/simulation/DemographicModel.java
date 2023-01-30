package daxing.v2.localAncestryInfer.simulation;

import java.util.ArrayList;
import java.util.List;

public class DemographicModel {

    String description;
    String time_units = "generations"; // generations or years

    double generation_time = 1; // 25 for human, 1 for wheat
    List<String> doi;

    List<Deme> demes;

    List<Migration> migrations;

    List<Pulse> pulses;

    List<Default> defaults;

    /**
     * A value other than 0 will cause a bug in Msprime
     */
    double selfing_rate = 0;

    /**
     * A value other than 0 will cause a bug in Msprime
     */
    double cloning_rate = 0;

    public DemographicModel(){

    }

    public DemographicModel(String description, String time_units, double generation_time, List<String> doi,
                            List<Default> defaults, List<Deme> demes, List<Migration> migrations, List<Pulse> pulses,
                            double selfing_rate, double cloning_rate){
        this.description=description;
        this.time_units=time_units;
        this.generation_time=generation_time;
        this.doi=doi;
        this.defaults=defaults;
        this.demes=demes;
        this.migrations=migrations;
        this.pulses=pulses;
        this.selfing_rate=selfing_rate;
        this.cloning_rate=cloning_rate;
    }

    private DemographicModel(Builder builder){
        this.time_units=builder.time_units;
        this.demes=builder.demes;
        this.description=builder.description;
        this.generation_time=builder.generation_time;
        this.doi=builder.doi;
        this.migrations=builder.migrations;
        this.pulses=builder.pulses;
        this.defaults=builder.defaults;
        this.selfing_rate=builder.selfing_rate;
        this.cloning_rate=builder.cloning_rate;
    }

    public static class Builder{

        /**
         * Required parameters
         */
        private String time_units;
        private List<Deme> demes;


        /**
         * Optional parameters
         */
        private String description;

        private double generation_time;
        private List<String> doi;
        private List<Migration> migrations;
        private List<Pulse> pulses;
        private List<Default> defaults;
        private double selfing_rate = 0;
        private double cloning_rate = 0;

        public Builder(String time_units, List<Deme> demes){
            this.time_units= time_units;
            this.demes=demes;
        }
        public Builder description(String description){
            this.description=description;
            return this;
        }

        public Builder generation_time(double generation_time){
            this.generation_time=generation_time;
            return this;
        }

        public Builder doi(List<String> doi){
            this.doi=doi;
            return this;
        }

        public Builder migrations(List<Migration> migrations){
            this.migrations=migrations;
            return this;
        }

        public Builder pulses(List<Pulse> pulses){
            this.pulses=pulses;
            return this;
        }

        public Builder defaults(List<Default> defaults){
            this.defaults=defaults;
            return this;
        }

        public Builder selfing_rate(double selfing_rate){
            this.selfing_rate=selfing_rate;
            return this;
        }

        public Builder cloning_rate(double cloning_rate){
            this.cloning_rate=cloning_rate;
            return this;
        }

        public DemographicModel build(){
            return new DemographicModel(this);
        }
    }

    /**
     * When you read a model into DemographicModel object, trim_default() need to be run.
     * Because some value may be missing in your model, and these value will be infer reasonably
     */
    public void trim_default(){
        List<Double> proportion = new ArrayList<>();
        proportion.add(1.0);
        String ancestor;
        int demeIndex;
        for (int i = 1; i < this.demes.size(); i++) {
            if (this.demes.get(i).ancestors.size()==1){
                this.demes.get(i).proportions=proportion;
                ancestor = this.demes.get(i).ancestors.get(0);
                demeIndex = this.getDemeIndex(ancestor);
                assert demeIndex >=0 : "error, check deme name";
                int epochCount = this.demes.get(demeIndex).epochs.size();
                this.demes.get(i).start_time = this.demes.get(demeIndex).epochs.get(epochCount-1).end_time;
            }

        }
    }

    public int getDemeIndex(String demeName){
        for (int i = 0; i < this.demes.size(); i++) {
            if (this.demes.get(i).name.equals(demeName)){
                return i;
            }
        }
        return -1;
    }

    public double getCloning_rate() {
        return cloning_rate;
    }

    public double getSelfing_rate() {
        return selfing_rate;
    }

    public double getGeneration_time() {
        return generation_time;
    }

    public List<Default> getDefaults() {
        return defaults;
    }

    public List<Deme> getDemes() {
        return demes;
    }

    public List<Migration> getMigrations() {
        return migrations;
    }

    public List<Pulse> getPulses() {
        return pulses;
    }

    public List<String> getDoi() {
        return doi;
    }

    public String getTime_units() {
        return time_units;
    }

    public String getDescription() {
        return description;
    }

    public void setCloning_rate(double cloning_rate) {
        this.cloning_rate = cloning_rate;
    }

    public void setDefaults(List<Default> defaults) {
        this.defaults = defaults;
    }

    public void setDemes(List<Deme> demes) {
        this.demes = demes;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public void setDoi(List<String> doi) {
        this.doi = doi;
    }

    public void setGeneration_time(double generation_time) {
        this.generation_time = generation_time;
    }

    public void setMigrations(List<Migration> migrations) {
        this.migrations = migrations;
    }

    public void setPulses(List<Pulse> pulses) {
        this.pulses = pulses;
    }

    public void setSelfing_rate(double selfing_rate) {
        this.selfing_rate = selfing_rate;
    }

    public void setTime_units(String time_units) {
        this.time_units = time_units;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("DemographicModel{");
        sb.append("description='").append(description).append('\'');
        sb.append(", time_units='").append(time_units).append('\'');
        sb.append(", generation_time=").append(generation_time);
        sb.append(", doi=").append(doi);
        sb.append(", demes=").append(demes);
        sb.append(", migrations=").append(migrations);
        sb.append(", pulses=").append(pulses);
        sb.append(", defaults=").append(defaults);
        sb.append(", selfing_rate=").append(selfing_rate);
        sb.append(", cloning_rate=").append(cloning_rate);
        sb.append('}');
        return sb.toString();
    }
}
