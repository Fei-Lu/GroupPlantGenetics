package daxing.v2.localAncestryInfer.simulation;

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
