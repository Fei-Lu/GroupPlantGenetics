package daxing.temp.individual;

public enum IfIntrogression {

    NNN(Donor.NONE), YNN(Donor.WE), NYN(Donor.DE),
    NNY(Donor.FTT), YYN(Donor.NONE), NYY(Donor.NONE),
    YNY(Donor.NONE), YYY(Donor.NONE), Y(Donor.AT), N(Donor.NONE);

    Donor donor;
    IfIntrogression(Donor donor) {
        this.donor=donor;
    }

    public Donor getDonor() {
        return donor;
    }
}
