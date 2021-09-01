package daxing.load.esg;

public enum IfPseudoHexaploid {

    HEXAPLOID(0), PSEUDO_HEXAPLOID(1);

    int value;

    IfPseudoHexaploid(int value) {
        this.value=value;
    }

    public int getValue() {
        return value;
    }
}
