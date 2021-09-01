package daxing.load.esg;

public enum ESG {
    abd(0, "000"), abD(1, "001"), aBd(2, "010"), Abd(3, "100"),
    aBD(4, "011"), AbD(5, "101"), ABd(6, "110"), ABD(7, "111");

    int value;
    String character;

    ESG(int value, String character) {
        this.value=value;
        this.character = character;
    }

    public int getValue() {
        return value;
    }

    public String getCharacter() {
        return character;
    }

    public static ESG newInstanceFrom(String esgInfo){
        switch (esgInfo){
            case "000":
                return ESG.abd;
            case "001":
                return ESG.abD;
            case "010":
                return ESG.aBd;
            case "100":
                return ESG.Abd;
            case "011":
                return ESG.aBD;
            case "101":
                return ESG.AbD;
            case "110":
                return ESG.ABd;
            case "111":
                return ESG.ABD;
        }
        return null;
    }
}
