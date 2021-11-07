package daxing.common.maf;

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

public enum Species {

    traes(0, "Traes"), thelo(1, "Thelo"), secer(2, "Secer"),
    hovul(3, "Hovul"), avatl(4, "Avatl"), brdis(5, "Brdis"),
    phedu(6, "Phedu"), orsat(7, "Orsat"), ercur(8, "Ercur"),
    ortho(9, "Ortho"), sobic(10, "Sobic"), coaqu(11, "Coaqu"),
    zemayB73(12, "ZemayB73"), pahal(13, "Pahal"), seita(14, "Seita"),
    ancom(15, "Ancom");

    int index;
    String str;
    Species(int index, String str){
        this.index=index;
        this.str=str;
    }

    public int getIndex() {
        return index;
    }

    public String getStr() {
        return str;
    }

    private static final Map<String,Species> strToSpeciesMap= Arrays.stream(values()).collect(Collectors.toMap(Species::getStr, e->e));

    public static Species getInstanceFromStr(Species str){
        return strToSpeciesMap.get(str);
    }
}
