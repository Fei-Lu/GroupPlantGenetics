package xiaohan.utils;

import java.util.HashMap;

public class chrUtils {

    /**
    * from chr1,2,3, to chr1A,chr1B,chr1D
     * input : String like 1,2,3,4,5,6
     * output : String like 1A,1B,1D
    * */
    public static String getChrtoChrABD(String chr){
        int chrnumber = Integer.parseInt(chr) / 6 + 1;
        String subgenome = null;
        if(Integer.parseInt(chr)%6 == 1 || Integer.parseInt(chr)%6 == 2){
            subgenome = "A";
        }
        if(Integer.parseInt(chr)%6 == 3 || Integer.parseInt(chr)%6 == 4){
            subgenome = "B";
        }
        if(Integer.parseInt(chr)%6 == 5 || Integer.parseInt(chr)%6 == 6){
            subgenome = "D";
        }
        String chrABD = String.valueOf(chrnumber) + subgenome;
        return chrABD;
    }

    /**
     * from chr1A,chr1B,chr1D to chr1,2,3
     * input : String like 1A,1B,1D & Interger of position
     * output : String like 1,2,3,4,5,6
     * */
    public static String getChrABDtoChr(String chrABD,int position){
        String length = "471304005\n" +
                "438720154\n" +
                "452179604\n" +
                "462376173\n" +
                "453218924\n" +
                "462216879\n" +
                "454103970\n" +
                "448155269\n" +
                "476235359\n" +
                "452555092\n" +
                "451014251\n" +
                "451004620\n" +
                "453230519\n" +
                "451372872\n" +
                "451901030\n" +
                "452440856\n" +
                "452077197\n" +
                "450509124\n" +
                "450046986\n" +
                "453822637\n" +
                "453812268";
        String[] lengths = length.split("\n");
        String subgenome = chrABD.substring(1,2);
        int chrnumber = Integer.parseInt(chrABD.substring(0,1));
        int chr42 = -1;
        if(subgenome.equals("A")){
            chr42 = (chrnumber-1)*6 + 1;
        }
        if(subgenome.equals("B")){
            chr42 = (chrnumber-1)*6 + 3;
        }
        if(subgenome.equals("D")){
            chr42 = (chrnumber-1)*6 + 5;
        }
        if(position > Integer.parseInt(lengths[chr42/2])){
            chr42 ++;
        }
        String chr = String.valueOf(chr42);
        return chr;
    }

    /**
     * from chrindex to chrABDindex
     * indexlist like : 1,2,3,4,5,6 of 1A,1B,1D,2A,2B,2D
     * input : 1,2,3,4,5,6
     * output : 1,1,2,2,3,3
    * */
    public static int getChrIndextoChrABDIndex(int chrindex){
        int chrABDindex = (chrindex -1 )/2 + 1;
        return chrABDindex;
    }

    /**
     * from chrABDindex to chrindex
     * indexlist like : 1,2,3,4,5,6 of 1A,1B,1D,2A,2B,2D
     * input : 1A,1B,1D & Interger of position
     * output : 1,2,3,4,5,6
     * */

    public static int getChrABDIndextoChrIndex(int chrABDindex,int position){
        int chrindex = 0;
        String length = "471304005\n" +
                "438720154\n" +
                "452179604\n" +
                "462376173\n" +
                "453218924\n" +
                "462216879\n" +
                "454103970\n" +
                "448155269\n" +
                "476235359\n" +
                "452555092\n" +
                "451014251\n" +
                "451004620\n" +
                "453230519\n" +
                "451372872\n" +
                "451901030\n" +
                "452440856\n" +
                "452077197\n" +
                "450509124\n" +
                "450046986\n" +
                "453822637\n" +
                "453812268";
        String[] lengths = length.split("\n");
        int chr42 = -1;
        chr42 = chrABDindex*2;
        if(position < Integer.parseInt(lengths[chrABDindex-1])){
            chr42 --;
        }
        return chr42;
    }

    public static String getNamechr(int chrindex){
        String chrName = "chr"+chrindex;
        return chrName;
    }

    public static String getNamechrABD(int chrABDindex){
        int chrnumber = -1;
        if(chrABDindex%3 == 0){
            chrnumber = chrABDindex/3;
        }
        else chrnumber = chrABDindex/3+1;
        String subgenome = "D";
        if(chrABDindex%3 == 1){
            subgenome = "A";
        }
        if(chrABDindex%3 == 2){
            subgenome = "B";
        }
        String chrABDName = "chr"+String.valueOf(chrnumber)+subgenome;
        return chrABDName;
    }

}
