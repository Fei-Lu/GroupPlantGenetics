package xiaohan.utils;

import com.ibm.icu.text.UForwardCharacterIterator;
import javafx.beans.property.IntegerProperty;

import java.util.HashMap;

public class chrUtils {

    public static String input = "1\t1A\t1\t471304005\n" +
            "3\t1B\t2\t438720154\n" +
            "5\t1D\t3\t452179604\n" +
            "7\t2A\t4\t462376173\n" +
            "9\t2B\t5\t453218924\n" +
            "11\t2D\t6\t462216879\n" +
            "13\t3A\t7\t454103970\n" +
            "15\t3B\t8\t448155269\n" +
            "17\t3D\t9\t476235359\n" +
            "19\t4A\t10\t452555092\n" +
            "21\t4B\t11\t451014251\n" +
            "23\t4D\t12\t451004620\n" +
            "25\t5A\t13\t453230519\n" +
            "27\t5B\t14\t451372872\n" +
            "29\t5D\t15\t451901030\n" +
            "31\t6A\t16\t452440856\n" +
            "33\t6B\t17\t452077197\n" +
            "35\t6D\t18\t450509124\n" +
            "37\t7A\t19\t450046986\n" +
            "39\t7B\t20\t453822637\n" +
            "41\t7D\t21\t453812268";

    public static String[][] map = new String[21][4];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////                  1,2,3,4,5,6, to 1A,1B,1D                                     ///////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * input: 1,2,3,4
     * output: 1A,1B
     */
    public static String getChr42toChrABD(String chr) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        int CHR = (Integer.parseInt(chr) + 1) / 2 - 1;
        String chrABD = map[CHR][1];
        return chrABD;
    }

    public static String getChr42toChrABD(int chr) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        int CHR = (chr + 1) / 2 - 1;
        String chrABD = map[CHR][1];
        return chrABD;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////                  1,2,3,4,5,6, to 1A,1B,1D,2A,2B,2D                            ///////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * input: 1,2,3,4
     * output: 1A,1B,1D,2A
     */
    public static String getChr21toChrABD(String chr) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        int CHR = Integer.parseInt(chr) - 1;
        String chrABD = map[CHR][1];
        return chrABD;
    }


    public static String getChr21toChrABD(int chr) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        int CHR = chr - 1;
        String chrABD = map[CHR][1];
        return chrABD;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////                    1A,1B,1D     to    1,2,3,4,5,6,                            ///////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * input: 1A,1B
     * output: 1,2,3,4
     */
    public static String getChrABDtoChr42(String chr, int position) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        HashMap<String, Integer> ABDindex = new HashMap<>();
        for (int i = 0; i < map.length; i++) {
            ABDindex.put(map[i][1], i);
        }
        int CHR = ABDindex.get(chr);
        if (Integer.parseInt(map[CHR][3]) < position) {
            CHR = Integer.parseInt(map[CHR][0])+1;
        } else {
            CHR = Integer.parseInt(map[CHR][0]);
        }
        return String.valueOf(CHR);
    }

    public static int getChrABDtoChr42Int(String chr, int position) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        HashMap<String, Integer> ABDindex = new HashMap<>();
        for (int i = 0; i < map.length; i++) {
            ABDindex.put(map[i][1], i);
        }
        int CHR = ABDindex.get(chr);
        if (Integer.parseInt(map[CHR][3]) < position) {
            CHR = Integer.parseInt(map[CHR][0])+1;
        } else {
            CHR = Integer.parseInt(map[CHR][0]);
        }
        return CHR;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////                    1A,1B,1D     to    1,2,3,                                  ///////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * input: 1A,1B
     * output: 1,2
     */
    public static String getChrABDtoChr21(String chr, int position) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        HashMap<String, Integer> ABDindex = new HashMap<>();
        for (int i = 0; i < map.length; i++) {
            ABDindex.put(map[i][1], i);
        }
        String CHR = map[ABDindex.get(chr)][2];
        return CHR;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////                    1A,1B,1D  pos   to    1,2,3,4,5,6,   pos                   ///////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * input: 1A,1B
     * output: 1,2,3,4
     */
    public static int getChrABDpostoChr42pos(String chr, int position) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        HashMap<String, Integer> ABDindex = new HashMap<>();
        for (int i = 0; i < map.length; i++) {
            ABDindex.put(map[i][1], i);
        }
        int CHR = ABDindex.get(chr);
        int out = 0;
        if (Integer.parseInt(map[CHR][3]) < position) {
            out = position - Integer.parseInt(map[CHR][3]);
        } else {
            out = position;
        }
        return out;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////                   1,2,3,4,5,6, pos   to   1A,1B,1D  pos                       ///////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public static int getChr42postoChrABDpos(String chr, int position) {
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                map[i][j] = input.split("\n")[i].split("\t")[j];
            }
        }
        int CHR = (Integer.parseInt(chr) + 1) / 2 - 1;
        if (Integer.parseInt(chr) % 2 == 0) {
            position += Integer.parseInt(map[CHR][3]);
        }
        return position;
    }
}
