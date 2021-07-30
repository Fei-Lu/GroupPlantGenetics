package xiaohan.utils;

import java.util.HashMap;

/**
 * @ author: yxh
 * @ created: 2021-07-26 : 1:47 PM
 */
public class FastqUtils {

    public static String getphred(int score) {
        int[] scores = new int[41];
        for (int i = 0; i < 41; i++) {
            scores[i] = i;
        }
        String[] phreds = {"!", "\"", "#", "$", "%", "&", " ", "(", ")", "*", "+", ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I"};
        HashMap<Integer, String> scoreQMap = new HashMap<>();
        for (int i = 0; i < 41; i++) {
            scoreQMap.put(scores[i], phreds[i]);
        }
        String phred = scoreQMap.get(score);
        return phred;
    }

    public static int getscore(String phred) {
        int[] scores = new int[41];
        for (int i = 0; i < 41; i++) {
            scores[i] = i;
        }
        String[] phreds = {"!", "\"", "#", "$", "%", "&", " ", "(", ")", "*", "+", ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I"};
        HashMap<String,Integer> scoreQMap = new HashMap<>();
        for (int i = 0; i < 41; i++) {
            scoreQMap.put(phreds[i],scores[i]);
        }
        int score = scoreQMap.get(phred);
        return score;
    }
}
