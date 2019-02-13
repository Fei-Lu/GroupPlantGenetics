/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import cern.colt.Arrays;
import static java.lang.Math.random;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author kanglipeng
 */
public class HashmapPutTest {

    public static void main(String[] args) {
        HashMap<String, Set<String>> ma = new HashMap();
        String[] a = {"A", "B", "C", "D", "E", "F", "G", "H"};
        String[] b = {"a", "b", "c", "d", "a", "b", "c", "d"};
        Set<String> s = new HashSet();
        for (int i = 0; i < 8; i++) {
            s = new HashSet();            
            s = ma.get(b[i]);
            s.add(a[i]);
        }
        int aaa = 100;
        System.out.print(ma.get("a"));

    }

}
