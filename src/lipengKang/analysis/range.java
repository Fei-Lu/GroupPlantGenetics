/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author kanglipeng
 */
public class range {
   /**
     * Insert a {@link format.range.Range} into the list
     * @param <T> {@link format.range.Range}
     * @param rangeIndex
     * @param r 
     * @return true if it is successful, otherwise false
     */
    public static ArrayList getCommonElements(ArrayList <Integer> listA, ArrayList<Integer> listB) {
      ArrayList <Integer>smallestList = listA.size() < listB.size() ? listA : listB;
     ArrayList <Integer>biggestList = listA.size() > listB.size() ? listA : listB;
      Set<Integer> commonSet = new HashSet<>();
      commonSet.addAll(smallestList);
      Set<Integer> bigSet = new HashSet<>();
      bigSet.addAll(biggestList);
      commonSet.retainAll(bigSet);
      ArrayList<Integer> commonList = new ArrayList<>(commonSet);
      return commonList;
  }

    
}
