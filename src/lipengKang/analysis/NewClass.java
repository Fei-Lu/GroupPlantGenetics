/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.RandomArray;
import static utils.RandomArray.randomCommon;

/**
 *
 * @author kanglipeng
 */
public class NewClass {

    public static void main(String args[]) {
        /* String wheatA="AAAAAAAA";
  StringBuilder sorghum=new StringBuilder();
for(int i=0;i<wheatA.length();i++){sorghum.append('-');}
  System.out.println(sorghum.toString());*/

     /*   int x = 0;
        int[] reult1 = {20, 21, 2};
        for (int i : reult1) {

            if (i == 2) {
                x = 2;
                continue;
            }
            System.out.println(i);

        }*/
     StringBuilder test=new StringBuilder();
      int[] randomPos = randomCommon(0,7326564 , 10);
  for (int i : randomPos) {
  test.append(i+" ");
  
  }
         System.out.println( test.toString()); 



    }
}
