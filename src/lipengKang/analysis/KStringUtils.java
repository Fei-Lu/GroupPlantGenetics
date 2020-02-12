/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import pgl.infra.utils.PStringUtils;
import java.util.List;
import static pgl.infra.utils.PStringUtils.fastSplit;

/**
 *
 * @author kanglipeng
 */
public class KStringUtils extends PStringUtils {
     /**
     * Return a list of split String using Guava splitter, blank space
     * @param line
     * @return 
     */
    public static List<String> fastSplit (String line) {
        List<String> ls = fastSplit(line, " ");
        return ls;
    }
    
     /**
     * Return a list of split String using Guava splitter, comma
     * @param line
     * @return 
     */
       public static List<String> fastSplitComma (String line) {
        List<String> ls = fastSplit(line, ",");
        return ls;
    } 
       /**
     * Return a list of split String using Guava splitter, strigula
     * @param line
     * @return 
     */
    
      public static List<String> fastSplitStrigula (String line) {
        List<String> ls = fastSplit(line, "-");
        return ls;
      }
       /**
     * Return a list of split String using Guava splitter, dot
     * @param line
     * @return 
     */
    
      public static List<String> fastSplitDot (String line) {
        List<String> ls = fastSplit(line, ".");
        return ls;
      }
      /**
     * Return a list of split String using Guava splitter, dot
     * @param line
     * @return 
     */
    
      public static List<String> fastSplitSemicolon (String line) {
        List<String> ls = fastSplit(line, ";");
        return ls;
      }
}
