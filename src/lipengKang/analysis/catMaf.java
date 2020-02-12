/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import pgl.infra.utils.IOUtils;
import lipengKang.analysis.KStringUtils;

/**
 *
 * @author kanglipeng
 */
public class catMaf {

    public catMaf(String mafDir) {
        this.catSpecificMaf(mafDir);
    }

    public void catSpecificMaf(String mafDir) {

        try {
            BufferedReader br = IOUtils.getTextReader(mafDir);
            String temp = null;
            String[] tem = null;
            long Num = 0;
            while ((temp = br.readLine()) != null) {

                List<String> tList = KStringUtils.fastSplit(temp);
                List<String> tListNew = new ArrayList<>();

                if (temp.startsWith("s wheatA.chr")) {

                    for (int j = 0; j < tList.size(); j++) {

                        if (tList.get(j) != null && !tList.get(j).equals("")) {

                            tListNew.add(tList.get(j));
                        }
                    }
                    tem = tListNew.toArray(new String[tListNew.size()]);
                    if (tem == null || tem.length == 0) {
                        continue;
                    }
                    int alignedSeqSize = Integer.parseInt(tem[3]);
                    Num = Num + alignedSeqSize;
                }
            }
            System.out.println(Num);
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static void main(String[] args) {
        String mafDir = "/Users/kanglipeng/Desktop/gerptest/25.maf";
        new catMaf(mafDir);

    }

}
