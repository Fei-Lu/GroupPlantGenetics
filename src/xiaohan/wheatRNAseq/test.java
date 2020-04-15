package xiaohan.wheatRNAseq;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

import pgl.infra.utils.IOUtils;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author yxh
 */
public class test {

    public test() throws IOException {
//        this.mkPosGeneMap();
        //this.wordList();
        this.writecode();


    }

    /**
     * @throws IOException
     */


    public void writecode(){
        for(int i = 18;i<45;i++){
            System.out.print("nohup vcftools --gzvcf /data3/wgs/vcf/GATK/vmap3/1.SNP/");
            System.out.print(i);
            System.out.print(".snp.vcf.gz --maf 0 --max-maf 0.05 --out ");
            System.out.print(i+".snp.maf005 --recode && tabix -p "+i+".snp.maf005.recode.vcf.gz >log1.txt 2>&1 &");
        }
    }

    public void mkPosGeneMap() {
        String geneNameS = null;
        int gfIndex = 0;
        ArrayList<String> geneNameList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/refer/wheat_v1.1_Lulab.gtf");
            String temp = null;
            Set<String> geneSet = new HashSet<String>();
            Set<Integer> chrSet = new HashSet<>();
            String[] tem = null;
            String geneName = null;
            HashMap<String, Integer> geneChrMap = new HashMap();
            HashMap<String, Integer> geneMinMap = new HashMap();
            HashMap<String, Integer> geneMaxMap = new HashMap();
            HashMap<String, Byte> geneStrandMap = new HashMap();
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                if (!tem[2].startsWith("exon")) continue;
                String[] te = tem[8].split(";");
                geneName = te[1].split("\"")[1];
                if (!geneSet.contains(geneName)) {
                    geneMinMap.put(geneName, Integer.MAX_VALUE);
                    geneMaxMap.put(geneName, Integer.MIN_VALUE);
                }
                geneSet.add(geneName);
                int min = Integer.parseInt(tem[3]);
                int max = Integer.parseInt(tem[4]);
                if (geneMinMap.get(geneName) > min) geneMinMap.put(geneName, min);
                if (geneMaxMap.get(geneName) < max) geneMaxMap.put(geneName, max);
                int chr = Integer.parseInt(tem[0]);
                chrSet.add(chr);
                geneChrMap.put(geneName, chr);
                if (tem[6].startsWith("-")) geneStrandMap.put(geneName, (byte) 1);
                else geneStrandMap.put(geneName, (byte) 0);
            }

            for (String s : geneSet) {
                geneNameList.add(s);
            }
            Collections.sort(geneNameList);
            for (String s : geneNameList) {
                System.out.println(geneChrMap.get(s) + "\t" + geneMinMap.get(s) + "\t" + geneMaxMap.get(s) + "\t" + s);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
            System.out.println(geneNameS + "\t" + gfIndex);
        }

    }

    public void findNumber() throws IOException {
        String infile1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/table/ZX-ZX.txt";
        HashMap<String, String> trans = new HashMap<String, String>();
        String table = "ZX393\tZX1497\n" +
                "ZX394\tZX1498\n" +
                "ZX395\tZX1499\n" +
                "ZX396\tZX1500\n" +
                "ZX397\tZX1501\n" +
                "ZX398\tZX1502\n" +
                "ZX399\tZX1503\n" +
                "ZX400\tZX1504\n" +
                "ZX401\tZX1505\n" +
                "ZX402\tZX1506\n" +
                "ZX403\tZX1507\n" +
                "ZX404\tZX1508\n" +
                "ZX405\tZX1509\n" +
                "ZX406\tZX1510\n" +
                "ZX407\tZX1511\n" +
                "ZX408\tZX1512\n" +
                "ZX409\tZX1513\n" +
                "ZX410\tZX1514\n" +
                "ZX411\tZX1515\n" +
                "ZX412\tZX1516\n" +
                "ZX413\tZX1517\n" +
                "ZX414\tZX1518\n" +
                "ZX415\tZX1519\n" +
                "ZX416\tZX1520\n" +
                "ZX417\tZX1521\n" +
                "ZX418\tZX1522\n" +
                "ZX419\tZX1523\n" +
                "ZX420\tZX1524\n" +
                "ZX421\tZX1525\n" +
                "ZX422\tZX1526\n" +
                "ZX423\tZX1527\n" +
                "ZX424\tZX1528\n" +
                "ZX425\tZX1529\n" +
                "ZX426\tZX1530\n" +
                "ZX427\tZX1531\n" +
                "ZX428\tZX1532\n" +
                "ZX429\tZX1533\n" +
                "ZX430\tZX1534\n" +
                "ZX431\tZX1535\n" +
                "ZX432\tZX1536\n" +
                "ZX433\tZX1537\n" +
                "ZX434\tZX1538\n" +
                "ZX435\tZX1539\n" +
                "ZX436\tZX1540\n" +
                "ZX437\tZX1541\n" +
                "ZX438\tZX1542\n" +
                "ZX439\tZX1543\n" +
                "ZX440\tZX1544\n" +
                "ZX441\tZX1545\n" +
                "ZX442\tZX1546\n" +
                "ZX443\tZX1547\n" +
                "ZX444\tZX1548\n" +
                "ZX445\tZX1549\n" +
                "ZX446\tZX1550\n" +
                "ZX447\tZX1551\n" +
                "ZX448\tZX1552\n" +
                "ZX449\tZX1553\n" +
                "ZX450\tZX1554\n" +
                "ZX451\tZX1555\n" +
                "ZX452\tZX1556\n" +
                "ZX453\tZX1557\n" +
                "ZX454\tZX1558\n" +
                "ZX455\tZX1559\n" +
                "ZX456\tZX1560\n" +
                "ZX457\tZX1561\n" +
                "ZX458\tZX1562\n" +
                "ZX459\tZX1563\n" +
                "ZX460\tZX1564\n" +
                "ZX461\tZX1565\n" +
                "ZX462\tZX1566\n" +
                "ZX463\tZX1567\n" +
                "ZX464\tZX1568\n" +
                "ZX465\tZX1569\n" +
                "ZX466\tZX1570\n" +
                "ZX467\tZX1571\n" +
                "ZX468\tZX1572\n" +
                "ZX469\tZX1573\n" +
                "ZX470\tZX1574\n" +
                "ZX471\tZX1575\n" +
                "ZX472\tZX1576\n" +
                "ZX473\tZX1577\n" +
                "ZX474\tZX1578\n" +
                "ZX475\tZX1579\n" +
                "ZX476\tZX1580\n" +
                "ZX477\tZX1581\n" +
                "ZX478\tZX1582\n" +
                "ZX479\tZX1583\n" +
                "ZX480\tZX1584\n" +
                "ZX481\tZX1585\n" +
                "ZX482\tZX1586\n" +
                "ZX483\tZX1587\n" +
                "ZX484\tZX1588\n" +
                "ZX485\tZX1589\n" +
                "ZX486\tZX1590\n" +
                "ZX487\tZX1591\n" +
                "ZX488\tZX1592\n" +
                "ZX489\tZX1593\n" +
                "ZX490\tZX1594\n" +
                "ZX491\tZX1595\n" +
                "ZX492\tZX1596\n" +
                "ZX493\tZX1597\n" +
                "ZX494\tZX1598\n" +
                "ZX495\tZX1599\n" +
                "ZX496\tZX1600\n" +
                "ZX497\tZX1601\n" +
                "ZX498\tZX1602\n" +
                "ZX499\tZX1603\n" +
                "ZX500\tZX1604\n" +
                "ZX501\tZX1605\n" +
                "ZX502\tZX1606\n" +
                "ZX503\tZX1607\n" +
                "ZX504\tZX1608\n" +
                "ZX505\tZX1609\n" +
                "ZX506\tZX1610\n" +
                "ZX507\tZX1611\n" +
                "ZX508\tZX1612\n" +
                "ZX509\tZX1613\n" +
                "ZX510\tZX1614\n" +
                "ZX511\tZX1615\n" +
                "ZX512\tZX1616\n" +
                "ZX513\tZX1617\n" +
                "ZX514\tZX1618\n" +
                "ZX515\tZX1619\n" +
                "ZX516\tZX1620\n" +
                "ZX517\tZX1621\n" +
                "ZX518\tZX1622\n" +
                "ZX519\tZX1623\n" +
                "ZX520\tZX1624\n" +
                "ZX521\tZX1625\n" +
                "ZX522\tZX1626\n" +
                "ZX523\tZX1627\n" +
                "ZX524\tZX1628\n" +
                "ZX525\tZX1629\n" +
                "ZX526\tZX1630\n" +
                "ZX527\tZX1631\n" +
                "ZX528\tZX1632\n" +
                "ZX529\tZX1633\n" +
                "ZX530\tZX1634\n" +
                "ZX531\tZX1635\n" +
                "ZX532\tZX1636\n" +
                "ZX533\tZX1637\n" +
                "ZX534\tZX1638\n" +
                "ZX535\tZX1639\n" +
                "ZX536\tZX1640\n" +
                "ZX537\tZX1641\n" +
                "ZX538\tZX1642\n" +
                "ZX539\tZX1643\n" +
                "ZX540\tZX1644\n" +
                "ZX541\tZX1645\n" +
                "ZX542\tZX1646\n" +
                "ZX543\tZX1647\n" +
                "ZX544\tZX1648\n" +
                "ZX545\tZX1649\n" +
                "ZX546\tZX1650\n" +
                "ZX547\tZX1651\n" +
                "ZX548\tZX1652\n" +
                "ZX549\tZX1653\n" +
                "ZX550\tZX1654\n" +
                "ZX551\tZX1655\n" +
                "ZX552\tZX1656\n" +
                "ZX553\tZX1657\n" +
                "ZX554\tZX1658\n" +
                "ZX555\tZX1659\n" +
                "ZX556\tZX1660\n" +
                "ZX557\tZX1661\n" +
                "ZX558\tZX1662\n" +
                "ZX559\tZX1663\n" +
                "ZX560\tZX1664\n" +
                "ZX561\tZX1665\n" +
                "ZX562\tZX1666\n" +
                "ZX563\tZX1667\n" +
                "ZX564\tZX1668\n" +
                "ZX565\tZX1669\n" +
                "ZX566\tZX1670\n" +
                "ZX567\tZX1671\n" +
                "ZX568\tZX1672\n" +
                "ZX569\tZX1673\n" +
                "ZX570\tZX1674\n" +
                "ZX571\tZX1675\n" +
                "ZX572\tZX1676\n" +
                "ZX573\tZX1677\n" +
                "ZX574\tZX1678\n" +
                "ZX575\tZX1679\n" +
                "ZX576\tZX1680\n" +
                "ZX577\tZX1681\n" +
                "ZX578\tZX1682\n" +
                "ZX579\tZX1683\n" +
                "ZX580\tZX1684\n" +
                "ZX581\tZX1685\n" +
                "ZX582\tZX1686\n" +
                "ZX583\tZX1687\n" +
                "ZX584\tZX1688\n" +
                "ZX585\tZX1689\n" +
                "ZX586\tZX1690\n" +
                "ZX587\tZX1691\n" +
                "ZX588\tZX1692\n" +
                "ZX589\tZX1693\n" +
                "ZX590\tZX1694\n" +
                "ZX591\tZX1695\n" +
                "ZX592\tZX1696\n" +
                "ZX593\tZX1697\n" +
                "ZX594\tZX1698\n" +
                "ZX595\tZX1699\n" +
                "ZX596\tZX1700\n" +
                "ZX597\tZX1701\n" +
                "ZX598\tZX1702\n" +
                "ZX599\tZX1703\n" +
                "ZX600\tZX1704\n" +
                "ZX601\tZX1705\n" +
                "ZX602\tZX1706\n" +
                "ZX603\tZX1707\n" +
                "ZX604\tZX1708\n" +
                "ZX605\tZX1709\n" +
                "ZX606\tZX1710\n" +
                "ZX607\tZX1711\n" +
                "ZX608\tZX1712\n" +
                "ZX609\tZX1713\n" +
                "ZX610\tZX1714\n" +
                "ZX611\tZX1715\n" +
                "ZX612\tZX1716\n" +
                "ZX613\tZX1717\n" +
                "ZX614\tZX1718\n" +
                "ZX615\tZX1719\n" +
                "ZX616\tZX1720\n" +
                "ZX617\tZX1721\n" +
                "ZX618\tZX1722\n" +
                "ZX619\tZX1723\n" +
                "ZX620\tZX1724\n" +
                "ZX621\tZX1725\n" +
                "ZX622\tZX1726\n" +
                "ZX623\tZX1727\n" +
                "ZX624\tZX1728\n" +
                "ZX625\tZX1729\n" +
                "ZX626\tZX1730\n" +
                "ZX627\tZX1731\n" +
                "ZX628\tZX1732\n" +
                "ZX629\tZX1733\n" +
                "ZX630\tZX1734\n" +
                "ZX631\tZX1735\n" +
                "ZX632\tZX1736\n" +
                "ZX633\tZX1737\n" +
                "ZX634\tZX1738\n" +
                "ZX635\tZX1739\n" +
                "ZX636\tZX1740\n" +
                "ZX637\tZX1741\n" +
                "ZX638\tZX1742\n" +
                "ZX639\tZX1743\n" +
                "ZX640\tZX1744\n" +
                "ZX641\tZX1745\n" +
                "ZX642\tZX1746\n" +
                "ZX643\tZX1747\n" +
                "ZX644\tZX1748\n" +
                "ZX645\tZX1749\n" +
                "ZX646\tZX1750\n" +
                "ZX647\tZX1751\n" +
                "ZX648\tZX1752\n" +
                "ZX649\tZX1753\n" +
                "ZX650\tZX1754\n" +
                "ZX651\tZX1755\n" +
                "ZX652\tZX1756\n" +
                "ZX653\tZX1757\n" +
                "ZX654\tZX1758\n" +
                "ZX655\tZX1759\n" +
                "ZX656\tZX1760\n" +
                "ZX657\tZX1761\n" +
                "ZX658\tZX1762\n" +
                "ZX659\tZX1763\n" +
                "ZX660\tZX1764\n" +
                "ZX661\tZX1765\n" +
                "ZX662\tZX1766\n" +
                "ZX663\tZX1767\n" +
                "ZX664\tZX1768\n" +
                "ZX665\tZX1769\n" +
                "ZX666\tZX1770\n" +
                "ZX667\tZX1771\n" +
                "ZX668\tZX1772\n" +
                "ZX669\tZX1773\n" +
                "ZX670\tZX1774\n" +
                "ZX671\tZX1775\n" +
                "ZX672\tZX1776\n" +
                "ZX673\tZX1777\n" +
                "ZX674\tZX1778\n" +
                "ZX675\tZX1779\n" +
                "ZX676\tZX1780\n" +
                "ZX677\tZX1781\n" +
                "ZX678\tZX1782\n" +
                "ZX679\tZX1783\n" +
                "ZX680\tZX1784\n" +
                "ZX681\tZX1785\n" +
                "ZX682\tZX1786\n" +
                "ZX683\tZX1787\n" +
                "ZX684\tZX1788\n" +
                "ZX685\tZX1789\n" +
                "ZX686\tZX1790\n" +
                "ZX687\tZX1791\n" +
                "ZX688\tZX1792\n" +
                "ZX689\tZX1793\n" +
                "ZX690\tZX1794\n" +
                "ZX691\tZX1795\n" +
                "ZX692\tZX1796\n" +
                "ZX693\tZX1797\n" +
                "ZX694\tZX1798\n" +
                "ZX695\tZX1799\n" +
                "ZX696\tZX1800\n" +
                "ZX697\tZX1801\n" +
                "ZX698\tZX1802\n" +
                "ZX699\tZX1803\n" +
                "ZX700\tZX1804\n" +
                "ZX701\tZX1805\n" +
                "ZX702\tZX1806\n" +
                "ZX703\tZX1807\n" +
                "ZX704\tZX1808\n" +
                "ZX705\tZX1809\n" +
                "ZX706\tZX1810\n" +
                "ZX707\tZX1811\n" +
                "ZX708\tZX1812\n" +
                "ZX709\tZX1813\n" +
                "ZX710\tZX1814\n" +
                "ZX711\tZX1815\n" +
                "ZX712\tZX1816\n" +
                "ZX713\tZX1817\n" +
                "ZX714\tZX1818\n" +
                "ZX715\tZX1819\n" +
                "ZX716\tZX1820\n" +
                "ZX717\tZX1821\n" +
                "ZX718\tZX1822\n" +
                "ZX719\tZX1823\n" +
                "ZX720\tZX1824\n" +
                "ZX721\tZX1825\n" +
                "ZX722\tZX1826\n" +
                "ZX723\tZX1827\n" +
                "ZX724\tZX1828\n" +
                "ZX725\tZX1829\n" +
                "ZX726\tZX1830\n" +
                "ZX727\tZX1831\n" +
                "ZX728\tZX1832\n" +
                "ZX729\tZX1833\n" +
                "ZX730\tZX1834\n" +
                "ZX731\tZX1835\n" +
                "ZX732\tZX1836\n" +
                "ZX733\tZX1837\n" +
                "ZX734\tZX1838\n" +
                "ZX735\tZX1839\n" +
                "ZX736\tZX1840\n" +
                "ZX737\tZX1841\n" +
                "ZX738\tZX1842\n" +
                "ZX739\tZX1843\n" +
                "ZX740\tZX1844\n" +
                "ZX741\tZX1845\n" +
                "ZX742\tZX1846\n" +
                "ZX743\tZX1847\n" +
                "ZX744\tZX1848\n" +
                "ZX745\tZX1849\n" +
                "ZX746\tZX1850\n" +
                "ZX747\tZX1851\n" +
                "ZX748\tZX1852\n" +
                "ZX749\tZX1853\n" +
                "ZX750\tZX1854\n" +
                "ZX751\tZX1855\n" +
                "ZX752\tZX1856\n" +
                "ZX753\tZX1857\n" +
                "ZX754\tZX1858\n" +
                "ZX755\tZX1859\n" +
                "ZX756\tZX1860\n" +
                "ZX757\tZX1861\n" +
                "ZX758\tZX1862\n" +
                "ZX759\tZX1863\n" +
                "ZX760\tZX1864\n" +
                "ZX761\tZX1865\n" +
                "ZX762\tZX1866\n" +
                "ZX763\tZX1867\n" +
                "ZX764\tZX1868\n" +
                "ZX765\tZX1869\n" +
                "ZX766\tZX1870\n" +
                "ZX767\tZX1871\n" +
                "ZX768\tZX1872\n" +
                "ZX769\tZX1873\n" +
                "ZX770\tZX1874\n" +
                "ZX771\tZX1875\n" +
                "ZX772\tZX1876\n" +
                "ZX773\tZX1877\n" +
                "ZX774\tZX1878\n" +
                "ZX775\tZX1879\n" +
                "ZX776\tZX1880\n" +
                "ZX777\tZX1881\n" +
                "ZX778\tZX1882\n" +
                "ZX779\tZX1883\n" +
                "ZX780\tZX1884\n" +
                "ZX781\tZX1885\n" +
                "ZX782\tZX1886";
        String[] temps = table.split("\n");
        String[] temp = null;
        for (int i = 1; i < temps.length; i++) {
            temp = temps[i].split("\t");
            trans.put(temp[0], temp[1]);
        }
        try {
            BufferedReader br = IOUtils.getTextReader(infile1);
            for (int i = 1; i <= 96; i++) {
                String count = br.readLine();
                System.out.println(count + "\t");
                System.out.println(trans.get(count) + "\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static void main(String args[]) throws IOException {
        new test();
    }
}

