package xiaohan.rareallele;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;

public class VCFsplit {

    public VCFsplit(){
        this.split();
    }

    public void split(){
        String positionFileS =
                "1	0	471304005	chr1A	0	471304005\n" +
                        "2	0	122798051	chr1A	471304005	594102056\n" +
                        "3	0	438720154	chr1B	0	438720154\n" +
                        "4	0	251131716	chr1B	438720154	689851870\n" +
                        "5	0	452179604	chr1D	0	452179604\n" +
                        "6	0	43273582	chr1D	452179604	495453186\n" +
                        "7	0	462376173	chr2A	0	462376173\n" +
                        "8	0	318422384	chr2A	462376173	780798557\n" +
                        "9	0	453218924	chr2B	0	453218924\n" +
                        "10	0	348037791	chr2B	453218924	801256715\n" +
                        "11	0	462216879	chr2D	0	462216879\n" +
                        "12	0	189635730	chr2D	462216879	651852609\n" +
                        "13	0	454103970	chr3A	0	454103970\n" +
                        "14	0	296739669	chr3A	454103970	750843639\n" +
                        "15	0	448155269	chr3B	0	448155269\n" +
                        "16	0	382674495	chr3B	448155269	830829764\n" +
                        "17	0	476235359	chr3D	0	476235359\n" +
                        "18	0	139317064	chr3D	476235359	615552423\n" +
                        "19	0	452555092	chr4A	0	452555092\n" +
                        "20	0	292033065	chr4A	452555092	744588157\n" +
                        "21	0	451014251	chr4B	0	451014251\n" +
                        "22	0	222603248	chr4B	451014251	673617499\n" +
                        "23	0	451004620	chr4D	0	451004620\n" +
                        "24	0	58852447	chr4D	451004620	509857067\n" +
                        "25	0	453230519	chr5A	0	453230519\n" +
                        "26	0	256543224	chr5A	453230519	709773743\n" +
                        "27	0	451372872	chr5B	0	451372872\n" +
                        "28	0	261776885	chr5B	451372872	713149757\n" +
                        "29	0	451901030	chr5D	0	451901030\n" +
                        "30	0	114179647	chr5D	451901030	566080677\n" +
                        "31	0	452440856	chr6A	0	452440856\n" +
                        "32	0	165638404	chr6A	452440856	618079260\n" +
                        "33	0	452077197	chr6B	0	452077197\n" +
                        "34	0	268911281	chr6B	452077197	720988478\n" +
                        "35	0	450509124	chr6D	0	450509124\n" +
                        "36	0	23083594	chr6D	450509124	473592718\n" +
                        "37	0	450046986	chr7A	0	450046986\n" +
                        "38	0	286659250	chr7A	450046986	736706236\n" +
                        "39	0	453822637	chr7B	0	453822637\n" +
                        "40	0	296797748	chr7B	453822637	750620385\n" +
                        "41	0	453812268	chr7D	0	453812268\n" +
                        "42	0	184873787	chr7D	453812268	638686055";

        HashMap<Integer,Integer> ChromosomeMap = new HashMap<Integer, Integer>();
        HashMap<Integer,Integer> ChromosomeLength = new HashMap<Integer, Integer>();
        String[] temps = positionFileS.split("\n");
        String[] temp = null;
        for (int i = 0; i < temps.length; i ++){
            temp = temps[i].split("\t");
            ChromosomeMap.put(Integer.parseInt(temp[0]),Integer.parseInt(temp[5]));
            ChromosomeLength.put(Integer.parseInt(temp[0]),Integer.parseInt(temp[4]));
        }
        try {
            String infile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/AA.txt";
            BufferedReader br = xiaohan.rareallele.IOUtils.getTextReader(infile);
            String temp1 = null;
            while ((temp1 = br.readLine())!= null){
                if(temp1.startsWith("#")) {
                    continue;
                }
                String[] temp2 = temp1.split("\t");
                int chr = Integer.parseInt(temp2[0]);
                int chrNumber = chr*2 - 1;
                int position = Integer.parseInt(temp2[1]);
                if (ChromosomeMap.get(chrNumber) <= position)
                {
                    continue;
                }
                if(ChromosomeMap.get(chrNumber)>position){
                    position = position-ChromosomeLength.get(chrNumber);
                    chr = chrNumber;
                    String chr1 = String.valueOf(chr);
                    String position1 = String.valueOf(position);
                    System.out.println(chr1+"\t"+position1);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main (String[] args){
        new VCFsplit();
    }
}
