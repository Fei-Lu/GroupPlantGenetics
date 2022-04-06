package daxing.v2.TREE;

import com.google.common.base.Enums;
import daxing.common.factors.WheatLineage;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import daxing.common.vmap2Group.GroupBySubcontinent;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang3.EnumUtils;
import org.apache.commons.lang3.StringUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;

public class Misc {

    public static void extractPopFromPhylip(String phylipDir, String taxaInfoFile, String outDir){
        List<File> fileList = IOTool.getFileListInDirEndsWith(phylipDir, ".phy");
        Map<String, String> taxaGroupbyContinentMap= RowTableTool.getMap(taxaInfoFile, 0, 36);
        IntStream.range(0,fileList.size()).forEach(e->{
            try (BufferedReader br = IOTool.getReader(fileList.get(e))) {
                int chrID = Integer.parseInt(fileList.get(e).getName().substring(3,6));
                WheatLineage subgenome = WheatLineage.valueOf(RefV1Utils.getSubgenomeFromChrID((short)chrID));
                EnumSet<GroupBySubcontinent> enumSet = GroupBySubcontinent.getSubgenomeGroupByContinent(subgenome);
                EnumMap<GroupBySubcontinent, BufferedWriter> groupBWMap= new EnumMap<>(GroupBySubcontinent.class);
                BufferedWriter bw;
                for (GroupBySubcontinent group : enumSet){
                    bw = IOTool.getWriter(new File(outDir, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_vmap2.1_"+group.name()+".phy"));
                    groupBWMap.put(group, bw);
                }
                String line;
                GroupBySubcontinent group;
                String[] temp;
                line=br.readLine();
                temp = StringUtils.split(line, " ");
                StringBuilder sb = new StringBuilder();
                Map<GroupBySubcontinent,Integer> groupCountMap=getTotalNumTaxa();
                for (Map.Entry<GroupBySubcontinent,BufferedWriter> entry : groupBWMap.entrySet()){
                    sb.setLength(0);
                    group = entry.getKey();
                    bw = entry.getValue();
                    sb.append(groupCountMap.get(group)).append(" ").append(temp[1]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                while ((line=br.readLine())!=null){
                    temp = StringUtils.split(line.substring(0,20)," ");
                    if (! EnumUtils.isValidEnum(GroupBySubcontinent.class, taxaGroupbyContinentMap.get(temp[0]))) continue;
                    group=GroupBySubcontinent.valueOf(taxaGroupbyContinentMap.get(temp[0]));
                    bw = groupBWMap.get(group);
                    bw.write(line);
                    bw.newLine();
                }
                for (Map.Entry<GroupBySubcontinent,BufferedWriter> entry : groupBWMap.entrySet()){
                    bw = entry.getValue();
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });

    }

    private static Map<GroupBySubcontinent, Integer> getTotalNumTaxa(){
        Map<GroupBySubcontinent, Integer> map = new HashMap<>();
        map.put(GroupBySubcontinent.WE, 91);
        map.put(GroupBySubcontinent.DE,52);
        map.put(GroupBySubcontinent.FTT,50);
        map.put(GroupBySubcontinent.AT,36);
        map.put(GroupBySubcontinent.LR_AF,18);
        map.put(GroupBySubcontinent.LR_AM,26);
        map.put(GroupBySubcontinent.LR_CSA,19);
        map.put(GroupBySubcontinent.LR_WA,56);
        map.put(GroupBySubcontinent.LR_EU,58);
        map.put(GroupBySubcontinent.LR_EA,127);
        map.put(GroupBySubcontinent.CL,223);
        return map;
    }

    public static void extractFTT(String v2Dir, String outDir, String taxaInfoFile){
        Map<String, String> taxonMap=RowTableTool.getMap(taxaInfoFile, 0, 36);
        List<File> files =IOTool.getFileListInDirEndsWith(v2Dir, ".vcf.gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz","_FTT.vcf.gz")).toArray(String[]::new);
        IntStream.range(0,files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                List<String> temp, tem,te;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                temp = PStringUtils.fastSplit(line);
                TIntArrayList fttIndexList = new TIntArrayList();
                StringBuilder sb = new StringBuilder();
                sb.append(String.join("\t", temp.subList(0,9))).append("\t");
                for (int i = 9; i < temp.size(); i++) {
                    if (!taxonMap.get(temp.get(i)).equals("FTT")) continue;
                    fttIndexList.add(i);
                    sb.append(temp.get(i)).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
                String genotypeField;
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    sb.setLength(0);
                    sb.append(String.join("\t", temp.subList(0,9))).append("\t");
                    int sumGenotype=0;
                    for (int i = 0; i < fttIndexList.size(); i++) {
                        genotypeField = temp.get(fttIndexList.get(i));
                        sb.append(genotypeField).append("\t");
                        if (genotypeField.startsWith("./.")) continue;
                        tem = PStringUtils.fastSplit(genotypeField,":");
                        te = PStringUtils.fastSplit(tem.get(0), "/");
                        sumGenotype = sumGenotype + Integer.parseInt(te.get(0)) + Integer.parseInt(te.get(1));
                    }
                    if (sumGenotype == 0) continue;
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });

    }
}
