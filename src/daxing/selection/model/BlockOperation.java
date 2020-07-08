package daxing.selection.model;

import daxing.common.IOTool;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;

public class BlockOperation {

    List<Block> blockList;

    public BlockOperation(String triadsMafFile){
        this.blockList=new ArrayList<>();
        List<Block> blockList;
        String line=null;
        try (BufferedReader br = IOTool.getReader(triadsMafFile)) {
            br.readLine();
            Block block;
            List<String> lines=new ArrayList<>();
            while ((line=br.readLine())!=null){
                if (line.startsWith("s")){
                    lines.add(line);
                }else if (line.length()==0){
                    block=new Block(lines);
                    blockList=block.subsetBlockWithoutIndel();
                    this.blockList.addAll(blockList);
                    lines=new ArrayList<>();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public List<Block> getBlockList() {
        return blockList;
    }

    public void write(String outFile){
        try (BufferedWriter bw = IOTool.getTextWriter(outFile)) {
            for (int i = 0; i < blockList.size(); i++) {
                bw.write(blockList.get(i).toString());
                bw.newLine();
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * If any block lacks corresponding subgenome alignment, then remove it
     * @param subgenome A B D
     */
    public void removeLackingSubgenome(String subgenome){
        Predicate<Block> p = block -> block.haveAlignment(subgenome);
        this.blockList.removeIf(p.negate());
    }

    /**
     * first removed, then sort
     * @param subgenome
     */
    public void sortBy(String subgenome){
        this.removeLackingSubgenome(subgenome);
        Comparator<Block> chrComparator=Comparator.comparing(block -> block.getChr(subgenome));
        Comparator<Block> chrPlusComparator=chrComparator.thenComparing(block -> block.getIfMinus(subgenome));
        Comparator<Block> chrPlusPosComparator=chrPlusComparator.thenComparing(block -> block.getSeqStart(subgenome));
        Collections.sort(this.blockList, chrPlusPosComparator);
    }

    public static void initialize(String triadsMafFile, String outFile){
        List<Block> blockList;
        String line=null;
        try (BufferedReader br = IOTool.getReader(triadsMafFile);
             BufferedWriter bw = IOTool.getTextWriter(outFile)) {
            br.readLine();
            Block block;
            List<String> lines=new ArrayList<>();
            while ((line=br.readLine())!=null){
                if (line.startsWith("s")){
                    lines.add(line);
                }else if (line.length()==0){
                    block=new Block(lines);
                    blockList=block.subsetBlockWithoutIndel();
                    for (int i = 0; i < blockList.size(); i++) {
                        bw.write(blockList.get(i).toString());
                        bw.newLine();
                        bw.newLine();
                    }
                    lines=new ArrayList<>();
                }
            }
            bw.flush();
        } catch (IOException e) {
            System.out.println(line);
            e.printStackTrace();
        }
    }

}
