package daxing.selection.model;

import daxing.common.ChrPos;
import daxing.common.IOTool;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.Predicate;

/**
 * multiple alignment records without indel (blocks) in a maf file
 */
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
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
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

    public static Comparator<Block> getComparator(String subgenome){
        Comparator<Block> chrComparator=Comparator.comparing(block -> block.getChr(subgenome));
        Comparator<Block> chrPlusComparator=chrComparator.thenComparing(block -> block.getIfMinus(subgenome));
        Comparator<Block> chrPlusPosComparator=chrPlusComparator.thenComparing(block -> block.getSeqStart(subgenome));
        return chrPlusPosComparator;
    }

    /**
     * sort by subgenome A
     * @param chrA
     * @param posA
     * @return
     */
    public int getBlockIndex(String chrA, int posA){
        this.sortBy("A");
        Block query=new Block(chrA, posA);
        int hit=Collections.binarySearch(this.blockList, query, getComparator("A"));
        int index=hit;
        if (hit < -1){
            index=-hit-2;
            index= isInThisBlock(index, chrA, posA) ? index : hit;
        }
        return index;
    }

    /**
     *
     * @param chrA
     * @param posA one-based position in subgenome A
     * @return posB one-based position in subgenome B
     */
    public ChrPos getChrPosB(String chrA, int posA){
        int index=getBlockIndex(chrA, posA);
        if (index < 0) return null;
        if (!this.blockList.get(index).haveAlignment("B")) return null;
        return this.blockList.get(index).getChrPosB(posA);
    }

    /**
     *
     * @param chrA
     * @param posA one-based position in subgenome A
     * @return posD one-based position in subgenome D
     */
    public ChrPos getChrPosD(String chrA, int posA){
        int index=getBlockIndex(chrA, posA);
        if (index < 0) return null;
        if (!this.blockList.get(index).haveAlignment("D")) return null;
        return this.blockList.get(index).getChrPosD(posA);
    }

    private boolean isInThisBlock(int blockIndex, String chrA, int posA){
        Block block=this.blockList.get(blockIndex);
        String chr=block.getChr("A");
        if (!chr.equals(chrA)) return false;
        int seqStart=block.getSeqStart("A");
        int seqLen=block.getSeqLen("A");
        int zeroBasedPosA=posA - 1;
        if (zeroBasedPosA < seqStart) return false;
        if (zeroBasedPosA >= seqStart + seqLen) return false;
        return true;
    }

    /**
     * first removed, then sort
     * @param subgenome
     */
    public void sortBy(String subgenome){
        this.removeLackingSubgenome(subgenome);
        Collections.sort(this.blockList, BlockOperation.getComparator(subgenome));
    }

    public static void initialize(String triadsMafFile, String outFile){
        List<Block> blockList;
        String line=null;
        try (BufferedReader br = IOTool.getReader(triadsMafFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
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
