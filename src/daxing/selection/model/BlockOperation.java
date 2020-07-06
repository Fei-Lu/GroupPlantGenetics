package daxing.selection.model;

import daxing.common.IOTool;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class BlockOperation {

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
