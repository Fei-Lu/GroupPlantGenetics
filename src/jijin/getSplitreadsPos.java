/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jijin;

import htsjdk.samtools.*;
import pgl.infra.range.Range;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author feilu
 */

public class getSplitreadsPos {
    public getSplitreadsPos() {

    }

    public static void test(String file, String outfile) throws IOException {
        List<SplitReads> info =new LinkedList<>();
        File bamFile = new File(file);
        FileWriter fw = new FileWriter(outfile);
        BufferedWriter bw = new BufferedWriter(fw);
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
        //System.out.println(sr);
        //BAMFileReader br = new BAMFileReader(bamFile, null, false, false, ValidationStringency.SILENT,new DefaultSAMRecordFactory());
        SAMRecordIterator reads = sr.iterator();
        LinkedList range = new LinkedList();
        //int i=0;
        //StringBuilder sb = new StringBuilder(); //字符串缓冲生成器利用
        while (reads.hasNext()) {
            SAMRecord eachreads = reads.next();
            //int flag = eachreads.getFlags();
            //String flags = Integer.toBinaryString(flag);
            //int x = getBinaryValueAtNthBit(Integer.valueOf(flags), 11);
            //System.out.println(eachreads.getCigar());
            //System.out.println(eachreads.getReferenceIndex());
//            System.out.println(eachreads.getAttribute("SA"));
//            bw.write(String.valueOf(eachreads.getAttribute("SA")));
//            bw.newLine();
            if (eachreads.getAttribute("SA")!=null){
                //System.out.println(String.valueOf(eachreads.getReferenceName()));
                //System.out.println(String.valueOf(eachreads.getAttribute("SA").toString().split(",")[0]));
                bw.write(String.valueOf(eachreads.getReferenceName()));
                bw.write("\t");
                bw.write(String.valueOf(eachreads.getAttribute("SA").toString().split(",")[0]));
                bw.write("\t");
                bw.write(String.valueOf(eachreads.getStart()));
                bw.write("\t");
                bw.write(String.valueOf(eachreads.getEnd()));
                bw.newLine();
                if(String.valueOf(eachreads.getReferenceName()).equals(String.valueOf(eachreads.getAttribute("SA").toString().split(",")[0]))){
                    //System.out.println(String.valueOf(eachreads.getReferenceName()));
                    //info.add(new SplitReads(String.valueOf(eachreads.getAttribute("SA")),String.valueOf(eachreads.getReadName()),eachreads.getReferenceIndex(),eachreads.getStart(),eachreads.getEnd())) ;
                    System.out.println(eachreads.getReadName());

                }
                else{
                    continue;
                }


                //System.out.println(eachreads.getAttribute("SA").toString().split(",")[0]);
               //System.out.println(String.valueOf(eachreads.getAttribute("SA")));
//                System.out.println(String.valueOf(eachreads.getReadName()));
//                System.out.println(String.valueOf(eachreads.getStart()));
//                System.out.println(String.valueOf(eachreads.getEnd()));
                //info.add(String.valueOf(eachreads.getAttribute("SA"))+"-"+String.valueOf(eachreads.getReadName()));
                //info.add(i);

                //info.add(i);//将对象作为一个元素放到链表中
                //sb.append(String.valueOf(eachreads.getAttribute("SA"))).append("aaaaaaaaaaa");
                //System.out.println(sb.toString());

            }
            //i++;



            //info.get(1).getchrom();
//            if (getBinaryValueAtNthBit(flag, 1) == 1 &&
//                    getBinaryValueAtNthBit(flag, 2) == 0 &&
//                    getBinaryValueAtNthBit(flag, 3) == 0 &&
//                    getBinaryValueAtNthBit(flag, 4) == 0 &&
//                    getBinaryValueAtNthBit(flag, 10) == 0 &&
//                    getBinaryValueAtNthBit(flag, 11) == 0 &&
//                    getBinaryValueAtNthBit(flag, 12) == 1) {
//





//            }
        }

        System.out.println("一条也没有啊");

//        Iterator<SplitReads> i =info.iterator();
//        while (i.hasNext()){
//            System.out.println(i.next().gettags());
//        }
        bw.flush();
        bw.close();
        //return info;

    }

    public static void readsRange(int chr, int start, int end){
        List<Range> cdsList = new ArrayList();
        cdsList.add(new Range(chr, start, end));
    }

    public static void getInterval(int s1,int e1,int s2,int e2){
        if(e1 < s2 || e2 < s1) {
            System.out.println("两集合无交集");
        }else{
            System.out.println("两集合的交集为：[" + Math.max(s1,s2) + ","  +  Math.min(e1,e2) + "]");
        }
    }
    public static int getBinaryValueAtNthBit(int flag, int bitPosRightMost) {

        int shift = (bitPosRightMost - 1);

        return (flag >> shift) & 1;

    }





    public static void main(String[] args) throws IOException {

        long a=System.currentTimeMillis();
        getSplitreadsPos gs =new getSplitreadsPos();
        //gs.test(args[0],args[1]);
        gs.test("E:\\Desktop\\split_reads\\test2.bam","E:\\Desktop\\split_reads\\splitread2.txt");

        System.out.println("\r<br> 执行耗时 : "+(System.currentTimeMillis()-a)/1000f+" 秒 ");
//        int binaryValueAtNthBit = gs.getBinaryValueAtNthBit(97,2 );//调取二进制从右数第 参2 个位置的数
//        System.out.println(binaryValueAtNthBit);
//        System.out.println(Integer.toBinaryString(97));

//        ArrayList list = new ArrayList();
//        // 增加：add() 将指定对象存储到容器中
//        list.add("计算机网络");
//        list.add("现代操作系统");
//        list.add("java编程思想");
//        list.add("java核心技术");
//        list.add("java语言程序设计");
//        System.out.println(list);
//        Iterator it = list.iterator();
//        while (it.hasNext()) {
//            String next = (String) it.next();
//            System.out.println(next);} //迭代器的使用举例


    }
}
