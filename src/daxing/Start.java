package daxing;

import daxing.applets.ScriptMethods;
import daxing.common.*;
import daxing.load.GeneAnnotation;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import pgl.graphcis.tablesaw.TablesawUtils;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Table;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) {
//        System.out.println("0: "+       0b1);
//        System.out.println("1: "+      0b10);
//        System.out.println("2: "+     0b100);
//        System.out.println("3: "+    0b1000);
//        System.out.println("4: "+   0b10000);
//        System.out.println("5: "+  0b100000);
//        System.out.println("6: "+ 0b1000000);
//        System.out.println("7: "+0b10000000);
//        short chr=1;
//        int pos=10;
//        char ref='A';
//        char alt='T';
//        String info="test";
//        BiSNP snp1=new BiSNP(chr, pos, ref, alt, info);
//        byte refAlleleByte=snp1.getReferenceAlleleByte();
//        byte altAlleleByte=snp1.getAlternativeAlleleByte();
//        char refAllele=snp1.getReferenceAlleleBase();
//        char altAllele=snp1.getAlternativeAlleleBase();
//        byte refAlleleFuture=snp1.getReferenceAlleleFeature();
//        byte altAlleleFuture=snp1.getAlternativeAlleleFeature();
    }
}