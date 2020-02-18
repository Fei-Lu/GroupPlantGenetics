package feilu;

import java.util.BitSet;
import pgl.infra.utils.Benchmark;

public class BitSetTest {

    public BitSetTest () {
        //this.test1();
        this.test2();
    }
    
    public void test2 () {
        //Set is faster in OpenBitSet when the size is small;
        //And, cardinality is generally faster in BitSet
        int size = 65536;
        OpenBitSet ob = new OpenBitSet(size);
        BitSet bs = new BitSet(size);
        System.out.println("Set:");
        long start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            bs.set(i);
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            ob.set(i);
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            ob.fastSet(i);
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        System.out.println("Get:");
        start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            bs.get(i);
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            ob.get(i);
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            ob.fastGet(i);
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        System.out.println("Flip:");
        start = System.nanoTime();
        bs.flip(0, size);
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        ob.flip(0,size);
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        System.out.println("And:");
        start = System.nanoTime();
        bs.and(bs);
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        ob.and(ob);
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        System.out.println("Cardinality:");
        start = System.nanoTime();
        bs.cardinality();
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        ob.cardinality();
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();

    }
    
    public void test1 () {
        BitSet bs = new BitSet(10);
        bs.set(72);
        System.out.println(bs.length());
        System.out.println(bs.size());
    }

}
