package feilu;

import pgl.infra.utils.Benchmark;

import java.util.ArrayList;
import java.util.List;

public class ListAndArray {

    public ListAndArray() {
        this.test();
    }

    public void test () {
        int size = 100000000;
        List<Integer> l = new ArrayList(size);
        int[] a = new int[size];
        long start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            a[i] = i;
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            l.add(i);
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        int c;
        start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            c = a[i];
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
        start = System.nanoTime();
        for (int i = 0; i < size; i++) {
            l.get(i);
        }
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
    }
}
