package daxing.v2.localAncestryInfer.laidp;

import java.util.BitSet;
import java.util.concurrent.CountDownLatch;

public class MAFCalculator {

    private static final int BLOCK_SIZE = 1 << 20;
    private static final int MASK = BLOCK_SIZE - 1;

    public static double[] calculateMAF(BitSet[][] bs, int variants, int numThreads) {
        int[] counts0 = new int[variants];
        int[] counts1 = new int[variants];

        int blocks = (variants + BLOCK_SIZE - 1) / BLOCK_SIZE;
        CountDownLatch latch = new CountDownLatch(numThreads);

        for (int t = 0; t < numThreads; t++) {
            int start = t * blocks / numThreads;
            int end = (t + 1) * blocks / numThreads;
            new Thread(() -> {
                for (int i = start; i < end; i++) {
                    int blockStart = i * BLOCK_SIZE;
                    int blockEnd = Math.min(blockStart + BLOCK_SIZE, variants);
                    for (int j = 0; j < bs.length; j++) {
                        BitSet phasing1 = bs[j][0];
                        BitSet phasing2 = bs[j][1];
                        for (int k = blockStart; k < blockEnd; k++) {
                            int idx = k & MASK;
                            if (idx == 0) {
                                if (phasing1.get(k)) {
                                    counts0[k]++;
                                }
                                if (phasing2.get(k)) {
                                    counts1[k]++;
                                }
                            } else {
                                if (phasing1.get(k)) {
                                    counts0[k]++;
                                }
                                if (phasing2.get(k)) {
                                    counts1[k]++;
                                }
                            }
                        }
                    }
                }
                latch.countDown();
            }).start();
        }

        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        double[] mafs = new double[variants];
        for (int i = 0; i < variants; i++) {
            mafs[i] = (counts0[i] + counts1[i]) / (2.0 * bs.length);
        }

        return mafs;
    }
}

