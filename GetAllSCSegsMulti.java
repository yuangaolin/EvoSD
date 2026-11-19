import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class GetAllSCSegsMulti {

    static class Segment {
        int segIndex, start, end;
        String seqMasked, seqRemoved;
        int numN;
        int length;

        Segment(int segIndex, int start, int end,
                String seqMasked, String seqRemoved,
                int numN, int length) {
            this.segIndex = segIndex;
            this.start = start;
            this.end = end;
            this.seqMasked = seqMasked;
            this.seqRemoved = seqRemoved;
            this.numN = numN;
            this.length = length;
        }
    }

    static class ChrResult {
        int order;
        List<Segment> segments = new ArrayList<>();
        double scnLen = 0;
        int numN = 0;

        ChrResult(int order) {
            this.order = order;
        }
    }

    public static void main(String[] args) {
        try {
            long startTime = System.currentTimeMillis();

            int kmerLength = Integer.parseInt(args[5]);
            int extendSpan = 500;

            String chrFile = args[0];
            String scnFile = args[1];
            String genomeFile = args[2];
            String outMasked = args[3];
            String outRemoved = args[4];
            int threads = args.length >= 7
                    ? Integer.parseInt(args[6])
                    : Runtime.getRuntime().availableProcessors();

            Scanner inChr = new Scanner(new File(chrFile));
            Scanner inSCN = new Scanner(new File(scnFile));
            Scanner inGenome = new Scanner(new File(genomeFile));

            BufferedWriter writer = new BufferedWriter(new FileWriter(outMasked));
            BufferedWriter writerRe = new BufferedWriter(new FileWriter(outRemoved));

            ExecutorService pool = Executors.newFixedThreadPool(threads);
            List<Future<ChrResult>> futures = new ArrayList<>();

            int segIndex = -1;
            int order = 0;

            while (inChr.hasNextLine()) {
                inChr.nextLine();
                String chr = inChr.nextLine();

                inSCN.nextLine();
                String strSCN = inSCN.nextLine();
                segIndex++;

                inGenome.nextLine();
                String genome = inGenome.nextLine();

                final int curSegIndex = segIndex;
                final int curOrder = order++;

                futures.add(pool.submit(() ->
                        processOneChromosome(chr, strSCN, genome,
                                kmerLength, extendSpan,
                                curSegIndex, curOrder)
                ));
            }

            pool.shutdown();

            futures.sort(Comparator.comparingInt(f -> {
                try {
                    return f.get().order;
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }));

            double SCNLength = 0;
            int numN = 0;
            int indexW = 0;

            for (Future<ChrResult> f : futures) {
                ChrResult r = f.get();
                for (Segment s : r.segments) {
                    writer.write(">" + indexW + "  " +
                            s.segIndex + "  " +
                            s.start + "  " +
                            s.end);
                    writer.newLine();

                    writer.write(s.seqMasked);
                    writer.newLine();

                    writerRe.write(">" + indexW + "  " +
                            s.segIndex + "  " +
                            s.start + "  " +
                            s.end);
                    writerRe.newLine();

                    writerRe.write(s.seqRemoved);
                    writerRe.newLine();

                    indexW++;
                }
                SCNLength += r.scnLen;
                numN += r.numN;
            }

            writer.close();
            writerRe.close();
            inChr.close();
            inSCN.close();
            inGenome.close();

            long endTime = System.currentTimeMillis();
            System.out.println("SDquest identified putative SDs!");
            System.out.println("CR/SCNlength= " + numN / SCNLength +
                    "   runTime: " + (endTime - startTime));

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    // ================= 原算法逐字等价封装 =================
    private static ChrResult processOneChromosome(String chr, String strSCN, String genome,
                                                  int kmerLength, int extendSpan,
                                                  int segIndex, int order) {

        ChrResult result = new ChrResult(order);

        int start, end = 0, index;
        double SCNLength = 0;
        int numN = 0;
        char cur;
        int extendL = 0;
        boolean findN;

        start = 0;
        end = 0;
        index = 0;
        int length = strSCN.length();

        while (index < length) {
            if (strSCN.charAt(index) == '1') {
                start = index;
                index++;
                while ((index < length) && (strSCN.charAt(index) == '1')) {
                    index++;
                }
                end = index + (kmerLength - 1);

                extendL = 0;
                findN = false;
                while (extendL < extendSpan && start > 0 && !findN) {
                    start--;
                    if (genome.charAt(start) == 'N' || genome.charAt(start) == 'n') {
                        findN = true;
                        start++;
                        while (chr.charAt(start) == 'N') {
                            start++;
                        }
                    } else if (chr.charAt(start) != 'N') {
                        extendL++;
                    }
                }
                if (start == 0 && chr.charAt(start) == 'N') {
                    start++;
                    while (chr.charAt(start) == 'N') {
                        start++;
                    }
                }

                extendL = 0;
                findN = false;
                while (extendL < extendSpan && end < chr.length() && !findN) {
                    if (genome.charAt(end) == 'N' || genome.charAt(end) == 'n') {
                        findN = true;
                        while (chr.charAt(end) == 'N') {
                            end--;
                        }
                    } else if (chr.charAt(end) != 'N') {
                        extendL++;
                    }
                    end++;
                }
                if (end == chr.length() && chr.charAt(end - 1) == 'N') {
                    end--;
                    while (chr.charAt(end - 1) == 'N') {
                        end--;
                    }
                }

                StringBuilder maskedSeq = new StringBuilder();
                StringBuilder removedSeq = new StringBuilder();

                SCNLength += (end - start);

                for (int i = start; i < end; i++) {
                    cur = chr.charAt(i);
                    maskedSeq.append(cur);
                    if (cur != 'N') {
                        removedSeq.append(cur);
                    } else {
                        numN++;
                    }
                }

                result.segments.add(new Segment(
                        segIndex, start, end,
                        maskedSeq.toString(),
                        removedSeq.toString(),
                        numN, end - start
                ));
            }
            index++;
        }

        result.scnLen = SCNLength;
        result.numN = numN;
        return result;
    }
}
