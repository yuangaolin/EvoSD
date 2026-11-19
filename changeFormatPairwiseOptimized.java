import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class changeFormatPairwiseOptimized {

    // ---------------------------
    // Fast field parsing (no regex)
    // ---------------------------
    static int nextNonSpace(String s, int i) {
        int n = s.length();
        while (i < n && s.charAt(i) <= ' ') i++;
        return i;
    }

    static int nextSpace(String s, int i) {
        int n = s.length();
        while (i < n && s.charAt(i) > ' ') i++;
        return i;
    }

    // field k range [start,end)
    static int[] fieldRange(String s, int k) {
        int i = 0;
        for (int f = 0; f <= k; f++) {
            i = nextNonSpace(s, i);
            int j = nextSpace(s, i);
            if (f == k) return new int[]{i, j};
            i = j;
        }
        return new int[]{s.length(), s.length()};
    }

    static int parseIntField(String s, int k) {
        int[] r = fieldRange(s, k);
        int i = r[0], end = r[1];
        int sign = 1;
        if (i < end && s.charAt(i) == '-') {
            sign = -1;
            i++;
        }
        int val = 0;
        for (; i < end; i++) {
            val = val * 10 + (s.charAt(i) - '0');
        }
        return val * sign;
    }

    static double parseDoubleField(String s, int k) {
        int[] r = fieldRange(s, k);
        return Double.parseDouble(s.substring(r[0], r[1]));
    }

    static String parseStringField(String s, int k) {
        int[] r = fieldRange(s, k);
        return s.substring(r[0], r[1]);
    }

    // ---------------------------
    // Precompute distance via prefixN (count of N/n)
    // distance = (end-start) - Ncount in [start,end)
    // ---------------------------
    static int[] buildPrefixN(char[] seq) {
        int[] p = new int[seq.length + 1];
        int acc = 0;
        p[0] = 0;
        for (int i = 0; i < seq.length; i++) {
            char c = seq[i];
            if (c == 'N' || c == 'n') acc++;
            p[i + 1] = acc;
        }
        return p;
    }

    static int nonNDistance(int[] prefixN, int start, int end) {
        // assume valid bounds; caller ensures
        int nCount = prefixN[end] - prefixN[start];
        return (end - start) - nCount;
    }

    // (unused in main; kept for completeness)
    public static int existN(String str) {
        int dis = 0;
        for (int i = 0; i < str.length(); i++) {
            char c = str.charAt(i);
            if (c == 'N' || c == 'n') dis++;
        }
        return dis;
    }

    public static void main(String[] args) {
        if (args.length < 4) {
            System.err.println("Usage: java changeFormatPairwise_Optimized <genome.masked_all.fasta> <Genome_size_Indexes.txt> <SCN_LastzResult_...txt> <outdir>");
            return;
        }

        try {
            long startTime = System.currentTimeMillis();

            // args[0] genome.masked_all.fasta
            // args[1] Genome_size_Indexes.txt  (provides headers to read same count from fasta)
            ArrayList<String> chrs = new ArrayList<>();
            ArrayList<char[]> segs = new ArrayList<>();

            try (BufferedReader inIdx = new BufferedReader(new FileReader(args[1]), 1 << 20);
                 BufferedReader inFasta = new BufferedReader(new FileReader(args[0]), 1 << 20)) {

                boolean end = false;
                String str;
                while ((str = inIdx.readLine()) != null && !end) {
                    if (!str.isEmpty() && str.charAt(0) == '>') {
                        chrs.add(str.substring(1));
                        // fasta: header line + sequence line
                        inFasta.readLine(); // skip fasta header
                        String seq = inFasta.readLine();
                        if (seq == null) break;
                        segs.add(seq.toCharArray());
                    } else {
                        end = true;
                    }
                }
            }

            // Build prefixN arrays once for O(1) distance queries
            int[][] prefixNs = new int[segs.size()][];
            for (int i = 0; i < segs.size(); i++) {
                prefixNs[i] = buildPrefixN(segs.get(i));
            }

            // args[2] SCN_LastzResult_500NonCR_NewExtendPec50.txt
            // args[3] outdir
            String outdir = args[3];
            File outDirFile = new File(outdir);
            if (!outDirFile.exists() && !outDirFile.mkdirs()) {
                throw new IOException("Cannot create output directory: " + outdir);
            }

            int indexW = 0;
            int lineNo = 0;

            try (BufferedReader in = new BufferedReader(new FileReader(args[2]), 1 << 20);
                 BufferedWriter writer = new BufferedWriter(new FileWriter(new File(outDirFile, "Pairwise_SDs.txt")), 1 << 20);
                 BufferedWriter writerBG = new BufferedWriter(new FileWriter(new File(outDirFile, "BG_SDIndexes.fasta")), 1 << 20)) {

                // skip header
                in.readLine();

                writer.write("index  chrA  startA  endA  strandA  chrB  startB  endB  strandB  Match  Mismath  Indel  identity");
                writer.newLine();

                String line;
                while ((line = in.readLine()) != null) {
                    lineNo++;

                    int chr = parseIntField(line, 4);
                    int otherChr = parseIntField(line, 7);

                    // start1
                    int s16 = parseIntField(line, 16);
                    int start1 = (s16 == -1) ? parseIntField(line, 5) : s16;

                    // end1
                    int e17 = parseIntField(line, 17);
                    int end1 = (e17 == -1) ? parseIntField(line, 6) : e17;

                    // start2
                    int s18 = parseIntField(line, 18);
                    int start2 = (s18 == -1) ? parseIntField(line, 8) : s18;

                    // end2
                    int e19 = parseIntField(line, 19);
                    int end2 = (e19 == -1) ? parseIntField(line, 9) : e19;

                    // overlap filter (same as original)
                    if (!(chr != otherChr || start1 > end2 || start2 > end1)) {
                        continue;
                    }

                    double identity = parseDoubleField(line, 22);
                    if (identity < 0.7) {
                        continue;
                    }

                    // bounds safety (defensive): skip invalid
                    if (chr < 0 || chr >= segs.size() || otherChr < 0 || otherChr >= segs.size()) continue;
                    if (start1 < 0 || end1 > segs.get(chr).length || start1 >= end1) continue;
                    if (start2 < 0 || end2 > segs.get(otherChr).length || start2 >= end2) continue;

                    int length1 = nonNDistance(prefixNs[chr], start1, end1);
                    if (length1 < 500) continue;

                    int length2 = nonNDistance(prefixNs[otherChr], start2, end2);
                    if (length2 < 500) continue;

                    // strands and counts (keep as strings for writing)
                    String strandA = parseStringField(line, 14);
                    String strandB = parseStringField(line, 15);

                    String match = parseStringField(line, 10);
                    String mismatch = parseStringField(line, 11);
                    String indel = parseStringField(line, 13);

                    // write outputs (same format as original)
                    writer.write(indexW + "  " + chrs.get(chr) + "  " + start1 + "  " + end1 + "  " + strandA + "  "
                            + chrs.get(otherChr) + "  " + start2 + "  " + end2 + "  " + strandB + "  "
                            + match + "  " + mismatch + "  " + indel + "  " + identity);
                    writer.newLine();

                    writerBG.write(chr + "  " + start1 + "  " + end1);
                    writerBG.newLine();
                    writerBG.write(otherChr + "  " + start2 + "  " + end2);
                    writerBG.newLine();

                    indexW++;
                }
            }

            long endTime = System.currentTimeMillis();
            System.out.println("Pairwise SDs can be found in " + outdir + "/Pairwise_SDs.txt!");
            System.out.println("runTime:" + (endTime - startTime));

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
