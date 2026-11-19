import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

public class ElementSDsOptimized {

    // ---------------------------
    // Fast parsing (no regex)
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

    static int parseIntAt(String s, int start, int end) {
        int i = start;
        int sign = 1;
        if (i < end && s.charAt(i) == '-') { // safe
            sign = -1;
            i++;
        }
        int val = 0;
        for (; i < end; i++) val = val * 10 + (s.charAt(i) - '0');
        return val * sign;
    }

    // parse 3 ints: chr start end
    static void parse3Ints(String line, int[] out) {
        int i = 0;

        i = nextNonSpace(line, i);
        int j = nextSpace(line, i);
        out[0] = parseIntAt(line, i, j);

        i = nextNonSpace(line, j);
        j = nextSpace(line, i);
        out[1] = parseIntAt(line, i, j);

        i = nextNonSpace(line, j);
        j = nextSpace(line, i);
        out[2] = parseIntAt(line, i, j);
    }

    // ---------------------------
    // Prefix sum for non-'N'
    // prefix[k] = #nonN in [0, k)
    // distance inclusive [s..e] = prefix[e+1] - prefix[s]
    // ---------------------------
    static int[] buildPrefixNonN(char[] cr) {
        int[] prefix = new int[cr.length + 1];
        int acc = 0;
        prefix[0] = 0;
        for (int i = 0; i < cr.length; i++) {
            if (cr[i] != 'N') acc++;
            prefix[i + 1] = acc;
        }
        return prefix;
    }

    static int distInclusive(int[] prefixNonN, int s, int eInclusive) {
        if (eInclusive < s) return 0;
        return prefixNonN[eInclusive + 1] - prefixNonN[s];
    }

    // ---------------------------
    // int dynamic array to avoid Integer boxing
    // ---------------------------
    static final class IntList {
        int[] a;
        int size;

        IntList(int cap) {
            a = new int[Math.max(16, cap)];
            size = 0;
        }

        void clear() { size = 0; }

        void add(int v) {
            if (size == a.length) {
                int[] b = new int[a.length << 1];
                System.arraycopy(a, 0, b, 0, a.length);
                a = b;
            }
            a[size++] = v;
        }

        int get(int idx) { return a[idx]; }
    }

    public static void main(String[] args) {
        try {
            long startTime = System.currentTimeMillis();


            ArrayList<char[]> CRs = new ArrayList<>();
            ArrayList<int[]> prefixNonN = new ArrayList<>();
            ArrayList<char[]> ends = new ArrayList<>();
            ArrayList<boolean[]> lines = new ArrayList<>();

            // Read genome + endpoints together
            try (BufferedReader brChr = new BufferedReader(new FileReader(args[0]), 1 << 20);
                 BufferedReader brEnds = new BufferedReader(new FileReader(args[1]), 1 << 20)) {

                while (true) {
                    String h1 = brChr.readLine();
                    if (h1 == null) break;
                    String seq = brChr.readLine();
                    if (seq == null) break;

                    String h2 = brEnds.readLine();
                    if (h2 == null) break;
                    String endLine = brEnds.readLine();
                    if (endLine == null) break;

                    char[] cr = seq.toCharArray();
                    CRs.add(cr);
                    prefixNonN.add(buildPrefixNonN(cr));

                    ends.add(endLine.toCharArray());

                    lines.add(new boolean[cr.length]); // default false
                }
            }

            final int clusterLength = 100;
            int badNum = 0;
            int badInterval = 0;

            String outdir = args[3];
            File outDir = new File(outdir);

            try (BufferedReader br = new BufferedReader(new FileReader(args[2]), 1 << 20);
                 BufferedWriter writer = new BufferedWriter(new FileWriter(new File(outDir, "ElementSDs_100.fasta")), 1 << 20);
                 BufferedWriter writerC = new BufferedWriter(new FileWriter(new File(outDir, "Clusters_100.fasta")), 1 << 20);
                 BufferedWriter writerB = new BufferedWriter(new FileWriter(new File(outDir, "BadMosaicSD_100.fasta")), 1 << 20)) {

                writerB.write("chr  start  end  points  cluster  length  compactL");
                writerB.newLine();
                writer.write("index  chr  midStart  midEnd  midLength  startC  endC  lengthC");
                writer.newLine();
                writerC.write("index   chr   start   end   size   length  compactLength");
                writerC.newLine();

                // skip header of BG_MosaicSDs.fasta
                br.readLine();

                int[] triple = new int[3];
                int index = 0, indexC = 0;

                IntList points = new IntList(1024);
                IntList cluster = new IntList(1024);

                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isEmpty()) continue;

                    parse3Ints(line, triple);
                    int chr = triple[0];
                    int start = triple[1];
                    int end = triple[2];

                    char[] endsChr = ends.get(chr);
                    int[] pref = prefixNonN.get(chr);
                    boolean[] lineMask = lines.get(chr);

                    // collect points in [start, end)
                    points.clear();
                    for (int i = start; i < end; i++) {
                        if (endsChr[i] == '1') points.add(i);
                    }

                    // Original code would crash if no points; keep robust and count as bad
                    if (points.size == 0) {
                        badNum++;
                        writerB.write(chr + "  " + start + "  " + end + "  0  -1  " + (end - start) + "  "
                                + distInclusive(pref, start, end - 1));
                        writerB.newLine();
                        continue;
                    }

                    // clustering
                    cluster.clear();
                    int cur = points.get(0);
                    cluster.add(cur);

                    int lastPoint = 0;
                    int mPoint = -1;

                    for (int i = 1; i < points.size; i++) {
                        int pi = points.get(i);

                        // distance(cur, pi) in non-N bases (inclusive)
                        if (distInclusive(pref, cur, pi) > clusterLength) {
                            int cStart = cluster.get(0);
                            int cEnd = cur;
                            int cSize = cluster.size;
                            int spanLen = cEnd - cStart + 1;
                            int compactLen = distInclusive(pref, cStart, cEnd);

                            writerC.write(indexC + "   " + chr + "   " + cStart + "   " + cEnd + "   " + cSize
                                    + "   " + spanLen + "  " + compactLen);
                            writerC.newLine();
                            indexC++;

                            for (int k = cStart; k <= cEnd; k++) lineMask[k] = true;

                            if (mPoint != -1) {
                                int mid = cluster.get(cluster.size / 2);

                                int midLen = distInclusive(pref, mPoint, mid);
                                int lenC = distInclusive(pref, lastPoint, cStart);

                                writer.write(index + "  " + chr + "  " + mPoint + "  " + (mid + 1) + "  " + midLen
                                        + "  " + lastPoint + "  " + (cStart + 1) + "  " + lenC);
                                writer.newLine();
                                index++;
                            }

                            lastPoint = cur;
                            mPoint = cluster.get(cluster.size / 2);
                            cluster.clear();
                        }

                        cur = pi;
                        cluster.add(cur);
                    }

                    // flush last cluster
                    int cStart = cluster.get(0);
                    int cEnd = cur;
                    int cSize = cluster.size;
                    int spanLen = cEnd - cStart + 1;
                    int compactLen = distInclusive(pref, cStart, cEnd);

                    writerC.write(indexC + "   " + chr + "   " + cStart + "   " + cEnd + "   " + cSize
                            + "   " + spanLen + "  " + compactLen);
                    writerC.newLine();
                    indexC++;

                    for (int k = cStart; k <= cEnd; k++) lineMask[k] = true;

                    if (mPoint != -1) {
                        int mid = cluster.get(cluster.size / 2);

                        int midLen = distInclusive(pref, mPoint, mid);
                        int lenC = distInclusive(pref, lastPoint, cStart);

                        writer.write(index + "  " + chr + "  " + mPoint + "  " + (mid + 1) + "  " + midLen
                                + "  " + lastPoint + "  " + (cStart + 1) + "  " + lenC);
                        writer.newLine();
                        index++;
                    } else {
                        // bad mosaic SD
                        badNum++;
                        badInterval += points.size / 2;

                        writerB.write(chr + "  " + start + "  " + end + "  " + points.size + "  " + (indexC - 1)
                                + "  " + (end - start) + "  " + distInclusive(pref, start, end - 1));
                        writerB.newLine();
                    }
                }
            }

            // output ClustersLine.txt
            try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outdir, "ClustersLine.txt")), 1 << 20)) {
                for (int n = 0; n < lines.size(); n++) {
                    bw.write(">" + n);
                    bw.newLine();
                    boolean[] arr = lines.get(n);
                    for (int i = 0; i < arr.length; i++) bw.write(arr[i] ? '1' : '0');
                    bw.newLine();
                }
            }

            long endTime = System.currentTimeMillis();
            System.out.println("end! badNum=" + badNum + " badInterval=" + badInterval);
            System.out.println("runTime: " + (endTime - startTime));

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
