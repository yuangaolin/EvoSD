import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;

public class NewElementSDsPairOptimized {
    // ---------------------------
    // Fast IO parsing helpers
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
        if (i < end && s.charAt(i) == '-') {
            sign = -1;
            i++;
        }
        int val = 0;
        for (; i < end; i++) val = val * 10 + (s.charAt(i) - '0');
        return val * sign;
    }

    // parse k-th token (0-based) as int, no split
    static int parseIntToken(String line, int tokenIdx) {
        int i = 0;
        int cur = 0;
        int n = line.length();
        while (i < n) {
            i = nextNonSpace(line, i);
            if (i >= n) break;
            int j = nextSpace(line, i);
            if (cur == tokenIdx) return parseIntAt(line, i, j);
            cur++;
            i = j;
        }
        throw new IllegalArgumentException("Not enough tokens: need " + tokenIdx + " in line: " + line);
    }

    static String parseStringToken(String line, int tokenIdx) {
        int i = 0;
        int cur = 0;
        int n = line.length();
        while (i < n) {
            i = nextNonSpace(line, i);
            if (i >= n) break;
            int j = nextSpace(line, i);
            if (cur == tokenIdx) return line.substring(i, j);
            cur++;
            i = j;
        }
        throw new IllegalArgumentException("Not enough tokens: need " + tokenIdx + " in line: " + line);
    }

    // ---------------------------
    // Reverse complement (same behavior as your reverse)
    // ---------------------------
    static final char[] COMP = new char[128];
    static {
        for (int i = 0; i < COMP.length; i++) COMP[i] = (char) i;
        COMP['A'] = 'T'; COMP['C'] = 'G'; COMP['G'] = 'C'; COMP['T'] = 'A'; COMP['N'] = 'N';
        COMP['a'] = 't'; COMP['c'] = 'g'; COMP['g'] = 'c'; COMP['t'] = 'a'; COMP['n'] = 'n';
    }

    static String reverse(String str) {
        int L = str.length();
        char[] out = new char[L];
        for (int i = 0; i < L; i++) {
            char b = str.charAt(L - 1 - i);
            out[i] = b < 128 ? COMP[b] : b;
        }
        return new String(out);
    }

    // ---------------------------
    // CIGAR alignment (optimized: no substring/parseInt per op)
    // Keeps your behavior: skips 'N' in seg1/seg2 when consuming bases
    // ---------------------------
    static String[] cigarAlignment(String seg1, String seg2, String cigar) {
        // first pass: total alignment length
        int total = 0;
        int num = 0;
        for (int i = 0; i < cigar.length(); i++) {
            char c = cigar.charAt(i);
            if (c >= '0' && c <= '9') {
                num = num * 10 + (c - '0');
            } else {
                // M/I/D
                total += num;
                num = 0;
            }
        }

        char[] a1 = new char[total];
        char[] a2 = new char[total];

        int idx1 = 0, idx2 = 0, out = 0;
        num = 0;
        for (int i = 0; i < cigar.length(); i++) {
            char op = cigar.charAt(i);
            if (op >= '0' && op <= '9') {
                num = num * 10 + (op - '0');
                continue;
            }

            if (op == 'M') {
                for (int k = 0; k < num; k++) {
                    while (seg1.charAt(idx1) == 'N') idx1++;
                    while (seg2.charAt(idx2) == 'N') idx2++;
                    a1[out] = seg1.charAt(idx1++);
                    a2[out] = seg2.charAt(idx2++);
                    out++;
                }
            } else if (op == 'D') {
                for (int k = 0; k < num; k++) {
                    while (seg1.charAt(idx1) == 'N') idx1++;
                    a1[out] = seg1.charAt(idx1++);
                    a2[out] = '-';
                    out++;
                }
            } else if (op == 'I') {
                for (int k = 0; k < num; k++) {
                    while (seg2.charAt(idx2) == 'N') idx2++;
                    a1[out] = '-';
                    a2[out] = seg2.charAt(idx2++);
                    out++;
                }
            }
            num = 0;
        }

        return new String[]{ new String(a1), new String(a2) };
    }

    // ---------------------------
    // computeIndex: reduced allocations (no int[2][L])
    // Still follows your original semantics incl. gap fill.
    // Note: this version assumes alignment strings don't contain 'N' (true given cigarAlignment skips N).
    // ---------------------------
    static void computeIndex(String[] alignment,
                             int[] seg1To2, int[] seg2To1,
                             String seg1, String seg2,
                             String str1, String str2) {

        final String al1 = alignment[0];
        final String al2 = alignment[1];
        final int L = al1.length();

        int[] colToIdx1 = new int[L];
        int[] colToIdx2 = new int[L];

        // Map alignment columns -> seg1 index (skipping seg 'N' in the same way as original)
        if (str1.equalsIgnoreCase("+")) {
            int idx = 0;
            while (idx < seg1.length() && seg1.charAt(idx) == 'N') idx++;
            int cur = idx;

            // leading gaps
            int p = 0;
            while (p < L && al1.charAt(p) == '-') colToIdx1[p++] = cur;
            if (p < L) colToIdx1[p] = cur;

            for (int i = p + 1; i < L; i++) {
                if (al1.charAt(i) != '-') {
                    cur++;
                    while (cur < seg1.length() && seg1.charAt(cur) == 'N') cur++;
                }
                colToIdx1[i] = cur;
            }
        } else {
            int idx = 0;
            while (idx < seg1.length() && seg1.charAt(idx) == 'N') idx++;
            int cur = seg1.length() - 1 - idx;

            int p = 0;
            while (p < L && al1.charAt(p) == '-') colToIdx1[p++] = cur;
            if (p < L) colToIdx1[p] = cur;

            // move forward in "consumed non-N bases" means decreasing index on reverse strand
            int consumed = idx;
            for (int i = p + 1; i < L; i++) {
                if (al1.charAt(i) != '-') {
                    consumed++;
                    while (consumed < seg1.length() && seg1.charAt(consumed) == 'N') consumed++;
                    cur = seg1.length() - 1 - consumed;
                }
                colToIdx1[i] = cur;
            }
        }

        // Map alignment columns -> seg2 index
        if (str2.equalsIgnoreCase("+")) {
            int idx = 0;
            while (idx < seg2.length() && seg2.charAt(idx) == 'N') idx++;
            int cur = idx;

            int p = 0;
            while (p < L && al2.charAt(p) == '-') colToIdx2[p++] = cur;
            if (p < L) colToIdx2[p] = cur;

            for (int i = p + 1; i < L; i++) {
                if (al2.charAt(i) != '-') {
                    cur++;
                    while (cur < seg2.length() && seg2.charAt(cur) == 'N') cur++;
                }
                colToIdx2[i] = cur;
            }
        } else {
            int idx = 0;
            while (idx < seg2.length() && seg2.charAt(idx) == 'N') idx++;
            int cur = seg2.length() - 1 - idx;

            int p = 0;
            while (p < L && al2.charAt(p) == '-') colToIdx2[p++] = cur;
            if (p < L) colToIdx2[p] = cur;

            int consumed = idx;
            for (int i = p + 1; i < L; i++) {
                if (al2.charAt(i) != '-') {
                    consumed++;
                    while (consumed < seg2.length() && seg2.charAt(consumed) == 'N') consumed++;
                    cur = seg2.length() - 1 - consumed;
                }
                colToIdx2[i] = cur;
            }
        }

        // Fill direct mappings (non-gap columns)
        for (int i = 0; i < L; i++) {
            if (al1.charAt(i) != '-') seg1To2[colToIdx1[i]] = colToIdx2[i];
            if (al2.charAt(i) != '-') seg2To1[colToIdx2[i]] = colToIdx1[i];
        }

        // gap fill seg1To2
        int idx = 0, start = 0;
        while (idx < seg1To2.length && seg1To2[idx] == -1) idx++;
        if (idx < seg1To2.length) {
            for (int i = start; i < idx; i++) seg1To2[i] = seg1To2[idx];
            idx++; start = idx;
            while (idx < seg1To2.length) {
                while (idx < seg1To2.length && seg1To2[idx] == -1) idx++;
                for (int i = start; i < idx; i++) seg1To2[i] = seg1To2[start - 1];
                idx++; start = idx;
            }
        }

        // gap fill seg2To1
        idx = 0; start = 0;
        while (idx < seg2To1.length && seg2To1[idx] == -1) idx++;
        if (idx < seg2To1.length) {
            for (int i = start; i < idx; i++) seg2To1[i] = seg2To1[idx];
            idx++; start = idx;
            while (idx < seg2To1.length) {
                while (idx < seg2To1.length && seg2To1[idx] == -1) idx++;
                for (int i = start; i < idx; i++) seg2To1[i] = seg2To1[start - 1];
                idx++; start = idx;
            }
        }
    }

    // ---------------------------
    // exist(): keep your logic, but avoid repeated charAt via local refs
    // ---------------------------
    static boolean exist(char[] clusters2, String chr1, int[] seg1to2,
                         int segStart, int segEnd, int start1, int start2) {

        int index1 = segStart;
        int index2 = segEnd - 1;

        while (index1 != index2) {
            index1++;
            while (chr1.charAt(index1) == 'N') index1++;
            if (index1 != index2) {
                index2--;
                while (chr1.charAt(index2) == 'N') index2--;
            }
        }

        int mapped = seg1to2[index1 - start1];
        return clusters2[start2 + mapped] != '1';
    }

    // ---------------------------
    // Per-chr pairs storage: sort by start, binary search range
    // ---------------------------
    static final class PairBucket {
        int[] s = new int[1024];
        int[] e = new int[1024];
        int[] id = new int[1024];
        int size = 0;

        void add(int start, int end, int idx) {
            if (size == s.length) {
                int n = s.length << 1;
                s = Arrays.copyOf(s, n);
                e = Arrays.copyOf(e, n);
                id = Arrays.copyOf(id, n);
            }
            s[size] = start;
            e[size] = end;
            id[size] = idx;
            size++;
        }

        void sortByStart() {
            // sort indices by s; in-place quicksort on parallel arrays
            if (size <= 1) return;
            quickSort(0, size - 1);
        }

        void quickSort(int lo, int hi) {
            if (lo >= hi) return;
            int i = lo, j = hi;
            int pivot = s[(lo + hi) >>> 1];
            while (i <= j) {
                while (s[i] < pivot) i++;
                while (s[j] > pivot) j--;
                if (i <= j) {
                    swap(i, j);
                    i++; j--;
                }
            }
            if (lo < j) quickSort(lo, j);
            if (i < hi) quickSort(i, hi);
        }

        void swap(int i, int j) {
            int ts = s[i]; s[i] = s[j]; s[j] = ts;
            int te = e[i]; e[i] = e[j]; e[j] = te;
            int tid = id[i]; id[i] = id[j]; id[j] = tid;
        }

        // first index with s[idx] >= key
        int lowerBoundStart(int key) {
            int lo = 0, hi = size;
            while (lo < hi) {
                int mid = (lo + hi) >>> 1;
                if (s[mid] < key) lo = mid + 1;
                else hi = mid;
            }
            return lo;
        }
    }

    public static void main(String[] args) {
        try {
            long startTime = System.currentTimeMillis();

            // args[0] chroms (masked)
            // args[1] chromsN (unmasked or N-handling source)
            // args[2] ClustersLine.txt
            // args[3] ElementSDs_... (pairs source)  (used to build pairs buckets)
            // args[4] length file (one line per record)
            // args[5] SCN_LastzResult... (main align file)
            // args[6] outdir

            ArrayList<String> chroms = new ArrayList<>();
            ArrayList<String> chromsN = new ArrayList<>();
            ArrayList<char[]> clusters = new ArrayList<>();

            try (BufferedReader brChr = new BufferedReader(new FileReader(args[0]), 1 << 20);
                 BufferedReader brChrN = new BufferedReader(new FileReader(args[1]), 1 << 20);
                 BufferedReader brClu = new BufferedReader(new FileReader(args[2]), 1 << 20)) {

                while (true) {
                    String h = brChr.readLine();
                    if (h == null) break;
                    String seq = brChr.readLine();
                    if (seq == null) break;

                    String hN = brChrN.readLine();
                    if (hN == null) break;
                    String seqN = brChrN.readLine();
                    if (seqN == null) break;

                    String hC = brClu.readLine();
                    if (hC == null) break;
                    String clu = brClu.readLine();
                    if (clu == null) break;

                    chroms.add(seq);
                    chromsN.add(seqN);
                    clusters.add(clu.toCharArray());
                }
            }

            System.out.println("end reading " + clusters.size());

            // Build pairs buckets per chr
            PairBucket[] buckets = new PairBucket[chroms.size()];
            for (int i = 0; i < buckets.length; i++) buckets[i] = new PairBucket();

            try (BufferedReader brPairs = new BufferedReader(new FileReader(args[3]), 1 << 20)) {
                brPairs.readLine(); // skip header
                String line;
                int idx = 1;
                while ((line = brPairs.readLine()) != null) {
                    if (line.isEmpty()) continue;

                    // original: onepair[1] chr, onepair[5] start, onepair[6] end
                    int chr = parseIntToken(line, 1);
                    int s = parseIntToken(line, 5);
                    int e = parseIntToken(line, 6);
                    buckets[chr].add(s, e, idx);
                    idx++;
                }
                System.out.println(idx);
            }

            for (PairBucket b : buckets) b.sortByStart();

            // open length file + main align file
            String outdir = args[6];
            File outDir = new File(outdir);

            try (BufferedReader brLen = new BufferedReader(new FileReader(args[4]), 1 << 20);
                 BufferedReader br = new BufferedReader(new FileReader(args[5]), 1 << 20);
                 BufferedWriter writerE = new BufferedWriter(new FileWriter(new File(outDir, "ElementSDs_pairwiseEqual.fasta")), 1 << 20);
                 BufferedWriter writerU = new BufferedWriter(new FileWriter(new File(outDir, "ElementSDs_pairwiseUnequal.fasta")), 1 << 20);
                 BufferedWriter writerLE = new BufferedWriter(new FileWriter(new File(outDir, "ElementSDs_pairwiseLengthEqual.fasta")), 1 << 20);
                 BufferedWriter writerLU = new BufferedWriter(new FileWriter(new File(outDir, "ElementSDs_pairwiseLengthUequal.fasta")), 1 << 20);
                 BufferedWriter writerDistri = new BufferedWriter(new FileWriter(new File(outDir, "ElementSDs_pairwiseDistri.fasta")), 1 << 20);
                 BufferedWriter writerSpecial = new BufferedWriter(new FileWriter(new File(outDir, "ElementSDs_pairwiseSpecial_9.fasta")), 1 << 20);
                 BufferedWriter writerDif = new BufferedWriter(new FileWriter(new File(outDir, "ElementSDs_pairwiseDifference.fasta")), 1 << 20)) {

                int[][] distri = new int[12][12];

                br.readLine(); // skip header
                int good = 0, bad = 0, wrong = 0, miss = 0;
                int index = -1;

                // reuse buffers to reduce GC
                ArrayList<Integer> num1 = new ArrayList<>(256);
                ArrayList<Integer> num2 = new ArrayList<>(256);
                ArrayList<Integer> origin1 = new ArrayList<>(256);
                ArrayList<Integer> origin2 = new ArrayList<>(256);

                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isEmpty()) continue;
                    index++;

                    String lengthLine = brLen.readLine(); // one-to-one
                    if (lengthLine == null) break;

                    // Tokens used from original:
                    // 1 cigar, 4 chr,5 start,6 end,7 otherChr,8 otherStart,9 otherEnd, 14 strand1,15 strand2, 16-19 ext, 20 cigarS,21 cigarE,22 identity
                    String cigarMain = parseStringToken(line, 1);

                    int chr = parseIntToken(line, 4);
                    int start = parseIntToken(line, 5);
                    int end = parseIntToken(line, 6);

                    int otherChr = parseIntToken(line, 7);
                    int otherStart = parseIntToken(line, 8);
                    int otherEnd = parseIntToken(line, 9);

                    String strand1 = parseStringToken(line, 14);
                    String strand2 = parseStringToken(line, 15);

                    int extendS1 = parseIntToken(line, 16);
                    int extendE1 = parseIntToken(line, 17);
                    int extendS2 = parseIntToken(line, 18);
                    int extendE2 = parseIntToken(line, 19);

                    String cigarS = parseStringToken(line, 20);
                    String cigarE = parseStringToken(line, 21);

                    String identity = parseStringToken(line, 22);

                    // build segs
                    String seg1, seg2, segS1 = "", segS2 = "", segE1 = "", segE2 = "";

                    if (strand1.equalsIgnoreCase("+")) {
                        seg1 = chroms.get(chr).substring(start, end);
                        if (!cigarS.equalsIgnoreCase("-1")) {
                            segS1 = chromsN.get(chr).substring(extendS1, start).toUpperCase();
                            start = extendS1;
                        }
                        if (!cigarE.equalsIgnoreCase("-1")) {
                            segE1 = chromsN.get(chr).substring(end, extendE1).toUpperCase();
                            end = extendE1;
                        }
                    } else {
                        seg1 = reverse(chroms.get(chr).substring(start, end));
                        if (!cigarS.equalsIgnoreCase("-1")) {
                            segS1 = reverse(chromsN.get(chr).substring(end, extendE1).toUpperCase());
                            end = extendE1;
                        }
                        if (!cigarE.equalsIgnoreCase("-1")) {
                            segE1 = reverse(chromsN.get(chr).substring(extendS1, start).toUpperCase());
                            start = extendS1;
                        }
                    }

                    if (strand2.equalsIgnoreCase("+")) {
                        seg2 = chroms.get(otherChr).substring(otherStart, otherEnd);
                        if (!cigarS.equalsIgnoreCase("-1")) {
                            segS2 = chromsN.get(otherChr).substring(extendS2, otherStart).toUpperCase();
                            otherStart = extendS2;
                        }
                        if (!cigarE.equalsIgnoreCase("-1")) {
                            segE2 = chromsN.get(otherChr).substring(otherEnd, extendE2).toUpperCase();
                            otherEnd = extendE2;
                        }
                    } else {
                        seg2 = reverse(chroms.get(otherChr).substring(otherStart, otherEnd));
                        if (!cigarS.equalsIgnoreCase("-1")) {
                            segS2 = reverse(chromsN.get(otherChr).substring(otherEnd, extendE2).toUpperCase());
                            otherEnd = extendE2;
                        }
                        if (!cigarE.equalsIgnoreCase("-1")) {
                            segE2 = reverse(chromsN.get(otherChr).substring(extendS2, otherStart).toUpperCase());
                            otherStart = extendS2;
                        }
                    }

                    // alignment main
                    String[] alignment = cigarAlignment(seg1, seg2, cigarMain);

                    // prefix/suffix extension alignments (avoid repeated string concatenations with builders)
                    if (!cigarS.equalsIgnoreCase("-1")) {
                        String[] aS = cigarAlignment(segS1, segS2, cigarS);
                        StringBuilder b0 = new StringBuilder(aS[0].length() + alignment[0].length());
                        StringBuilder b1 = new StringBuilder(aS[1].length() + alignment[1].length());
                        b0.append(aS[0]).append(alignment[0]);
                        b1.append(aS[1]).append(alignment[1]);
                        alignment[0] = b0.toString();
                        alignment[1] = b1.toString();
                        seg1 = segS1 + seg1;
                        seg2 = segS2 + seg2;
                    }
                    if (!cigarE.equalsIgnoreCase("-1")) {
                        String[] aE = cigarAlignment(segE1, segE2, cigarE);
                        StringBuilder b0 = new StringBuilder(alignment[0].length() + aE[0].length());
                        StringBuilder b1 = new StringBuilder(alignment[1].length() + aE[1].length());
                        b0.append(alignment[0]).append(aE[0]);
                        b1.append(alignment[1]).append(aE[1]);
                        alignment[0] = b0.toString();
                        alignment[1] = b1.toString();
                        seg1 = seg1 + segE1;
                        seg2 = seg2 + segE2;
                    }

                    // build mapping arrays
                    int[] seg1to2 = new int[seg1.length()];
                    int[] seg2to1 = new int[seg2.length()];
                    Arrays.fill(seg1to2, -1);
                    Arrays.fill(seg2to1, -1);

                    computeIndex(alignment, seg1to2, seg2to1, seg1, seg2, strand1, strand2);

                    // query pairs via buckets (binary search on start)
                    num1.clear(); num2.clear(); origin1.clear(); origin2.clear();

                    // chr side
                    PairBucket b1 = buckets[chr];
                    int p = b1.lowerBoundStart(start);
                    for (int i = p; i < b1.size; i++) {
                        int ps = b1.s[i];
                        if (ps > end) break;
                        int pe = b1.e[i];
                        if (ps >= start && pe <= end) {
                            int pid = b1.id[i];
                            origin1.add(pid);
                            if (exist(clusters.get(otherChr), chroms.get(chr), seg1to2, ps, pe, start, otherStart)) {
                                num1.add(pid);
                            }
                        }
                    }

                    // otherChr side
                    PairBucket b2 = buckets[otherChr];
                    p = b2.lowerBoundStart(otherStart);
                    for (int i = p; i < b2.size; i++) {
                        int ps = b2.s[i];
                        if (ps > otherEnd) break;
                        int pe = b2.e[i];
                        if (ps >= otherStart && pe <= otherEnd) {
                            int pid = b2.id[i];
                            origin2.add(pid);
                            if (exist(clusters.get(chr), chroms.get(otherChr), seg2to1, ps, pe, otherStart, start)) {
                                num2.add(pid);
                            }
                        }
                    }

                    if (origin1.size() == 0 || origin2.size() == 0) miss++;

                    if (num1.size() == 0 || num2.size() == 0) {
                        wrong++;
                        if (origin1.size() > 0 && origin2.size() > 0) {
                            writerSpecial.write(index + ":  ");

                            if (strand1.equalsIgnoreCase("+")) {
                                for (int v : origin1) writerSpecial.write(v + "  ");
                            } else {
                                for (int i = origin1.size() - 1; i >= 0; i--) writerSpecial.write(origin1.get(i) + "  ");
                            }

                            writerSpecial.write("--  ");

                            if (strand2.equalsIgnoreCase("+")) {
                                for (int v : origin2) writerSpecial.write(v + "  ");
                            } else {
                                for (int i = origin2.size() - 1; i >= 0; i--) writerSpecial.write(origin2.get(i) + "  ");
                            }

                            writerSpecial.newLine();
                            writerSpecial.newLine();
                        }
                    } else {
                        writerDif.write(index + "  " + origin1.size() + "  " + origin2.size() + "  " + num1.size() + "  " + num2.size());
                        writerDif.newLine();

                        if (num1.size() == num2.size()) {
                            good++;
                            writerE.write(index + ":  ");

                            if (strand1.equalsIgnoreCase("+")) {
                                for (int v : num1) writerE.write(v + "  ");
                            } else {
                                for (int i = num1.size() - 1; i >= 0; i--) writerE.write("-" + num1.get(i) + "  ");
                            }

                            writerE.write("--  ");

                            if (strand2.equalsIgnoreCase("+")) {
                                for (int v : num2) writerE.write(v + "  ");
                            } else {
                                for (int i = num2.size() - 1; i >= 0; i--) writerE.write("-" + num2.get(i) + "  ");
                            }

                            writerE.write(identity);
                            writerE.newLine();
                            writerE.newLine();

                            writerLE.write(lengthLine);
                            writerLE.newLine();
                        } else {
                            bad++;
                            writerU.write(index + ":  ");

                            if (strand1.equalsIgnoreCase("+")) {
                                for (int v : num1) writerU.write(v + "  ");
                            } else {
                                for (int i = num1.size() - 1; i >= 0; i--) writerU.write(num1.get(i) + "  ");
                            }

                            writerU.write("--  ");

                            if (strand2.equalsIgnoreCase("+")) {
                                for (int v : num2) writerU.write(v + "  ");
                            } else {
                                for (int i = num2.size() - 1; i >= 0; i--) writerU.write(num2.get(i) + "  ");
                            }

                            writerU.newLine();
                            writerU.newLine();

                            writerLU.write(lengthLine);
                            writerLU.newLine();
                        }

                        int idx1;
                        if (num1.size() <= 10) idx1 = num1.size() - 1;
                        else if (num1.size() <= 20) idx1 = 10;
                        else idx1 = 11;

                        int idx2;
                        if (num2.size() <= 10) idx2 = num2.size() - 1;
                        else if (num2.size() <= 20) idx2 = 10;
                        else idx2 = 11;

                        distri[idx1][idx2]++;
                    }
                }

                for (int i = 0; i < 12; i++) {
                    for (int j = 0; j < 12; j++) {
                        writerDistri.write(distri[i][j] + "  ");
                    }
                    writerDistri.newLine();
                }

                long endTime = System.currentTimeMillis();
                System.out.println("end! " + index + "  " + good + "  " + bad + "  " + wrong + "  " + miss);
                System.out.println("runTime:  " + (endTime - startTime));
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
