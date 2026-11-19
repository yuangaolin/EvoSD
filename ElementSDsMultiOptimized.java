import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

public class ElementSDsMultiOptimized {

    // element record: start, end, midLen, index, multi
    static final class Elem {
        int s, e, midLen, idx, multi;
        Elem(int s, int e, int midLen, int idx) {
            this.s = s;
            this.e = e;
            this.midLen = midLen;
            this.idx = idx;
            this.multi = 0;
        }
    }

    // -------- fast line parsing (no regex) --------
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
        if (i < end && s.charAt(i) == '-') { sign = -1; i++; }
        int val = 0;
        for (; i < end; i++) val = val * 10 + (s.charAt(i) - '0');
        return val * sign;
    }

    static int getTokenInt(String line, int tokenIndex) {
        int i = 0, tok = 0, n = line.length();
        while (i < n) {
            i = nextNonSpace(line, i);
            if (i >= n) break;
            int j = nextSpace(line, i);
            if (tok == tokenIndex) return parseIntAt(line, i, j);
            tok++;
            i = j;
        }
        throw new IllegalArgumentException("Missing token[" + tokenIndex + "]: " + line);
    }

    // binary search: first element with s >= key
    static int lowerBoundByStart(Elem[] a, int key) {
        int lo = 0, hi = a.length;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (a[mid].s < key) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }

    // for query interval [qs, qe], increment multi for elements fully contained: s>=qs && e<=qe
    // Since array sorted by s, we start at first s>=qs, then scan forward until s>qe (can break)
    // But need e<=qe. This is O(K + tail) in worst-case; usually fast. If your intervals are huge
    // and elements dense, you can further index by end; but this is already a large improvement.
    static void addMultiplicity(Elem[] elemsSortedByStart, int qs, int qe) {
        if (elemsSortedByStart.length == 0) return;
        int pos = lowerBoundByStart(elemsSortedByStart, qs);
        for (int i = pos; i < elemsSortedByStart.length; i++) {
            Elem el = elemsSortedByStart[i];
            if (el.s > qe) break;           // no more possible
            if (el.e <= qe) el.multi++;     // fully contained
        }
    }

    public static void main(String[] args) {
        try {
            long startTime = System.currentTimeMillis();

            // -------- args[0] Genome_size_Indexes.txt -> GenomeSize (# of '>' lines until first non-header) --------
            int genomeSize = 0;
            try (BufferedReader br = new BufferedReader(new FileReader(args[0]), 1 << 20)) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isEmpty()) continue;
                    if (line.charAt(0) == '>') genomeSize++;
                    else break;
                }
            }

            // per chr elements
            @SuppressWarnings("unchecked")
            ArrayList<Elem>[] perChr = new ArrayList[genomeSize];
            for (int i = 0; i < genomeSize; i++) perChr[i] = new ArrayList<>();

            // -------- args[1] ElementSDs_100.fasta --------
            // header: skip line1
            // per your original: temp={onepair[5], onepair[6], onepair[4], index, 0}; pairs.get(onepair[1])
            int index = 0;
            try (BufferedReader br = new BufferedReader(new FileReader(args[1]), 1 << 20)) {
                br.readLine(); // header
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isEmpty()) continue;

                    int chr = getTokenInt(line, 1);
                    int s   = getTokenInt(line, 5);
                    int e   = getTokenInt(line, 6);
                    int mid = getTokenInt(line, 4);

                    perChr[chr].add(new Elem(s, e, mid, index));
                    index++;
                }
            }

            // Convert lists to arrays and sort by start for fast querying
            Elem[][] chrElems = new Elem[genomeSize][];
            for (int c = 0; c < genomeSize; c++) {
                ArrayList<Elem> list = perChr[c];
                Elem[] arr = list.toArray(new Elem[0]);
                Arrays.sort(arr, Comparator.comparingInt(a -> a.s));
                chrElems[c] = arr;
            }
            // allow GC of big lists
            for (int c = 0; c < genomeSize; c++) perChr[c].clear();

            // -------- args[2] SCN_LastzResult_500NonCR_NewExtendPec50.txt --------
            // For each line, update two intervals (chr1 and chr2) with optional extend overrides.
            try (BufferedReader br = new BufferedReader(new FileReader(args[2]), 1 << 20)) {
                br.readLine(); // header
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isEmpty()) continue;

                    // fields used (0-based):
                    // chr1=4 start1=5 end1=6 chr2=7 start2=8 end2=9
                    // extendS1=16 extendE1=17 extendS2=18 extendE2=19
                    int chr1 = getTokenInt(line, 4);
                    int s1   = getTokenInt(line, 5);
                    int e1   = getTokenInt(line, 6);
                    int exS1 = getTokenInt(line, 16);
                    int exE1 = getTokenInt(line, 17);
                    if (exS1 != -1) s1 = exS1;
                    if (exE1 != -1) e1 = exE1;
                    addMultiplicity(chrElems[chr1], s1, e1);

                    int chr2 = getTokenInt(line, 7);
                    int s2   = getTokenInt(line, 8);
                    int e2   = getTokenInt(line, 9);
                    int exS2 = getTokenInt(line, 18);
                    int exE2 = getTokenInt(line, 19);
                    if (exS2 != -1) s2 = exS2;
                    if (exE2 != -1) e2 = exE2;
                    addMultiplicity(chrElems[chr2], s2, e2);
                }
            }

            // -------- write output args[3] ElementSDs_LengthAndMulti.fasta --------
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File(args[3])), 1 << 20)) {
                writer.write("index   chr   start   end   midLength   multi");
                writer.newLine();
                for (int chr = 0; chr < genomeSize; chr++) {
                    Elem[] arr = chrElems[chr];
                    for (Elem el : arr) {
                        writer.write(el.idx + "  " + chr + "  " + el.s + "  " + el.e + "  " + el.midLen + "  " + el.multi);
                        writer.newLine();
                    }
                }
            }

            long endTime = System.currentTimeMillis();
            System.out.println("runTime: " + (endTime - startTime));

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
