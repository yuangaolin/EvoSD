import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;

public class SDblockOptimized {

    // ---------------------------
    // Fast token parsing helpers
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

    static int prevNonSpace(String s, int i) {
        while (i >= 0 && s.charAt(i) <= ' ') i--;
        return i;
    }

    static int prevSpace(String s, int i) {
        while (i >= 0 && s.charAt(i) > ' ') i--;
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

    static boolean isDashDashToken(String s, int start, int end) {
        return (end - start) == 2 && s.charAt(start) == '-' && s.charAt(start + 1) == '-';
    }

    // ---------------------------
    // DSU with parity
    // parity[x] = parity from x to parent[x]
    // constraint: parity(x) XOR parity(y) = w (0 same sign, 1 opposite sign)
    // ---------------------------
    static final class DSU {
        final int[] parent;
        final byte[] parity;
        final int[] size;

        DSU(int n) {
            parent = new int[n + 1];
            parity = new byte[n + 1];
            size = new int[n + 1];
            for (int i = 1; i <= n; i++) {
                parent[i] = i;
                parity[i] = 0;
                size[i] = 1;
            }
        }

        int find(int x) {
            if (parent[x] == x) return x;
            int p = parent[x];
            int r = find(p);
            parity[x] ^= parity[p];
            parent[x] = r;
            return r;
        }

        byte parityToRoot(int x) {
            find(x);
            return parity[x];
        }

        void union(int x, int y, int w) {
            int rx = find(x);
            int ry = find(y);
            byte px = parity[x];
            byte py = parity[y];
            if (rx == ry) return;

            // union by size
            if (size[rx] < size[ry]) {
                int tmp = rx; rx = ry; ry = tmp;
                byte t = px; px = py; py = t;
            }

            parent[ry] = rx;
            parity[ry] = (byte) (px ^ py ^ (w & 1));
            size[rx] += size[ry];
        }
    }

    // ---------------------------
    // int list (avoid boxing)
    // ---------------------------
    static final class IntList {
        int[] a;
        int size;

        IntList(int cap) {
            a = new int[Math.max(16, cap)];
        }

        void add(int v) {
            if (size == a.length) a = Arrays.copyOf(a, a.length << 1);
            a[size++] = v;
        }
    }

    public static void main(String[] args) {
        try {
            long startTime = System.currentTimeMillis();

            // args[0] lengths file (skip header, use token[4])
            // args[1] pairwiseEqual file: "0:  3  --  1  0.8685" (last token is identity float)
            // args[2] outdir

            // -------- read lengths --------
            int[] lengths;
            int n;
            try (BufferedReader br = new BufferedReader(new FileReader(args[0]), 1 << 20)) {
                br.readLine(); // header
                IntList lens = new IntList(1 << 20);
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isEmpty()) continue;

                    int i = 0, token = 0;
                    int L = line.length();
                    int val = 0;
                    boolean got = false;

                    while (i < L) {
                        i = nextNonSpace(line, i);
                        if (i >= L) break;
                        int j = nextSpace(line, i);
                        if (token == 4) {
                            val = parseIntAt(line, i, j);
                            got = true;
                            break;
                        }
                        token++;
                        i = j;
                    }
                    if (!got) throw new IllegalArgumentException("Cannot parse length token[4]: " + line);
                    lens.add(val);
                }
                n = lens.size;
                lengths = Arrays.copyOf(lens.a, n);
            }

            DSU dsu = new DSU(n);

            // -------- read pairwiseEqual and union edges --------
            try (BufferedReader br = new BufferedReader(new FileReader(args[1]), 1 << 20)) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isEmpty()) continue;

                    // cut off the last token (identity)
                    int end = prevNonSpace(line, line.length() - 1);
                    if (end < 0) continue;
                    int idTokStart = prevSpace(line, end) + 1;
                    int parseEnd = idTokStart; // parse [0, parseEnd)

                    // skip first token "idx:"
                    int i = 0;
                    i = nextNonSpace(line, i);
                    if (i >= parseEnd) continue;
                    int j = nextSpace(line, i);
                    i = j;

                    IntList left = new IntList(16);
                    IntList right = new IntList(16);

                    boolean sawSep = false;
                    // left ints until "--"
                    while (true) {
                        i = nextNonSpace(line, i);
                        if (i >= parseEnd) break;
                        j = nextSpace(line, i);
                        if (isDashDashToken(line, i, j)) {
                            sawSep = true;
                            i = j;
                            break;
                        }
                        left.add(parseIntAt(line, i, j));
                        i = j;
                    }
                    if (!sawSep) continue;

                    // right ints until parseEnd
                    while (true) {
                        i = nextNonSpace(line, i);
                        if (i >= parseEnd) break;
                        j = nextSpace(line, i);
                        right.add(parseIntAt(line, i, j));
                        i = j;
                    }

                    int m = Math.min(left.size, right.size);
                    for (int k = 0; k < m; k++) {
                        int a = left.a[k];
                        int b = right.a[k];
                        int ia = Math.abs(a);
                        int ib = Math.abs(b);
                        if (ia < 1 || ia > n || ib < 1 || ib > n) continue;

                        int w = ((a < 0) ^ (b < 0)) ? 1 : 0; // opposite sign -> 1
                        dsu.union(ia, ib, w);
                    }
                }
            }

            // -------- group by root --------
            IntList[] groups = new IntList[n + 1];
            for (int id = 1; id <= n; id++) {
                int r = dsu.find(id);
                IntList g = groups[r];
                if (g == null) groups[r] = g = new IntList(dsu.size[r]);
                g.add(id);
            }

            // -------- write outputs with Canonicalization (Scheme A) --------
            File outDir = new File(args[2]);
            try (BufferedWriter bwBlocks = new BufferedWriter(new FileWriter(new File(outDir, "blocks.fasta")), 1 << 20);
                 BufferedWriter bwSize = new BufferedWriter(new FileWriter(new File(outDir, "blocks_Size.fasta")), 1 << 20)) {

                bwSize.write("index   size   averageL   minL   maxL");
                bwSize.newLine();

                int blockIdx = 0;
                for (int r = 1; r <= n; r++) {
                    IntList g = groups[r];
                    if (g == null || g.size <= 1) continue;

                    // compute stats + decide flip so that smallest id is positive
                    int sum = 0;
                    int minL = Integer.MAX_VALUE;
                    int maxL = Integer.MIN_VALUE;

                    int minId = Integer.MAX_VALUE;     // smallest absolute id (here == smallest id)
                    int minIdSigned = 0;              // signed value of that minId under DSU parity

                    for (int k = 0; k < g.size; k++) {
                        int id = g.a[k];
                        if (id < minId) minId = id;

                        int L = lengths[id - 1];
                        sum += L;
                        if (L < minL) minL = L;
                        if (L > maxL) maxL = L;
                    }

                    // find signed value for minId
                    byte pMin = dsu.parityToRoot(minId);
                    minIdSigned = (pMin == 0) ? minId : -minId;

                    int flip = (minIdSigned < 0) ? -1 : 1;

                    // output members with flip applied
                    for (int k = 0; k < g.size; k++) {
                        int id = g.a[k];
                        byte p = dsu.parityToRoot(id);
                        int signed = (p == 0) ? id : -id;
                        signed *= flip;
                        bwBlocks.write(signed + ", ");
                    }
                    bwBlocks.newLine();
                    bwBlocks.newLine();

                    int avg = sum / g.size;
                    bwSize.write(blockIdx + "  " + g.size + "  " + avg + "  " + minL + "  " + maxL);
                    bwSize.newLine();
                    blockIdx++;
                }
            }

            long endTime = System.currentTimeMillis();
            System.out.println("end sdblock! runTime: " + (endTime - startTime));

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
