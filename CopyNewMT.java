import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLongArray;

public class CopyNewMT {

    // ========== 1) 互补碱基映射表 ==========
    private static final char[] COMP = new char[128];
    static {
        COMP['N'] = 'N';
        COMP['A'] = 'T';
        COMP['C'] = 'G';
        COMP['G'] = 'C';
        COMP['T'] = 'A';
        COMP['a'] = 't';
        COMP['c'] = 'g';
        COMP['g'] = 'c';
        COMP['t'] = 'a';
        COMP['n'] = 'n';
    }

    // ========== 2) 线程安全 Bitset（AtomicLongArray + CAS） ==========
    static final class ConcurrentBitSet {
        final AtomicLongArray words;
        final int size;

        ConcurrentBitSet(int size) {
            this.size = size;
            this.words = new AtomicLongArray((size + 63) >>> 6);
        }

        boolean get(int index) {
            long w = words.get(index >>> 6);
            return (w & (1L << (index & 63))) != 0L;
        }

        boolean set(int index) {
            int w = index >>> 6;
            long mask = 1L << (index & 63);
            while (true) {
                long old = words.get(w);
                long nw = old | mask;
                if (nw == old) return false;
                if (words.compareAndSet(w, old, nw)) return true;
            }
        }
    }

    // ========== 3) int[] 缓冲区复用（线程本地使用） ==========
    static final class IntBuffer {
        int[] a = new int[0];
        int[] b = new int[0];
        int[] t0 = new int[0];
        int[] t1 = new int[0];

        void ensureSegCapacity(int lenA, int lenB) {
            if (a.length < lenA) a = new int[lenA];
            if (b.length < lenB) b = new int[lenB];
        }

        void ensureTempCapacity(int len) {
            if (t0.length < len) t0 = new int[len];
            if (t1.length < len) t1 = new int[len];
        }
    }

    // ========== 4) 高速按空白分割 ==========
    static final class Splitter {
        static int splitByWhitespace(String line, String[] cols) {
            int n = line.length();
            int k = 0;
            int i = 0;
            while (i < n) {
                while (i < n && line.charAt(i) <= ' ') i++;
                if (i >= n) break;
                int j = i + 1;
                while (j < n && line.charAt(j) > ' ') j++;
                cols[k++] = line.substring(i, j);
                i = j;
            }
            return k;
        }
    }

    // ========== reverse：查表互补 ==========
    public static String reverse(String str) {
        int n = str.length();
        char[] temp = new char[n];
        for (int i = n - 1; i >= 0; i--) {
            char c = str.charAt(i);
            temp[n - 1 - i] = (c < 128 ? COMP[c] : 0);
        }
        return String.valueOf(temp);
    }

    // ========== CIGAR: 计算对齐长度 ==========
    private static int cigarAlignedLength(String cigar) {
        int len = 0, num = 0;
        for (int i = 0; i < cigar.length(); i++) {
            char c = cigar.charAt(i);
            if (c >= '0' && c <= '9') {
                num = num * 10 + (c - '0');
            } else {
                len += num;
                num = 0;
            }
        }
        return len;
    }

    // ========== CIGAR alignment ==========
    public static String[] cigarAlignment(String seg1, String seg2, String cigar) {
        int length = cigarAlignedLength(cigar);
        char[] str1 = new char[length];
        char[] str2 = new char[length];

        int index1 = 0, index2 = 0, index = 0;
        int num = 0;

        for (int i = 0; i < cigar.length(); i++) {
            char op = cigar.charAt(i);
            if (op >= '0' && op <= '9') {
                num = num * 10 + (op - '0');
                continue;
            }
            if (op == 'M') {
                while (num-- > 0) {
                    while (seg1.charAt(index1) == 'N') index1++;
                    str1[index] = seg1.charAt(index1++);
                    while (seg2.charAt(index2) == 'N') index2++;
                    str2[index] = seg2.charAt(index2++);
                    index++;
                }
            } else if (op == 'D') {
                while (num-- > 0) {
                    while (seg1.charAt(index1) == 'N') index1++;
                    str1[index] = seg1.charAt(index1++);
                    str2[index] = '-';
                    index++;
                }
            } else { // 'I'
                while (num-- > 0) {
                    str1[index] = '-';
                    while (seg2.charAt(index2) == 'N') index2++;
                    str2[index] = seg2.charAt(index2++);
                    index++;
                }
            }
            num = 0;
        }
        return new String[]{ String.valueOf(str1), String.valueOf(str2) };
    }

    public static String[] endsCigarAlignment(String seg1, String seg2, String cigar) {
        int length = cigarAlignedLength(cigar);
        char[] str1 = new char[length];
        char[] str2 = new char[length];

        int index1 = 0, index2 = 0, index = 0;
        int num = 0;

        for (int i = 0; i < cigar.length(); i++) {
            char op = cigar.charAt(i);
            if (op >= '0' && op <= '9') {
                num = num * 10 + (op - '0');
                continue;
            }
            if (op == 'M') {
                while (num-- > 0) {
                    str1[index] = seg1.charAt(index1++);
                    str2[index] = seg2.charAt(index2++);
                    index++;
                }
            } else if (op == 'D') {
                while (num-- > 0) {
                    str1[index] = seg1.charAt(index1++);
                    str2[index] = '-';
                    index++;
                }
            } else { // 'I'
                while (num-- > 0) {
                    str1[index] = '-';
                    str2[index] = seg2.charAt(index2++);
                    index++;
                }
            }
            num = 0;
        }
        return new String[]{ String.valueOf(str1), String.valueOf(str2) };
    }

    // ========== 拼接：避免链式 + 多次复制 ==========
    private static String concat3(String a, String b, String c) {
        int la = (a == null ? 0 : a.length());
        int lb = (b == null ? 0 : b.length());
        int lc = (c == null ? 0 : c.length());
        char[] out = new char[la + lb + lc];
        int p = 0;
        if (la > 0) { a.getChars(0, la, out, p); p += la; }
        if (lb > 0) { b.getChars(0, lb, out, p); p += lb; }
        if (lc > 0) { c.getChars(0, lc, out, p); }
        return String.valueOf(out);
    }

    // computeIndex：逻辑不变（线程本地 buf）
    public static void computeIndex(
            String[] alignment,
            int[] seg1To2,
            int[] seg2To1,
            String seg1,
            String seg2,
            String str1,
            String str2,
            IntBuffer buf
    ) {
        int len = alignment[0].length();
        buf.ensureTempCapacity(len);
        int[] temp0 = buf.t0;
        int[] temp1 = buf.t1;

        int index;

        if (str1.equalsIgnoreCase("+")) {
            index = 0;
            while (seg1.charAt(index) == 'N') index++;
            int p = 0;
            while (alignment[0].charAt(p) == '-' || alignment[0].charAt(p) == 'N') {
                temp0[p++] = index;
            }
            temp0[p] = index;
            for (int i = p + 1; i < len; i++) {
                if (alignment[0].charAt(i) != '-' && alignment[0].charAt(i) != 'N') {
                    index++;
                    while (seg1.charAt(index) == 'N') index++;
                }
                temp0[i] = index;
            }
        } else {
            index = 0;
            while (seg1.charAt(index) == 'N') index++;
            int p = 0;
            while (alignment[0].charAt(p) == '-' || alignment[0].charAt(p) == 'N') {
                temp0[p++] = seg1.length() - 1 - index;
            }
            temp0[p] = seg1.length() - 1 - index;
            for (int i = p + 1; i < len; i++) {
                if (alignment[0].charAt(i) != '-' && alignment[0].charAt(i) != 'N') {
                    index++;
                    while (seg1.charAt(index) == 'N') index++;
                }
                temp0[i] = seg1.length() - 1 - index;
            }
        }

        if (str2.equalsIgnoreCase("+")) {
            index = 0;
            while (seg2.charAt(index) == 'N') index++;
            int p = 0;
            while (alignment[1].charAt(p) == '-' || alignment[1].charAt(p) == 'N') {
                temp1[p++] = index;
            }
            temp1[p] = index;
            for (int i = p + 1; i < len; i++) {
                if (alignment[1].charAt(i) != '-' && alignment[1].charAt(i) != 'N') {
                    index++;
                    while (seg2.charAt(index) == 'N') index++;
                }
                temp1[i] = index;
            }
        } else {
            index = 0;
            while (seg2.charAt(index) == 'N') index++;
            int p = 0;
            while (alignment[1].charAt(p) == '-' || alignment[1].charAt(p) == 'N') {
                temp1[p++] = seg2.length() - 1 - index;
            }
            temp1[p] = seg2.length() - 1 - index;
            for (int i = p + 1; i < len; i++) {
                if (alignment[1].charAt(i) != '-' && alignment[1].charAt(i) != 'N') {
                    index++;
                    while (seg2.charAt(index) == 'N') index++;
                }
                temp1[i] = seg2.length() - 1 - index;
            }
        }

        for (int i = 0; i < len; i++) {
            if (alignment[0].charAt(i) != '-' && alignment[0].charAt(i) != 'N') {
                seg1To2[temp0[i]] = temp1[i];
            }
            if (alignment[1].charAt(i) != '-' && alignment[1].charAt(i) != 'N') {
                seg2To1[temp1[i]] = temp0[i];
            }
        }

        index = 0;
        int start = index;
        while (index < seg1To2.length && seg1To2[index] == -1) index++;
        for (int i = start; i < index; i++) seg1To2[i] = seg1To2[index];
        index++;
        start = index;
        while (index < seg1To2.length) {
            while (index < seg1To2.length && seg1To2[index] == -1) index++;
            for (int i = start; i < index; i++) seg1To2[i] = seg1To2[start - 1];
            index++;
            start = index;
        }

        index = 0;
        start = index;
        while (index < seg2To1.length && seg2To1[index] == -1) index++;
        for (int i = start; i < index; i++) seg2To1[i] = seg2To1[index];
        index++;
        start = index;
        while (index < seg2To1.length) {
            while (index < seg2To1.length && seg2To1[index] == -1) index++;
            for (int i = start; i < index; i++) seg2To1[i] = seg2To1[start - 1];
            index++;
            start = index;
        }
    }

    public static boolean copyEnds(
            ConcurrentBitSet line1,
            ConcurrentBitSet line2,
            int[] seg1to2,
            int[] seg2to1,
            int start1, int end1,
            int start2, int end2
    ) {
        boolean copy = false;

        for (int i = start1; i < end1; i++) {
            if (line1.get(i)) {
                int index1 = i - start1;
                int index2 = seg1to2[index1];
                if (line2.set(start2 + index2)) copy = true;
            }
        }

        for (int i = start2; i < end2; i++) {
            if (line2.get(i)) {
                int index2 = i - start2;
                int index1 = seg2to1[index2];
                if (line1.set(start1 + index1)) copy = true;
            }
        }

        return copy;
    }

    // ========== 输出：按块写 0/1 ==========
    private static void writeBitsetLine(BufferedOutputStream out, ConcurrentBitSet bs) throws IOException {
        int n = bs.size;
        byte[] buf = new byte[1 << 16]; // 64KB
        int bp = 0;

        for (int i = 0; i < n; i++) {
            buf[bp++] = (byte) (bs.get(i) ? '1' : '0');
            if (bp == buf.length) {
                out.write(buf, 0, bp);
                bp = 0;
            }
        }
        if (bp > 0) out.write(buf, 0, bp);
        out.write('\n');
    }

    // 单条记录处理：从 cols 里取值，逻辑与原来 while(in2.readLine()) 中一致
    private static boolean processOneRecord(
            String[] cols,
            int cnt,
            char[][] chroms,
            ConcurrentBitSet[] linesBS,
            IntBuffer buf
    ) {
        if (cnt <= 21) return false;

        String cigar  = cols[1];
        String cigarS = cols[20];
        String cigarE = cols[21];

        int chr = Integer.parseInt(cols[4]);
        int start = Integer.parseInt(cols[5]);
        int end = Integer.parseInt(cols[6]);

        int otherChr = Integer.parseInt(cols[7]);
        int otherStart = Integer.parseInt(cols[8]);
        int otherEnd = Integer.parseInt(cols[9]);

        int extendS1 = Integer.parseInt(cols[16]);
        int extendE1 = Integer.parseInt(cols[17]);
        int extendS2 = Integer.parseInt(cols[18]);
        int extendE2 = Integer.parseInt(cols[19]);

        String strand1 = cols[14];
        String strand2 = cols[15];

        String segS1 = "", segS2 = "", segE1 = "", segE2 = "";
        String[] alignment;
        String[] alignS;
        String[] alignE;

        String seg1, seg2;

        // seg1（逻辑不变）
        if (strand1.equalsIgnoreCase("+")) {
            seg1 = new String(chroms[chr], start, end - start);
            if (!cigarS.equalsIgnoreCase("-1")) {
                segS1 = new String(chroms[chr], extendS1, start - extendS1);
                start = extendS1;
            }
            if (!cigarE.equalsIgnoreCase("-1")) {
                segE1 = new String(chroms[chr], end, extendE1 - end);
                end = extendE1;
            }
        } else {
            seg1 = reverse(new String(chroms[chr], start, end - start));
            if (!cigarS.equalsIgnoreCase("-1")) {
                segS1 = reverse(new String(chroms[chr], end, extendE1 - end));
                end = extendE1;
            }
            if (!cigarE.equalsIgnoreCase("-1")) {
                segE1 = reverse(new String(chroms[chr], extendS1, start - extendS1));
                start = extendS1;
            }
        }

        // seg2（逻辑不变）
        if (strand2.equalsIgnoreCase("+")) {
            seg2 = new String(chroms[otherChr], otherStart, otherEnd - otherStart);
            if (!cigarS.equalsIgnoreCase("-1")) {
                segS2 = new String(chroms[otherChr], extendS2, otherStart - extendS2);
                otherStart = extendS2;
            }
            if (!cigarE.equalsIgnoreCase("-1")) {
                segE2 = new String(chroms[otherChr], otherEnd, extendE2 - otherEnd);
                otherEnd = extendE2;
            }
        } else {
            seg2 = reverse(new String(chroms[otherChr], otherStart, otherEnd - otherStart));
            if (!cigarS.equalsIgnoreCase("-1")) {
                segS2 = reverse(new String(chroms[otherChr], otherEnd, extendE2 - otherEnd));
                otherEnd = extendE2;
            }
            if (!cigarE.equalsIgnoreCase("-1")) {
                segE2 = reverse(new String(chroms[otherChr], extendS2, otherStart - extendS2));
                otherStart = extendS2;
            }
        }

        // 主对齐
        alignment = cigarAlignment(seg1, seg2, cigar);

        // 头尾对齐
        String pre0 = null, pre1 = null, suf0 = null, suf1 = null;

        if (!cigarS.equalsIgnoreCase("-1")) {
            alignS = endsCigarAlignment(segS1, segS2, cigarS);
            pre0 = alignS[0];
            pre1 = alignS[1];
            seg1 = segS1 + seg1;
            seg2 = segS2 + seg2;
        }

        if (!cigarE.equalsIgnoreCase("-1")) {
            alignE = endsCigarAlignment(segE1, segE2, cigarE);
            suf0 = alignE[0];
            suf1 = alignE[1];
            seg1 = seg1 + segE1;
            seg2 = seg2 + segE2;
        }

        if (pre0 != null || suf0 != null) {
            alignment[0] = concat3(pre0, alignment[0], suf0);
            alignment[1] = concat3(pre1, alignment[1], suf1);
        }

        // seg1to2/seg2to1：线程本地复用
        buf.ensureSegCapacity(seg1.length(), seg2.length());
        int[] seg1to2 = buf.a;
        int[] seg2to1 = buf.b;

        Arrays.fill(seg1to2, 0, seg1.length(), -1);
        Arrays.fill(seg2to1, 0, seg2.length(), -1);

        computeIndex(alignment, seg1to2, seg2to1, seg1, seg2, strand1, strand2, buf);

        return copyEnds(linesBS[chr], linesBS[otherChr],
                seg1to2, seg2to1, start, end, otherStart, otherEnd);
    }

    public static void main(String args[]) {
        try {
            long startTime = System.currentTimeMillis();

            // === 读 args[0] fasta：顺序不变（header -> seq） ===
            BufferedReader inChr = new BufferedReader(new FileReader(args[0]), 1 << 20);
            char[][] chromsTmp = new char[1024][];
            ConcurrentBitSet[] linesTmp = new ConcurrentBitSet[1024];
            int chromCount = 0;

            String line;
            while ((line = inChr.readLine()) != null) {
                String seqLine = inChr.readLine();
                if (seqLine == null) break;

                if (chromCount == chromsTmp.length) {
                    chromsTmp = Arrays.copyOf(chromsTmp, chromsTmp.length << 1);
                    linesTmp = Arrays.copyOf(linesTmp, linesTmp.length << 1);
                }
                char[] seq = seqLine.toCharArray();
                chromsTmp[chromCount] = seq;
                linesTmp[chromCount] = new ConcurrentBitSet(seq.length);
                chromCount++;
            }
            inChr.close();

            char[][] chroms = Arrays.copyOf(chromsTmp, chromCount);
            ConcurrentBitSet[] linesBS = Arrays.copyOf(linesTmp, chromCount);

            // === 读 args[1] endpoints：顺序不变 ===
            BufferedReader in1 = new BufferedReader(new FileReader(args[1]), 1 << 20);
            String[] cols1 = new String[8];
            while ((line = in1.readLine()) != null) {
                int cnt = Splitter.splitByWhitespace(line, cols1);
                if (cnt < 3) continue;

                int chr = Integer.parseInt(cols1[0]);
                int start = Integer.parseInt(cols1[1]);
                int end = Integer.parseInt(cols1[2]) - 1;

                char[] seq = chroms[chr];
                while (start < end && seq[start] == 'N') start++;
                while (start < end && seq[end] == 'N') end--;

                linesBS[chr].set(start);
                linesBS[chr].set(end);
            }
            in1.close();

            // === 关键优化：args[2] 只读一次，避免每轮重扫文件 ===
            final List<String> records = new ArrayList<>(1 << 20);
            try (BufferedReader in2 = new BufferedReader(new FileReader(args[2]), 1 << 20)) {
                in2.readLine(); // skip first line
                while ((line = in2.readLine()) != null) {
                    records.add(line);
                }
            }

            System.out.println("Thread: "+Runtime.getRuntime().availableProcessors());
            final int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
            final ExecutorService pool = Executors.newFixedThreadPool(threads);

            // 线程本地：避免共享对象导致错乱
            final ThreadLocal<IntBuffer> tlBuf = ThreadLocal.withInitial(IntBuffer::new);
            final ThreadLocal<String[]> tlCols = ThreadLocal.withInitial(() -> new String[32]);

            boolean copy;
            do {
                AtomicBoolean changed = new AtomicBoolean(false);

                final int n = records.size();
                final int chunk = (n + threads - 1) / threads;
                List<Future<?>> futures = new ArrayList<>(threads);

                for (int t = 0; t < threads; t++) {
                    final int from = t * chunk;
                    final int to = Math.min(n, (t + 1) * chunk);
                    if (from >= to) break;

                    futures.add(pool.submit(() -> {
                        IntBuffer buf = tlBuf.get();
                        String[] cols = tlCols.get();

                        for (int i = from; i < to; i++) {
                            String recLine = records.get(i);
                            int cnt = Splitter.splitByWhitespace(recLine, cols);
                            if (cnt <= 21) continue;

                            boolean tempCopy = processOneRecord(cols, cnt, chroms, linesBS, buf);
                            if (tempCopy) changed.set(true);
                        }
                    }));
                }

                for (Future<?> f : futures) {
                    f.get(); // 等待本轮全部处理结束
                }

                copy = changed.get();
            } while (copy);

            pool.shutdown();

            // === 输出顺序/格式不变 ===
            try (BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(args[3]), 1 << 20)) {
                for (int n = 0; n < linesBS.length; n++) {
                    out.write('>');
                    byte[] numBytes = Integer.toString(n).getBytes();
                    out.write(numBytes);
                    out.write('\n');

                    writeBitsetLine(out, linesBS[n]);
                }
                out.flush();
            }

            long endTime = System.currentTimeMillis();
            System.out.println("Running Time:" + (endTime - startTime));

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
