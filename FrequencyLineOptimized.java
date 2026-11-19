import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.RandomAccessFile;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;

public class FrequencyLineOptimized {

    // 稀疏索引间隔（每 BLOCK 条记录一个索引点）
    private static final int BLOCK = 4096;

    // 二进制记录格式：long key (8) + byte freq (1)
    private static final int RECORD_SIZE = 9;

    // 分段 mmap，避免 2GB/Integer 溢出
    private static final long SEGMENT_SIZE = 1L << 30; // 1GB

    private static int baseTo2Bits(char c) {
        switch (c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default:  return -1; // N or others
        }
    }

    // 从一段字符（0..len-1）编码成 long；遇到非法字符返回 -1
    private static long encodeKmerFromLine(final String line, final int len) {
        long v = 0L;
        for (int i = 0; i < len; i++) {
            int b = baseTo2Bits(line.charAt(i));
            if (b < 0) return -1L;
            v = (v << 2) | (long) b;
        }
        return v;
    }

    // 反向互补：对 2bit 编码的 k-mer 做 reverse + complement
    private static long revComp(long v, int k) {
        long rc = 0L;
        for (int i = 0; i < k; i++) {
            long b = v & 3L;
            rc = (rc << 2) | (b ^ 3L); // complement
            v >>>= 2;
        }
        return rc;
    }

    // 轻量动态数组（避免装箱）
    private static final class LongGrow {
        long[] a;
        int size;
        LongGrow(int cap) { a = new long[Math.max(1, cap)]; }
        void add(long v) {
            if (size == a.length) {
                long[] b = new long[a.length << 1];
                System.arraycopy(a, 0, b, 0, a.length);
                a = b;
            }
            a[size++] = v;
        }
        long[] toArray() {
            long[] r = new long[size];
            System.arraycopy(a, 0, r, 0, size);
            return r;
        }
        int size() { return size; }
    }

    // 解析 "KMER<space>FREQ"：返回 freq 起始下标；若格式不合法返回 -1
    private static int findFreqStart(final String line) {
        int n = line.length();
        int i = 0;
        while (i < n && line.charAt(i) != ' ' && line.charAt(i) != '\t') i++;
        if (i == n) return -1;
        while (i < n && (line.charAt(i) == ' ' || line.charAt(i) == '\t')) i++;
        if (i == n) return -1;
        return i;
    }

    /**
     * 将字典序 k-mer 文本构建为二进制索引 + 稀疏索引。
     * 输出文件名不依赖扩展名：在输入路径后追加固定后缀。
     */
    private static void buildIndexIfNeeded(String kmerTxtPath, String binPath, String sidxPath, int k) throws IOException {

        File txt = new File(kmerTxtPath);
        File bin = new File(binPath);
        File sidx = new File(sidxPath);

        // 若索引存在且比文本新，直接复用
        if (bin.exists() && sidx.exists()
                && bin.lastModified() >= txt.lastModified()
                && sidx.lastModified() >= txt.lastModified()) {
            System.out.println("Index up-to-date. Reusing: " + binPath + " and " + sidxPath);
            return;
        }

        System.out.println("Building index from: " + kmerTxtPath);
        System.out.println("Output bin : " + binPath);
        System.out.println("Output sidx: " + sidxPath);

        LongGrow sparseKeys = new LongGrow(1 << 20);
        LongGrow sparseOffsets = new LongGrow(1 << 20);

        long recordIndex = 0;

        // 大缓冲：降低 IO 开销
        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(new FileInputStream(txt), StandardCharsets.US_ASCII), 1 << 20);
             BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(bin), 1 << 20);
             DataOutputStream dos = new DataOutputStream(bos)) {

            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                // k-mer 应在行首；若行太短跳过
                if (line.length() < k + 2) continue;

                // 找 freq 起点
                int freqStart = findFreqStart(line);
                if (freqStart < 0) continue;

                // 编码 k-mer（不创建 substring）
                long key = encodeKmerFromLine(line, k);
                if (key < 0) continue; // 含 N 或非法字符则跳过（一般你的数据不会有）

                int freq;
                try {
                    freq = Integer.parseInt(line.substring(freqStart));
                } catch (NumberFormatException nfe) {
                    continue;
                }
                if (freq > 6) freq = 6;
                if (freq < 0) freq = 0;

                if ((recordIndex % BLOCK) == 0) {
                    sparseKeys.add(key);
                    sparseOffsets.add(recordIndex);
                }

                dos.writeLong(key);
                dos.writeByte((byte) freq);

                recordIndex++;
                if ((recordIndex % 5_000_000) == 0) {
                    System.out.println("Indexed records: " + recordIndex);
                }
            }

            dos.flush();
        }

        // 写稀疏索引文件：int count + (key, offset) * count
        try (DataOutputStream out = new DataOutputStream(
                new BufferedOutputStream(new FileOutputStream(sidx), 1 << 20))) {
            int n = sparseKeys.size();
            out.writeInt(n);
            long[] keys = sparseKeys.toArray();
            long[] offs = sparseOffsets.toArray();
            for (int i = 0; i < n; i++) {
                out.writeLong(keys[i]);
                out.writeLong(offs[i]);
            }
            out.flush();
        }

        System.out.println("Index build done. Total records: " + recordIndex + ", sparse points: " + sparseKeys.size());
    }

    private static final class KmerStore {
        final long totalRecords;
        final MappedByteBuffer[] segments;
        final long[] sKeys;
        final long[] sOffsets;

        final long segmentSizeAligned;

        KmerStore(String binFile, String idxFile) throws IOException {

            try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(idxFile), 1 << 20))) {
                int n = in.readInt();
                sKeys = new long[n];
                sOffsets = new long[n];
                for (int i = 0; i < n; i++) {
                    sKeys[i] = in.readLong();
                    sOffsets[i] = in.readLong();
                }
            }
        
            File f = new File(binFile);
            long size = f.length();
            this.totalRecords = size / RECORD_SIZE;
        
            long aligned = (SEGMENT_SIZE / RECORD_SIZE) * RECORD_SIZE; // 对齐到 9 的倍数
            if (aligned < RECORD_SIZE) aligned = RECORD_SIZE;          // 极端情况保护
            this.segmentSizeAligned = aligned;
        
            int numSegs = (int) ((size + segmentSizeAligned - 1) / segmentSizeAligned);
            segments = new MappedByteBuffer[numSegs];
        
            try (FileChannel ch = new RandomAccessFile(f, "r").getChannel()) {
                long off = 0;
                for (int i = 0; i < numSegs; i++) {
                    long segSize = Math.min(segmentSizeAligned, size - off);
                    MappedByteBuffer mbb = ch.map(FileChannel.MapMode.READ_ONLY, off, segSize);
                    mbb.order(ByteOrder.BIG_ENDIAN);
                    segments[i] = mbb;
                    off += segSize;
                }
            }
        }
        

        private long readKeyAt(long recordIndex) {
            long bytePos = recordIndex * RECORD_SIZE;
            int seg = (int) (bytePos / segmentSizeAligned);
            int off = (int) (bytePos % segmentSizeAligned);
            return segments[seg].getLong(off);
        }
        
        private byte readFreqAt(long recordIndex) {
            long bytePos = recordIndex * RECORD_SIZE + 8;
            int seg = (int) (bytePos / segmentSizeAligned);
            int off = (int) (bytePos % segmentSizeAligned);
            return segments[seg].get(off);
        }
        

        // 找到最后一个 sKeys[idx] <= key 的 idx（若全大于则返回 -1）
        private int findBlock(long key) {
            int lo = 0, hi = sKeys.length - 1;
            int ans = -1;
            while (lo <= hi) {
                int mid = (lo + hi) >>> 1;
                long v = sKeys[mid];
                if (v <= key) {
                    ans = mid;
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }
            return ans;
        }

        // 返回 0 表示找不到（你的语义是：找不到输出 '1'）
        byte lookup(long key) {
            if (totalRecords == 0) return 0;

            int b = findBlock(key);
            long start = (b < 0) ? 0 : sOffsets[b];
            long end = Math.min(start + BLOCK, totalRecords);

            long lo = start, hi = end - 1;
            while (lo <= hi) {
                long mid = (lo + hi) >>> 1;
                long mk = readKeyAt(mid);
                if (mk == key) return readFreqAt(mid);
                if (key < mk) hi = mid - 1;
                else lo = mid + 1;
            }
            return 0;
        }
    }

    public static void main(String[] args) {

        // 保持 4 个参数，与原程序一致
        // args[0]=kmer_text (25G, lex sorted)
        // args[1]=fasta (one line per chr)
        // args[2]=output
        // args[3]=k
        if (args.length < 4) {
            System.err.println("Usage: java FrequencyLineOptimized <kmer_text> <genome_fasta> <output> <k>");
            System.exit(1);
        }

        final String kmerTxtPath = args[0];
        final String fastaPath = args[1];
        final String outPath = args[2];
        final int k = Integer.parseInt(args[3]);

        if (k <= 1 || k > 31) {
            System.err.println("k must be in range 2..31");
            System.exit(1);
        }

        // 输出索引文件名：不依赖扩展名，只是追加固定后缀
        final String binPath = kmerTxtPath + ".kmerbin";
        final String sidxPath = kmerTxtPath + ".kmersidx";

        final long mask = (1L << (2 * k)) - 1;

        long startTime = System.currentTimeMillis();

        try {
            // 1) 构建/复用索引
            buildIndexIfNeeded(kmerTxtPath, binPath, sidxPath, k);

            // 2) 打开 mmap 索引
            KmerStore store = new KmerStore(binPath, sidxPath);

            // 3) 扫描 fasta 输出
            try (BufferedReader br = new BufferedReader(
                        new InputStreamReader(new FileInputStream(fastaPath), StandardCharsets.US_ASCII), 1 << 20);
                 BufferedWriter bw = new BufferedWriter(
                        new OutputStreamWriter(new FileOutputStream(outPath), StandardCharsets.US_ASCII), 1 << 20)) {

                String header;
                int chrCount = 0;

                while ((header = br.readLine()) != null) {

                    bw.write(header);
                    bw.newLine();

                    String seq = br.readLine(); // 你已确认序列单行
                    if (seq == null) break;
                    seq = seq.trim();

                    chrCount++;
                    int len = seq.length();
                    System.out.println("Handling chromosome " + chrCount + " length=" + len);

                    int outLen = Math.max(0, len - k + 1);
                    char[] out = new char[outLen];

                    // 用环形缓存维护窗口内是否含 N：code=-1 表示 N/非法
                    byte[] ring = new byte[k];
                    int ringPos = 0;
                    int nInWindow = 0;

                    long window = 0L;
                    int filled = 0;

                    for (int i = 0; i < len; i++) {
                        int code = baseTo2Bits(seq.charAt(i));

                        // 移出旧元素
                        if (filled >= k) {
                            byte old = ring[ringPos];
                            if (old < 0) nInWindow--;
                        } else {
                            filled++;
                        }

                        // 写入新元素
                        ring[ringPos] = (byte) code;
                        if (code < 0) nInWindow++;
                        ringPos++;
                        if (ringPos == k) ringPos = 0;

                        // rolling window
                        window = ((window << 2) & mask) | (code < 0 ? 0L : (long) code);

                        if (i >= k - 1) {
                            int pos = i - k + 1;

                            if (nInWindow == 0) {
                                long rc = revComp(window, k);

                                byte f = store.lookup(window);
                                if (f == 0) f = store.lookup(rc);

                                // 原逻辑：找不到写 "1"，找到写 freq(1..6)
                                out[pos] = (f == 0) ? '1' : (char) ('0' + f);

                            } else {
                                // 复刻你原来的 N 规则：
                                // 窗口含 N 且窗口首位 != 'N' -> '1'
                                // 否则（窗口首位为 N）-> '0'
                                char first = seq.charAt(pos);
                                out[pos] = (first != 'N') ? '1' : '0';
                            }
                        }
                    }

                    bw.write(out, 0, out.length);
                    bw.newLine();
                }

                bw.flush();
            }

            long endTime = System.currentTimeMillis();
            System.out.println("Done. elapsed(ms)=" + (endTime - startTime));

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(2);
        }
    }
}
