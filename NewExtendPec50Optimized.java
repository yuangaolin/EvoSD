import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class NewExtendPec50Optimized {

    // fixed extension window
    static final int EXTEND = 50;

    // ---------------------------
    // Fast line field parsing
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

    // Returns [start,end) of field k (0-based)
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
            char c = s.charAt(i);
            val = val * 10 + (c - '0');
        }
        return val * sign;
    }

    static char parseCharField(String s, int k) {
        int[] r = fieldRange(s, k);
        return (r[0] < r[1]) ? s.charAt(r[0]) : '\0';
    }

    static String parseStringField(String s, int k) {
        int[] r = fieldRange(s, k);
        return s.substring(r[0], r[1]);
    }

    // ---------------------------
    // Genome loading
    // ---------------------------
    static char[] toUpperCharArray(String s) {
        char[] a = s.toCharArray();
        for (int i = 0; i < a.length; i++) {
            char c = a[i];
            if (c >= 'a' && c <= 'z') a[i] = (char) (c - 32);
        }
        return a;
    }

    static ArrayList<char[]> readFastaSeqsRaw(String path) throws IOException {
        ArrayList<char[]> seqs = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path), 1 << 20)) {
            String header;
            while ((header = br.readLine()) != null) {
                String seq = br.readLine();
                if (seq == null) break;
                seqs.add(seq.toCharArray());
            }
        }
        return seqs;
    }

    static ArrayList<char[]> readFastaSeqsUpper(String path) throws IOException {
        ArrayList<char[]> seqs = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path), 1 << 20)) {
            String header;
            while ((header = br.readLine()) != null) {
                String seq = br.readLine();
                if (seq == null) break;
                seqs.add(toUpperCharArray(seq));
            }
        }
        return seqs;
    }

    // ---------------------------
    // Reverse-complement helpers
    // ---------------------------
    static char compUpper(char b) {
        switch (b) {
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            case 'N': return 'N';
            default:  return b;
        }
    }
    
    static char compPreserveCase(char b) {
        switch (b) {
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            case 'N': return 'N';
            case 'a': return 't';
            case 'c': return 'g';
            case 'g': return 'c';
            case 't': return 'a';
            case 'n': return 'n';
            default:  return b;
        }
    }
    

    static String reverseComplementRaw(char[] chromRaw, int start, int end) {
        int len = end - start;
        char[] out = new char[len];
        for (int i = 0; i < len; i++) {
            out[i] = compPreserveCase(chromRaw[end - 1 - i]);
        }
        return String.valueOf(out);
    }

    static String orientedSegmentRaw(char[] chromRaw, int start, int end, char strand) {
        if (strand == '+') {
            return new String(chromRaw, start, end - start);
        }
        return reverseComplementRaw(chromRaw, start, end);
    }

    static String orientedSegmentUpper(char[] chromUpper, int start, int end, char strand) {
        int len = end - start;
        char[] out = new char[len];
        if (strand == '+') {
            System.arraycopy(chromUpper, start, out, 0, len);
        } else {
            for (int i = 0; i < len; i++) {
                out[i] = compUpper(chromUpper[end - 1 - i]);
            }
        }
        return String.valueOf(out);
    }

    /**
     * Fill a 50bp window from genomic [winStart, winEnd) with orientation defined by strand.
     * Returns false if any base is 'N' or out of bounds.
     */
    static boolean fillWindow50Upper(char[] chromUpper, int winStart, int winEnd, char strand, char[] out50) {
        if (winEnd - winStart != EXTEND) return false;
        if (winStart < 0 || winEnd > chromUpper.length) return false;

        if (strand == '+') {
            for (int i = 0; i < EXTEND; i++) {
                char b = chromUpper[winStart + i];
                out50[i] = b;
                if (b == 'N') return false;
            }
        } else {
            // reverse-complement of the genomic interval
            for (int i = 0; i < EXTEND; i++) {
                char b = compUpper(chromUpper[winEnd - 1 - i]);
                out50[i] = b;
                if (b == 'N') return false;
            }
        }
        return true;
    }

    // ---------------------------
    // Original helper methods (mostly unchanged)
    // ---------------------------
    public static int firstCR(String str) {
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == 'N') return i;
        }
        return str.length();
    }

    public static int lastCR(String str) {
        for (int i = str.length() - 1; i >= 0; i--) {
            if (str.charAt(i) == 'N') return str.length() - 1 - i;
        }
        return str.length();
    }

    public static int startIndex(int length, String str) {
        int num = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) != 'N') num++;
            if (num == length + 1) return i;
        }
        return str.length();
    }

    public static int endIndex(int length, String str) {
        int num = 0;
        for (int i = str.length() - 1; i >= 0; i--) {
            if (str.charAt(i) != 'N') num++;
            if (num == length + 1) return str.length() - 1 - i;
        }
        return str.length();
    }

    public static int[] cutStart(int length, String seg, String strand, int start, int end) {
        int[] temp = new int[2];
        if (strand.equalsIgnoreCase("+")) {
            temp[0] = start + startIndex(length, seg);
            temp[1] = end;
        } else {
            temp[0] = start;
            temp[1] = end - startIndex(length, seg);
        }
        return temp;
    }

    public static int[] cutEnd(int length, String seg, String strand, int start, int end) {
        int[] temp = new int[2];
        if (strand.equalsIgnoreCase("+")) {
            temp[0] = start;
            temp[1] = end - endIndex(length, seg);
        } else {
            temp[0] = start + endIndex(length, seg);
            temp[1] = end;
        }
        return temp;
    }

    public static int[] cutCigarStart(int index, int length, String[] alignment,
                                      String seg1, String strand1, int start1, int end1,
                                      String seg2, String strand2, int start2, int end2) {
        int num = 0;
        if (index == 0) {
            while (length > 0) {
                if (alignment[0].charAt(num) != '-') length--;
                num++;
            }
        } else {
            while (length > 0) {
                if (alignment[1].charAt(num) != '-') length--;
                num++;
            }
        }

        while (num < alignment[0].length() && alignment[0].charAt(num) != alignment[1].charAt(num)) num++;

        int num1 = 0, num2 = 0;
        for (int i = 0; i < num; i++) {
            if (alignment[0].charAt(i) != '-') num1++;
            if (alignment[1].charAt(i) != '-') num2++;
        }

        alignment[0] = alignment[0].substring(num);
        alignment[1] = alignment[1].substring(num);

        int[] temp1 = cutStart(num1, seg1, strand1, start1, end1);
        int[] temp2 = cutStart(num2, seg2, strand2, start2, end2);
        int[] newIndex = new int[4];
        newIndex[0] = temp1[0]; newIndex[1] = temp1[1]; newIndex[2] = temp2[0]; newIndex[3] = temp2[1];
        return newIndex;
    }

    public static int[] cutCigarEnd(int index, int length, String[] alignment,
                                    String seg1, String strand1, int start1, int end1,
                                    String seg2, String strand2, int start2, int end2) {
        int num = alignment[0].length() - 1;
        if (index == 0) {
            while (length > 0) {
                if (alignment[0].charAt(num) != '-') length--;
                num--;
            }
        } else {
            while (length > 0) {
                if (alignment[1].charAt(num) != '-') length--;
                num--;
            }
        }

        while (num >= 0 && alignment[0].charAt(num) != alignment[1].charAt(num)) num--;

        int num1 = 0, num2 = 0;
        for (int i = num + 1; i < alignment[0].length(); i++) {
            if (alignment[0].charAt(i) != '-') num1++;
            if (alignment[1].charAt(i) != '-') num2++;
        }

        alignment[0] = alignment[0].substring(0, num + 1);
        alignment[1] = alignment[1].substring(0, num + 1);

        int[] temp1 = cutEnd(num1, seg1, strand1, start1, end1);
        int[] temp2 = cutEnd(num2, seg2, strand2, start2, end2);
        int[] newIndex = new int[4];
        newIndex[0] = temp1[0]; newIndex[1] = temp1[1]; newIndex[2] = temp2[0]; newIndex[3] = temp2[1];
        return newIndex;
    }

    // ---------------------------
    // CIGAR alignment (unchanged)
    // ---------------------------
    public static String[] cigarAlignment(String seg1, String seg2, String cigar) {
        int start = 0, length = 0;
        int index1 = 0, index2 = 0, index = 0, num;

        for (int i = 0; i < cigar.length(); i++) {
            char c = cigar.charAt(i);
            if (c == 'M' || c == 'I' || c == 'D') {
                length += Integer.parseInt(cigar.substring(start, i));
                start = i + 1;
            }
        }
        char[] str1 = new char[length];
        char[] str2 = new char[length];
        start = 0;

        for (int i = 0; i < cigar.length(); i++) {
            char c = cigar.charAt(i);
            if (c == 'M') {
                num = Integer.parseInt(cigar.substring(start, i));
                start = i + 1;
                while (num-- > 0) {
                    while (seg1.charAt(index1) == 'N') index1++;
                    str1[index] = seg1.charAt(index1++);
                    while (seg2.charAt(index2) == 'N') index2++;
                    str2[index] = seg2.charAt(index2++);
                    index++;
                }
            } else if (c == 'D') {
                num = Integer.parseInt(cigar.substring(start, i));
                start = i + 1;
                while (num-- > 0) {
                    while (seg1.charAt(index1) == 'N') index1++;
                    str1[index] = seg1.charAt(index1++);
                    str2[index] = '-';
                    index++;
                }
            } else if (c == 'I') {
                num = Integer.parseInt(cigar.substring(start, i));
                start = i + 1;
                while (num-- > 0) {
                    str1[index] = '-';
                    while (seg2.charAt(index2) == 'N') index2++;
                    str2[index] = seg2.charAt(index2++);
                    index++;
                }
            }
        }
        return new String[]{String.valueOf(str1), String.valueOf(str2)};
    }

    public static String[] computeCigar(String[] alignment) {
        String[] temp = new String[5];
        temp[0] = "";
        int match = 0, mismatch = 0, sgap = 0, ngap = 0, num;
        int i = 0;
        while (i < alignment[0].length()) {
            if (alignment[0].charAt(i) != '-' && alignment[1].charAt(i) != '-') {
                if (alignment[0].charAt(i) == alignment[1].charAt(i)) match++;
                else mismatch++;
                num = 1; i++;
                while (i < alignment[0].length() && alignment[0].charAt(i) != '-' && alignment[1].charAt(i) != '-') {
                    if (alignment[0].charAt(i) == alignment[1].charAt(i)) match++;
                    else mismatch++;
                    num++; i++;
                }
                temp[0] += (num + "M");
            } else if (alignment[0].charAt(i) == '-' && alignment[1].charAt(i) != '-') {
                sgap++; ngap++; num = 1; i++;
                while (i < alignment[0].length() && alignment[0].charAt(i) == '-' && alignment[1].charAt(i) != '-') {
                    i++; sgap++; num++;
                }
                temp[0] += (num + "I");
            } else if (alignment[0].charAt(i) != '-' && alignment[1].charAt(i) == '-') {
                sgap++; ngap++; num = 1; i++;
                while (i < alignment[0].length() && alignment[0].charAt(i) != '-' && alignment[1].charAt(i) == '-') {
                    i++; sgap++; num++;
                }
                temp[0] += (num + "D");
            }
        }
        temp[1] = String.valueOf(match);
        temp[2] = String.valueOf(mismatch);
        temp[3] = String.valueOf(ngap);
        temp[4] = String.valueOf(sgap);
        return temp;
    }

    public static String[] updateAlign(int start, ArrayList<String> strA, ArrayList<String> strB) {
        StringBuilder temp1 = new StringBuilder();
        StringBuilder temp2 = new StringBuilder();
        if (start == 0) {
            for (int i = strA.size() - 1; i >= 0; i--) {
                temp1.append(strA.get(i));
                temp2.append(strB.get(i));
            }
        } else {
            for (int i = 0; i < strA.size(); i++) {
                temp1.append(strA.get(i));
                temp2.append(strB.get(i));
            }
        }
        return computeCigar(new String[]{temp1.toString(), temp2.toString()});
    }

    public static int Nlength(String str) {
        int dis = 0;
        for (int i = 0; i < str.length(); i++) {
            char c = str.charAt(i);
            if (c == 'N' || c == 'n') dis++;
        }
        return dis;
    }

    // ---------------------------
    // Alignment (optimized)
    // ---------------------------

    /**
     * Fast gate: compute match count and alignment length (steps) for two 50bp strings (as char[50]),
     * without producing alignment strings.
     */
    static int[] alignGate50(char[] a50, char[] b50, short[][] score, byte[][] ptr) {
        // ptr: 0 diag, 1 up, -1 left
        score[0][0] = 0; ptr[0][0] = 2;
        for (int j = 1; j <= EXTEND; j++) { score[0][j] = (short) (score[0][j - 1] - 1); ptr[0][j] = -1; }
        for (int i = 1; i <= EXTEND; i++) { score[i][0] = (short) (score[i - 1][0] - 1); ptr[i][0] = 1; }

        for (int i = 1; i <= EXTEND; i++) {
            char ca = a50[i - 1];
            for (int j = 1; j <= EXTEND; j++) {
                char cb = b50[j - 1];
                int best = score[i - 1][j - 1] + ((ca == cb) ? 1 : -1);
                byte bp = 0;
                int up = score[i - 1][j] - 1;
                if (up > best) { best = up; bp = 1; }
                int left = score[i][j - 1] - 1;
                if (left > best) { best = left; bp = -1; }
                score[i][j] = (short) best;
                ptr[i][j] = bp;
            }
        }

        int i = EXTEND, j = EXTEND;
        int match = 0;
        int steps = 0;
        while (!(i == 0 && j == 0)) {
            byte p = ptr[i][j];
            steps++;
            if (p == 0) {
                if (a50[i - 1] == b50[j - 1]) match++;
                i--; j--;
            } else if (p == 1) {
                i--;
            } else {
                j--;
            }
        }
        return new int[]{match, steps};
    }

    /**
     * Full alignment for two 50bp strings: returns {aln1, aln2, score, matchNum}.
     * Uses StringBuilder (no O(L^2) concatenation).
     */
    static String[] alignFull50(char[] a50, char[] b50, short[][] score, byte[][] ptr) {
        // reuse DP already computed? For simplicity, recompute deterministically.
        score[0][0] = 0; ptr[0][0] = 2;
        for (int j = 1; j <= EXTEND; j++) { score[0][j] = (short) (score[0][j - 1] - 1); ptr[0][j] = -1; }
        for (int i = 1; i <= EXTEND; i++) { score[i][0] = (short) (score[i - 1][0] - 1); ptr[i][0] = 1; }

        for (int i = 1; i <= EXTEND; i++) {
            char ca = a50[i - 1];
            for (int j = 1; j <= EXTEND; j++) {
                char cb = b50[j - 1];
                int best = score[i - 1][j - 1] + ((ca == cb) ? 1 : -1);
                byte bp = 0;
                int up = score[i - 1][j] - 1;
                if (up > best) { best = up; bp = 1; }
                int left = score[i][j - 1] - 1;
                if (left > best) { best = left; bp = -1; }
                score[i][j] = (short) best;
                ptr[i][j] = bp;
            }
        }

        StringBuilder t1 = new StringBuilder(80);
        StringBuilder t2 = new StringBuilder(80);
        int i = EXTEND, j = EXTEND;
        int matchNum = 0;

        while (!(i == 0 && j == 0)) {
            byte p = ptr[i][j];
            if (p == 0) {
                char ca = a50[i - 1];
                char cb = b50[j - 1];
                t1.append(ca);
                t2.append(cb);
                if (ca == cb) matchNum++;
                i--; j--;
            } else if (p == -1) {
                t1.append('-');
                t2.append(b50[j - 1]);
                j--;
            } else { // up
                t1.append(a50[i - 1]);
                t2.append('-');
                i--;
            }
        }

        t1.reverse();
        t2.reverse();
        return new String[]{t1.toString(), t2.toString(), String.valueOf(score[EXTEND][EXTEND]), String.valueOf(matchNum)};
    }

    // ---------------------------
    // Main
    // ---------------------------
    public static void main(String[] args) {
        if (args.length < 4) {
            System.err.println("Usage: java NewExtendPec50_Optimized <masked_fasta> <unmasked_fasta> <pairs.txt> <out.txt>");
            return;
        }

        try {
            long startTime = System.currentTimeMillis();
            int maxN = 1_000_000;

            // args[0]: masked genome (raw, may contain N / lowercases)
            ArrayList<char[]> segsRaw = readFastaSeqsRaw(args[0]);

            // args[1]: unmasked genome (used as CR reference, stored uppercase)
            ArrayList<char[]> segsCRUpper = readFastaSeqsUpper(args[1]);

            // DP buffers (reused)
            short[][] score = new short[EXTEND + 1][EXTEND + 1];
            byte[][] ptr = new byte[EXTEND + 1][EXTEND + 1];

            // 50bp reusable windows
            char[] win1 = new char[EXTEND];
            char[] win2 = new char[EXTEND];

            int index = 0;
            int indexW = 0;

            try (BufferedReader brPairs = new BufferedReader(new FileReader(args[2]), 1 << 20);
                 BufferedWriter writer = new BufferedWriter(new FileWriter(new File(args[3])), 1 << 20)) {

                String header = brPairs.readLine();
                if (header == null) {
                    System.err.println("Empty pairs file.");
                    return;
                }

                writer.write(header + "  extendS1  extendE1  extendS2  extendE2  cigarS  cigarE  Identity");
                writer.newLine();

                String line;
                while ((line = brPairs.readLine()) != null) {
                    index++;

                    // Parse fields as per your format:
                    // 0 index, 1 cigar, 2 seg1, 3 seg2, 4 chr1, 5 start1, 6 end1, 7 chr2, 8 start2, 9 end2,
                    // 10 Match, 11 Mismatch, 12 nGap, 13 sGap, 14 strand1, 15 strand2
                    int chr1 = parseIntField(line, 4);
                    int chr2 = parseIntField(line, 7);
                    int start1 = parseIntField(line, 5);
                    int end1 = parseIntField(line, 6);
                    int start2 = parseIntField(line, 8);
                    int end2 = parseIntField(line, 9);

                    // Keep original overlap/condition filter
                    if (!(chr1 != chr2 || start1 > end2 || start2 > end1)) {
                        continue;
                    }

                    char strand1 = parseCharField(line, 14);
                    char strand2 = parseCharField(line, 15);

                    String idxField = parseStringField(line, 0);
                    String cigar = parseStringField(line, 1);
                    String seg1Id = parseStringField(line, 2);
                    String seg2Id = parseStringField(line, 3);
                    String chr2Field = parseStringField(line, 7);

                    // Load chromosomes
                    char[] chrRaw = segsRaw.get(chr1);
                    char[] chrCR = segsCRUpper.get(chr1);
                    int chrLen = chrCR.length;

                    char[] otherRaw = segsRaw.get(chr2);
                    char[] otherCR = segsCRUpper.get(chr2);
                    int otherLen = otherCR.length;

                    // Build oriented segments (raw + uppercase CR) for trimming/cigar ops (same as original intent)
                    String seg1 = orientedSegmentRaw(chrRaw, start1, end1, strand1);
                    String segCR1 = orientedSegmentUpper(chrCR, start1, end1, strand1);

                    String seg2 = orientedSegmentRaw(otherRaw, start2, end2, strand2);
                    String segCR2 = orientedSegmentUpper(otherCR, start2, end2, strand2);

                    // Initial alignment per original cigar
                    String[] alignment = cigarAlignment(seg1, seg2, cigar);

                    // Backups
                    int s1 = start1, e1 = end1, s2 = start2, e2 = end2;
                    int extendS1 = -1, extendE1 = -1, extendS2 = -1, extendE2 = -1;

                    String[] cigarS = new String[]{"-1", "-1", "-1", "-1", "-1"};
                    String[] cigarE = new String[]{"-1", "-1", "-1", "-1", "-1"};

                    ArrayList<String> align1 = new ArrayList<>();
                    ArrayList<String> align2 = new ArrayList<>();
                    ArrayList<Double> Pec = new ArrayList<>();

                    // -----------------------
                    // Extend START
                    // -----------------------
                    boolean start = true;
                    boolean change = true;
                    int time1 = 0;
                    align1.clear(); align2.clear(); Pec.clear();

                    while (start) {
                        start = false;
                        change = true;

                        // Determine 50bp upstream windows in oriented space (same as original)
                        // Original:
                        //  strand '+' : window is [start-50, start)
                        //  strand '-' : window is reverse-complement of [end, end+50)
                        int w1s, w1e, w2s, w2e;

                        if (strand1 == '+') { w1s = start1 - EXTEND; w1e = start1; }
                        else { w1s = end1; w1e = end1 + EXTEND; }
                        if (!fillWindow50Upper(chrCR, w1s, w1e, strand1, win1)) change = false;

                        if (change) {
                            if (strand2 == '+') { w2s = start2 - EXTEND; w2e = start2; }
                            else { w2s = end2; w2e = end2 + EXTEND; }
                            if (!fillWindow50Upper(otherCR, w2s, w2e, strand2, win2)) change = false;
                        }

                        if (change) {
                            int[] gate = alignGate50(win1, win2, score, ptr);
                            double idGate = (double) gate[0] / (double) gate[1];

                            if (idGate > 0.5) {
                                start = true;

                                // advance coords like original
                                if (strand1 == '+') start1 -= EXTEND;
                                else end1 += EXTEND;
                                if (strand2 == '+') start2 -= EXTEND;
                                else end2 += EXTEND;

                                if (chr1 != chr2 || start1 > end2 || start2 > end1) {
                                    // Only now compute full alignment strings for cigar building
                                    String[] full = alignFull50(win1, win2, score, ptr);
                                    align1.add(full[0]);
                                    align2.add(full[1]);
                                    Pec.add(Double.parseDouble(full[3]) / (double) full[0].length());
                                } else {
                                    start = false;
                                }
                            }
                        }
                    }

                    if (!align1.isEmpty()) {
                        for (int i = align1.size() - 1; i >= 0 && time1 == 0; i--) {
                            if (Pec.get(i) >= 0.7) {
                                double sumPec = 0;
                                for (int j = 0; j <= i; j++) sumPec += Pec.get(j);
                                if (sumPec / (i + 1) >= 0.7) time1 = i + 1;
                            }
                        }
                    }

                    // restore
                    start1 = s1; end1 = e1; start2 = s2; end2 = e2;

                    // if start can't extend, try trim start (original logic)
                    if (time1 == 0) {
                        change = true;
                        int CR1 = firstCR(seg1);
                        int CR2 = firstCR(seg2);

                        while ((CR1 < EXTEND || CR2 < EXTEND) && change) {
                            if (seg1.length() > EXTEND && seg2.length() > EXTEND) {
                                // build 50bp windows from segCR1/segCR2 prefixes (they are already oriented & uppercase strings)
                                String a = segCR1.substring(0, EXTEND);
                                String b = segCR2.substring(0, EXTEND);

                                // run full alignment on these 50bp strings (small, acceptable)
                                char[] a50 = a.toCharArray();
                                char[] b50 = b.toCharArray();
                                String[] full = alignFull50(a50, b50, score, ptr);
                                double id = Double.parseDouble(full[3]) / (double) full[0].length();

                                if (id < 0.7) {
                                    int[] newIndex;
                                    if (CR1 < EXTEND && CR2 < EXTEND) {
                                        if (CR1 < CR2) newIndex = cutCigarStart(1, CR2, alignment, seg1, String.valueOf(strand1), start1, end1, seg2, String.valueOf(strand2), start2, end2);
                                        else newIndex = cutCigarStart(0, CR1, alignment, seg1, String.valueOf(strand1), start1, end1, seg2, String.valueOf(strand2), start2, end2);
                                    } else if (CR1 < EXTEND) {
                                        newIndex = cutCigarStart(0, CR1, alignment, seg1, String.valueOf(strand1), start1, end1, seg2, String.valueOf(strand2), start2, end2);
                                    } else {
                                        newIndex = cutCigarStart(1, CR2, alignment, seg1, String.valueOf(strand1), start1, end1, seg2, String.valueOf(strand2), start2, end2);
                                    }

                                    start1 = newIndex[0]; end1 = newIndex[1];
                                    start2 = newIndex[2]; end2 = newIndex[3];

                                    // rebuild segments after trimming
                                    seg1 = orientedSegmentRaw(chrRaw, start1, end1, strand1);
                                    segCR1 = orientedSegmentUpper(chrCR, start1, end1, strand1);
                                    seg2 = orientedSegmentRaw(otherRaw, start2, end2, strand2);
                                    segCR2 = orientedSegmentUpper(otherCR, start2, end2, strand2);

                                    CR1 = firstCR(seg1);
                                    CR2 = firstCR(seg2);
                                    time1++;
                                } else {
                                    change = false;
                                }
                            } else {
                                change = false;
                            }
                        }

                        // if start trimmed, re-run start extension on trimmed coords (original structure)
                        if (time1 > 0) {
                            s1 = start1; e1 = end1; s2 = start2; e2 = end2;
                            align1.clear(); align2.clear(); Pec.clear();
                            start = true;
                            change = true;
                            int time3 = 0;

                            while (start) {
                                start = false;
                                change = true;

                                int w1s, w1e, w2s, w2e;
                                if (strand1 == '+') { w1s = start1 - EXTEND; w1e = start1; }
                                else { w1s = end1; w1e = end1 + EXTEND; }
                                if (!fillWindow50Upper(chrCR, w1s, w1e, strand1, win1)) change = false;

                                if (change) {
                                    if (strand2 == '+') { w2s = start2 - EXTEND; w2e = start2; }
                                    else { w2s = end2; w2e = end2 + EXTEND; }
                                    if (!fillWindow50Upper(otherCR, w2s, w2e, strand2, win2)) change = false;
                                }

                                if (change) {
                                    int[] gate = alignGate50(win1, win2, score, ptr);
                                    double idGate = (double) gate[0] / (double) gate[1];

                                    if (idGate > 0.5) {
                                        start = true;
                                        if (strand1 == '+') start1 -= EXTEND; else end1 += EXTEND;
                                        if (strand2 == '+') start2 -= EXTEND; else end2 += EXTEND;

                                        if (chr1 != chr2 || start1 > end2 || start2 > end1) {
                                            String[] full = alignFull50(win1, win2, score, ptr);
                                            align1.add(full[0]);
                                            align2.add(full[1]);
                                            Pec.add(Double.parseDouble(full[3]) / (double) full[0].length());
                                        } else {
                                            start = false;
                                        }
                                    }
                                }
                            }

                            if (!align1.isEmpty()) {
                                for (int i = align1.size() - 1; i >= 0 && time3 == 0; i--) {
                                    if (Pec.get(i) >= 0.7) {
                                        double sumPec = 0;
                                        for (int j = 0; j <= i; j++) sumPec += Pec.get(j);
                                        if (sumPec / (i + 1) >= 0.7) time3 = i + 1;
                                    }
                                }
                            }

                            start1 = s1; end1 = e1; start2 = s2; end2 = e2;
                            if (time3 > 0) {
                                while (align1.size() > time3) { align1.remove(align1.size() - 1); align2.remove(align2.size() - 1); }

                                if (strand1 == '+') extendS1 = s1 - EXTEND * time3;
                                else extendE1 = e1 + EXTEND * time3;

                                if (strand2 == '+') extendS2 = s2 - EXTEND * time3;
                                else extendE2 = e2 + EXTEND * time3;

                                cigarS = updateAlign(0, align1, align2);
                            }
                        }
                    } else {
                        // extend start directly
                        while (align1.size() > time1) { align1.remove(align1.size() - 1); align2.remove(align2.size() - 1); }

                        if (strand1 == '+') extendS1 = s1 - EXTEND * time1;
                        else extendE1 = e1 + EXTEND * time1;

                        if (strand2 == '+') extendS2 = s2 - EXTEND * time1;
                        else extendE2 = e2 + EXTEND * time1;

                        cigarS = updateAlign(0, align1, align2);
                    }

                    // -----------------------
                    // Extend END
                    // -----------------------
                    s1 = start1; e1 = end1; s2 = start2; e2 = end2;
                    align1.clear(); align2.clear(); Pec.clear();

                    boolean end = true;
                    change = true;
                    int time2 = 0;

                    // apply start extensions to coords before extending end (same as original)
                    if (extendS1 != -1) start1 = extendS1;
                    if (extendS2 != -1) start2 = extendS2;
                    if (extendE1 != -1) end1 = extendE1;
                    if (extendE2 != -1) end2 = extendE2;

                    while (end) {
                        end = false;
                        change = true;

                        // Original:
                        //  strand '+' : window is [end, end+50)
                        //  strand '-' : window is reverse-complement of [start-50, start)
                        int w1s, w1e, w2s, w2e;

                        if (strand1 == '+') { w1s = end1; w1e = end1 + EXTEND; }
                        else { w1s = start1 - EXTEND; w1e = start1; }
                        if (!fillWindow50Upper(chrCR, w1s, w1e, strand1, win1)) change = false;

                        if (change) {
                            if (strand2 == '+') { w2s = end2; w2e = end2 + EXTEND; }
                            else { w2s = start2 - EXTEND; w2e = start2; }
                            if (!fillWindow50Upper(otherCR, w2s, w2e, strand2, win2)) change = false;
                        }

                        if (change) {
                            int[] gate = alignGate50(win1, win2, score, ptr);
                            double idGate = (double) gate[0] / (double) gate[1];

                            if (idGate > 0.5) {
                                end = true;
                                if (strand1 == '+') end1 += EXTEND; else start1 -= EXTEND;
                                if (strand2 == '+') end2 += EXTEND; else start2 -= EXTEND;

                                if (chr1 != chr2 || start1 > end2 || start2 > end1) {
                                    String[] full = alignFull50(win1, win2, score, ptr);
                                    align1.add(full[0]);
                                    align2.add(full[1]);
                                    Pec.add(Double.parseDouble(full[3]) / (double) full[0].length());
                                } else {
                                    end = false;
                                }
                            }
                        }
                    }

                    if (!align1.isEmpty()) {
                        for (int i = align1.size() - 1; i >= 0 && time2 == 0; i--) {
                            if (Pec.get(i) >= 0.7) {
                                double sumPec = 0;
                                for (int j = 0; j <= i; j++) sumPec += Pec.get(j);
                                if (sumPec / (i + 1) >= 0.7) time2 = i + 1;
                            }
                        }
                    }

                    // restore to pre-end-extension coords (original behavior)
                    start1 = s1; end1 = e1; start2 = s2; end2 = e2;

                    // If end can't extend, trim end (original structure)
                    if (time2 == 0) {
                        int CR1 = lastCR(seg1);
                        int CR2 = lastCR(seg2);
                        change = true;

                        while ((CR1 < EXTEND || CR2 < EXTEND) && change) {
                            if (seg1.length() > EXTEND && seg2.length() > EXTEND) {
                                String a = segCR1.substring(segCR1.length() - EXTEND);
                                String b = segCR2.substring(segCR2.length() - EXTEND);

                                char[] a50 = a.toCharArray();
                                char[] b50 = b.toCharArray();
                                String[] full = alignFull50(a50, b50, score, ptr);
                                double id = Double.parseDouble(full[3]) / (double) full[0].length();

                                if (id < 0.7) {
                                    int[] newIndex;
                                    if (CR1 < EXTEND && CR2 < EXTEND) {
                                        if (CR1 < CR2) newIndex = cutCigarEnd(1, CR2, alignment, seg1, String.valueOf(strand1), start1, end1, seg2, String.valueOf(strand2), start2, end2);
                                        else newIndex = cutCigarEnd(0, CR1, alignment, seg1, String.valueOf(strand1), start1, end1, seg2, String.valueOf(strand2), start2, end2);
                                    } else if (CR1 < EXTEND) {
                                        newIndex = cutCigarEnd(0, CR1, alignment, seg1, String.valueOf(strand1), start1, end1, seg2, String.valueOf(strand2), start2, end2);
                                    } else {
                                        newIndex = cutCigarEnd(1, CR2, alignment, seg1, String.valueOf(strand1), start1, end1, seg2, String.valueOf(strand2), start2, end2);
                                    }

                                    start1 = newIndex[0]; end1 = newIndex[1];
                                    start2 = newIndex[2]; end2 = newIndex[3];

                                    seg1 = orientedSegmentRaw(chrRaw, start1, end1, strand1);
                                    segCR1 = orientedSegmentUpper(chrCR, start1, end1, strand1);
                                    seg2 = orientedSegmentRaw(otherRaw, start2, end2, strand2);
                                    segCR2 = orientedSegmentUpper(otherCR, start2, end2, strand2);

                                    CR1 = lastCR(seg1);
                                    CR2 = lastCR(seg2);
                                    time2++;
                                } else {
                                    change = false;
                                }
                            } else {
                                change = false;
                            }
                        }

                        // extend end after trim (original structure)
                        if (time2 > 0) {
                            s1 = start1; e1 = end1; s2 = start2; e2 = end2;
                            align1.clear(); align2.clear(); Pec.clear();

                            end = true;
                            change = true;
                            int time3 = 0;

                            // re-apply start extensions before extending end
                            if (extendS1 != -1) start1 = extendS1;
                            if (extendS2 != -1) start2 = extendS2;
                            if (extendE1 != -1) end1 = extendE1;
                            if (extendE2 != -1) end2 = extendE2;

                            while (end) {
                                end = false;
                                change = true;

                                int w1s, w1e, w2s, w2e;
                                if (strand1 == '+') { w1s = end1; w1e = end1 + EXTEND; }
                                else { w1s = start1 - EXTEND; w1e = start1; }
                                if (!fillWindow50Upper(chrCR, w1s, w1e, strand1, win1)) change = false;

                                if (change) {
                                    if (strand2 == '+') { w2s = end2; w2e = end2 + EXTEND; }
                                    else { w2s = start2 - EXTEND; w2e = start2; }
                                    if (!fillWindow50Upper(otherCR, w2s, w2e, strand2, win2)) change = false;
                                }

                                if (change) {
                                    int[] gate = alignGate50(win1, win2, score, ptr);
                                    double idGate = (double) gate[0] / (double) gate[1];

                                    if (idGate > 0.5) {
                                        end = true;
                                        if (strand1 == '+') end1 += EXTEND; else start1 -= EXTEND;
                                        if (strand2 == '+') end2 += EXTEND; else start2 -= EXTEND;

                                        if (chr1 != chr2 || start1 > end2 || start2 > end1) {
                                            String[] full = alignFull50(win1, win2, score, ptr);
                                            align1.add(full[0]);
                                            align2.add(full[1]);
                                            Pec.add(Double.parseDouble(full[3]) / (double) full[0].length());
                                        } else {
                                            end = false;
                                        }
                                    }
                                }
                            }

                            if (!align1.isEmpty()) {
                                for (int i = align1.size() - 1; i >= 0 && time3 == 0; i--) {
                                    if (Pec.get(i) >= 0.7) {
                                        double sumPec = 0;
                                        for (int j = 0; j <= i; j++) sumPec += Pec.get(j);
                                        if (sumPec / (i + 1) >= 0.7) time3 = i + 1;
                                    }
                                }
                            }

                            start1 = s1; end1 = e1; start2 = s2; end2 = e2;
                            if (time3 > 0) {
                                while (align1.size() > time3) { align1.remove(align1.size() - 1); align2.remove(align2.size() - 1); }

                                if (strand1 == '+') extendE1 = e1 + EXTEND * time3;
                                else extendS1 = s1 - EXTEND * time3;

                                if (strand2 == '+') extendE2 = e2 + EXTEND * time3;
                                else extendS2 = s2 - EXTEND * time3;

                                cigarE = updateAlign(1, align1, align2);
                            }
                        }
                    } else {
                        // extend end directly
                        while (align1.size() > time2) { align1.remove(align1.size() - 1); align2.remove(align2.size() - 1); }

                        if (strand1 == '+') extendE1 = e1 + EXTEND * time2;
                        else extendS1 = s1 - EXTEND * time2;

                        if (strand2 == '+') extendE2 = e2 + EXTEND * time2;
                        else extendS2 = s2 - EXTEND * time2;

                        cigarE = updateAlign(1, align1, align2);
                    }

                    // -----------------------
                    // Final identity and output (same as original)
                    // -----------------------
                    String[] newCigar = computeCigar(alignment);

                    String finalSeg1;
                    int finalM, finalMis, finalNgap, finalSgap;

                    if (extendS1 == -1 && extendE1 == -1) {
                        finalSeg1 = orientedSegmentRaw(chrRaw, start1, end1, '+'); // output uses raw forward segment; original uses segs.get(curChr).substring(...)
                        finalM = Integer.parseInt(newCigar[1]);
                        finalMis = Integer.parseInt(newCigar[2]);
                        finalNgap = Integer.parseInt(newCigar[3]);
                        finalSgap = Integer.parseInt(newCigar[4]);
                    } else if (extendS1 == -1) {
                        finalSeg1 = new String(chrRaw, start1, extendE1 - start1);
                        finalM = Integer.parseInt(newCigar[1]) + Integer.parseInt(cigarE[1]);
                        finalMis = Integer.parseInt(newCigar[2]) + Integer.parseInt(cigarE[2]);
                        finalNgap = Integer.parseInt(newCigar[3]) + Integer.parseInt(cigarE[3]);
                        finalSgap = Integer.parseInt(newCigar[4]) + Integer.parseInt(cigarE[4]);
                    } else if (extendE1 == -1) {
                        finalSeg1 = new String(chrRaw, extendS1, end1 - extendS1);
                        finalM = Integer.parseInt(newCigar[1]) + Integer.parseInt(cigarS[1]);
                        finalMis = Integer.parseInt(newCigar[2]) + Integer.parseInt(cigarS[2]);
                        finalNgap = Integer.parseInt(newCigar[3]) + Integer.parseInt(cigarS[3]);
                        finalSgap = Integer.parseInt(newCigar[4]) + Integer.parseInt(cigarS[4]);
                    } else {
                        finalSeg1 = new String(chrRaw, extendS1, extendE1 - extendS1);
                        finalM = Integer.parseInt(newCigar[1]) + Integer.parseInt(cigarS[1]) + Integer.parseInt(cigarE[1]);
                        finalMis = Integer.parseInt(newCigar[2]) + Integer.parseInt(cigarS[2]) + Integer.parseInt(cigarE[2]);
                        finalNgap = Integer.parseInt(newCigar[3]) + Integer.parseInt(cigarS[3]) + Integer.parseInt(cigarE[3]);
                        finalSgap = Integer.parseInt(newCigar[4]) + Integer.parseInt(cigarS[4]) + Integer.parseInt(cigarE[4]);
                    }

                    double identity = (double) finalM / (double) (finalM + finalMis + finalSgap);

                    String finalSeg2;
                    if (extendS2 == -1 && extendE2 == -1) {
                        finalSeg2 = new String(otherRaw, start2, end2 - start2);
                    } else if (extendS2 == -1) {
                        finalSeg2 = new String(otherRaw, start2, extendE2 - start2);
                    } else if (extendE2 == -1) {
                        finalSeg2 = new String(otherRaw, extendS2, end2 - extendS2);
                    } else {
                        finalSeg2 = new String(otherRaw, extendS2, extendE2 - extendS2);
                    }

                    int numCR1 = Nlength(finalSeg1);
                    int numCR2 = Nlength(finalSeg2);

                    if (finalSeg1.length() - numCR1 >= 500
                            && finalSeg2.length() - numCR2 >= 500
                            && identity >= 0.7
                            && Math.abs(finalSeg1.length() - finalSeg2.length()) < maxN) {

                        writer.write(idxField + "  " + newCigar[0] + "  " + seg1Id + "  " + seg2Id + "  " + chr1 + "  ");
                        writer.write(start1 + "  " + end1 + "  " + chr2Field + "  " + start2 + "  " + end2 + "  ");
                        writer.write(finalM + "  " + finalMis + "  " + finalNgap + "  " + finalSgap + "  " + strand1 + "  " + strand2 + "  ");
                        writer.write(extendS1 + "  " + extendE1 + "  " + extendS2 + "  " + extendE2 + "  " + cigarS[0] + "  " + cigarE[0] + "  " + identity);
                        writer.newLine();
                        indexW++;
                    }
                }
            }

            long endTime = System.currentTimeMillis();
            System.out.println("end!  totalNum: " + index + "   " + indexW);
            System.out.println("runTime: " + (endTime - startTime));

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
