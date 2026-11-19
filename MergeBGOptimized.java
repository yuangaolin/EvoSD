import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;

public class MergeBGOptimized {

    private static int nextNonSpace(String s, int i) {
        int n = s.length();
        while (i < n && s.charAt(i) <= ' ') i++;
        return i;
    }

    private static int nextSpace(String s, int i) {
        int n = s.length();
        while (i < n && s.charAt(i) > ' ') i++;
        return i;
    }

    private static int parseIntAt(String s, int start, int end) {
        int i = start;
        int sign = 1;
        if (i < end && s.charAt(i) == '-') { // not expected but safe
            sign = -1;
            i++;
        }
        int val = 0;
        for (; i < end; i++) {
            val = val * 10 + (s.charAt(i) - '0');
        }
        return val * sign;
    }

    // parse three ints: chr start end (space separated)
    private static void parse3Ints(String line, int[] out) {
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

    public static void main(String[] args) {
        try {
            long startTime = System.currentTimeMillis();

            try (BufferedReader br = new BufferedReader(new FileReader(args[0]), 1 << 20);
                 BufferedWriter bw = new BufferedWriter(new FileWriter(args[1]), 1 << 20)) {

                bw.write("chr  start  end");
                bw.newLine();

                String line = br.readLine();
                if (line == null) return;

                int[] v = new int[3];
                parse3Ints(line, v);
                int chr = v[0], start = v[1], end = v[2];
                int num = 0;

                while ((line = br.readLine()) != null) {
                    parse3Ints(line, v);
                    int c = v[0], s = v[1], e = v[2];

                    // keep your original rule: ignore if start==end
                    if (s == e) continue;

                    if (c == chr) {
                        if (s > end) {
                            bw.write(chr + "  " + start + "  " + end);
                            bw.newLine();
                            num++;
                            start = s;
                            end = e;
                        } else {
                            if (e > end) end = e;
                        }
                    } else {
                        bw.write(chr + "  " + start + "  " + end);
                        bw.newLine();
                        num++;
                        chr = c;
                        start = s;
                        end = e;
                    }
                }

                bw.write(chr + "  " + start + "  " + end);
                bw.newLine();
                num++;
                // 如需输出 num，可打开：
                // System.out.println("Merged intervals: " + num);
            }

            long endTime = System.currentTimeMillis();
            System.out.println("runTime:" + (endTime - startTime));

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
