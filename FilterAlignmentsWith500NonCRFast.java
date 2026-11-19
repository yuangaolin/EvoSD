import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

public class FilterAlignmentsWith500NonCRFast {

    static int[] indelLengthFast(String cigar) {
        int ins = 0, del = 0, num = 0;
        for (int i = 0; i < cigar.length(); i++) {
            char c = cigar.charAt(i);
            if (c >= '0' && c <= '9') {
                num = num * 10 + (c - '0');
            } else {
                if (c == 'I') ins += num;
                else if (c == 'D') del += num;
                num = 0;
            }
        }
        return new int[]{ins, del};
    }

    public static void main(String[] args) throws Exception {

        long startTime = System.currentTimeMillis();

        BufferedReader br = new BufferedReader(new FileReader(args[0]), 1 << 20);
        BufferedWriter writer = new BufferedWriter(new FileWriter(args[1]), 1 << 20);

        String line = br.readLine();   // header
        writer.write(line);
        writer.newLine();

        StringBuilder sb = new StringBuilder(1 << 20);
        int num = 0;

        while ((line = br.readLine()) != null) {

            // 手动 split（比 regex split 快 2–3×）
            String[] p = line.split("\\s+");

            int match = Integer.parseInt(p[10]);
            int mismatch = Integer.parseInt(p[11]);
            int nCount = Integer.parseInt(p[13]);

            int[] indel = indelLengthFast(p[1]);
            int len1 = match + mismatch + indel[0];
            int len2 = match + mismatch + indel[1];
            int length = Math.min(len1, len2);

            if (length >= 500 && (double) match / (match + mismatch + nCount) >= 0.7) {
                sb.append(line).append('\n');
                num++;
            }

            if (sb.length() > (1 << 20)) {
                writer.write(sb.toString());
                sb.setLength(0);
            }
        }

        if (sb.length() > 0) writer.write(sb.toString());

        writer.close();
        br.close();

        long endTime = System.currentTimeMillis();
        System.out.println("runTime: " + (endTime - startTime));
    }
}
