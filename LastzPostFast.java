import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;

public class LastzPostFast {

    static int gapLengthFast(String cigar) {
        int num = 0, sum = 0;
        for (int i = 0; i < cigar.length(); i++) {
            char c = cigar.charAt(i);
            if (c >= '0' && c <= '9') {
                num = num * 10 + (c - '0');
            } else {
                if (c == 'I' || c == 'D') sum += num;
                num = 0;
            }
        }
        return sum;
    }

    static boolean wrongGapFast(String s) {
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            if (!(c == 'M' || c == 'I' || c == 'D' || (c >= '0' && c <= '9')))
                return true;
        }
        return false;
    }

    static int[][] buildIndexMaps(ArrayList<String> segs) {
        int[][] maps = new int[segs.size()][];
        for (int i = 0; i < segs.size(); i++) {
            String s = segs.get(i);
            int[] tmp = new int[s.length()];
            int cnt = 0;
            for (int j = 0; j < s.length(); j++) {
                if (s.charAt(j) != 'N') tmp[cnt++] = j;
            }
            maps[i] = Arrays.copyOf(tmp, cnt);
        }
        return maps;
    }

    public static void main(String[] args) throws Exception {

        long start = System.currentTimeMillis();

        BufferedWriter writer = new BufferedWriter(new FileWriter(args[2]), 1 << 20);
        writer.write("index  cigar seg1  seg2  chr1  start1  end1  chr2  start2  end2  Match  Mismatch  nGap  sGap  strand1  strand2\n");

        ArrayList<String> segs = new ArrayList<>();
        ArrayList<int[]> segIndex = new ArrayList<>();

        BufferedReader brSeg = new BufferedReader(new FileReader(args[0]), 1 << 20);
        String line;
        while ((line = brSeg.readLine()) != null) {
            String[] sp = line.trim().split("\\s+");
            segIndex.add(new int[]{Integer.parseInt(sp[1]), Integer.parseInt(sp[2]), Integer.parseInt(sp[3])});
            segs.add(brSeg.readLine());
        }
        brSeg.close();

        int[][] indexMap = buildIndexMaps(segs);

        BufferedReader br = new BufferedReader(new FileReader(args[1]), 1 << 20);
        br.readLine(); // skip header

        StringBuilder sb = new StringBuilder(1 << 20);
        int index = 0, wrong = 0;

        while ((line = br.readLine()) != null) {
            String[] p = line.split("\\s+");
            if (p.length != 20 || p[0].equalsIgnoreCase("#cigar")) continue;

            String cigar = p[0];
            if (wrongGapFast(cigar)) {
                wrong++;
                continue;
            }

            int seg1 = Integer.parseInt(p[2]);
            int seg2 = Integer.parseInt(p[7]);

            if (seg1 != seg2 || Integer.parseInt(p[10]) > Integer.parseInt(p[6]) ||
                Integer.parseInt(p[5]) > Integer.parseInt(p[11])) {

                int[] map1 = indexMap[seg1];
                int[] map2 = indexMap[seg2];

                sb.append(index++).append("  ").append(cigar).append("  ")
                  .append(seg1).append("  ").append(seg2).append("  ");

                if (p[3].equals("+")) {
                    sb.append(segIndex.get(seg1)[0]).append("  ")
                      .append(segIndex.get(seg1)[1] + map1[Integer.parseInt(p[5])]).append("  ")
                      .append(segIndex.get(seg1)[1] + map1[Integer.parseInt(p[6]) - 1] + 1).append("  ");
                } else {
                    int len = map1.length;
                    sb.append(segIndex.get(seg1)[0]).append("  ")
                      .append(segIndex.get(seg1)[1] + map1[len - Integer.parseInt(p[6])]).append("  ")
                      .append(segIndex.get(seg1)[1] + map1[len - Integer.parseInt(p[5]) - 1] + 1).append("  ");
                }

                if (p[8].equals("+")) {
                    sb.append(segIndex.get(seg2)[0]).append("  ")
                      .append(segIndex.get(seg2)[1] + map2[Integer.parseInt(p[10])]).append("  ")
                      .append(segIndex.get(seg2)[1] + map2[Integer.parseInt(p[11]) - 1] + 1).append("  ");
                } else {
                    int len = map2.length;
                    sb.append(segIndex.get(seg2)[0]).append("  ")
                      .append(segIndex.get(seg2)[1] + map2[len - Integer.parseInt(p[11])]).append("  ")
                      .append(segIndex.get(seg2)[1] + map2[len - Integer.parseInt(p[10]) - 1] + 1).append("  ");
                }

                sb.append(p[16]).append("  ")
                  .append(p[17]).append("  ")
                  .append(p[18]).append("  ")
                  .append(gapLengthFast(cigar)).append("  ")
                  .append(p[3]).append("  ")
                  .append(p[8]).append('\n');

                if (sb.length() > (1 << 20)) {
                    writer.write(sb.toString());
                    sb.setLength(0);
                }
            }
        }

        if (sb.length() > 0) writer.write(sb.toString());
        writer.close();
        br.close();

        System.out.println("end! " + index + " wrongGap: " + wrong +
                " runTime: " + (System.currentTimeMillis() - start));
    }
}
