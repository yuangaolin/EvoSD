import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

public class SortPairsFast {

    static final String SEP = "  ";

    public static void main(String[] args) throws Exception {

        long startTime = System.currentTimeMillis();

        BufferedReader br = new BufferedReader(new FileReader(args[0]), 1 << 20);
        BufferedWriter writer = new BufferedWriter(new FileWriter(args[1]), 1 << 20);

        String header = br.readLine();
        writer.write(header);
        writer.newLine();

        StringBuilder sb = new StringBuilder(1 << 20);

        String line;
        int index = 0, indexR = 0;

        while ((line = br.readLine()) != null) {

            String[] p = line.split("\\s+");

            int a4 = Integer.parseInt(p[4]);
            int a7 = Integer.parseInt(p[7]);
            int a6 = Integer.parseInt(p[6]);
            int a8 = Integer.parseInt(p[8]);
            int a5 = Integer.parseInt(p[5]);
            int a9 = Integer.parseInt(p[9]);

            index++;

            if (a4 < a7 || (a4 == a7 && a6 < a8)) {
                sb.append(line).append(SEP).append('0').append('\n');
            } else {
                indexR++;
                sb.append(p[0]).append(SEP).append(p[1]).append(SEP)
                  .append(p[2]).append(SEP).append(p[3]).append(SEP)
                  .append(p[7]).append(SEP).append(p[8]).append(SEP)
                  .append(p[9]).append(SEP)
                  .append(p[4]).append(SEP).append(p[5]).append(SEP)
                  .append(p[6]).append(SEP)
                  .append(p[10]).append(SEP).append(p[11]).append(SEP)
                  .append(p[12]).append(SEP).append(p[13]).append(SEP)
                  .append(p[14]).append(SEP).append(p[15]).append(SEP)
                  .append('1').append('\n');
            }

            if (sb.length() > (1 << 20)) {
                writer.write(sb.toString());
                sb.setLength(0);
            }
        }

        if (sb.length() > 0) writer.write(sb.toString());

        br.close();
        writer.close();

        long endTime = System.currentTimeMillis();
        System.out.println("runTime:" + (endTime - startTime));
    }
}
