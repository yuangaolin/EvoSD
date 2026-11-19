import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

public class FilterRepetitivePairsFast {

    static final String SEP = "  ";

    public static void main(String[] args) throws Exception {

        long startTime = System.currentTimeMillis();

        BufferedReader br = new BufferedReader(new FileReader(args[0]), 1 << 20);
        BufferedWriter writer = new BufferedWriter(new FileWriter(args[1]), 1 << 20);

        String header = br.readLine();
        writer.write(header);
        writer.newLine();

        StringBuilder sb = new StringBuilder(1 << 20);

        String line = br.readLine();
        if (line == null) return;

        String[] prev = line.split("\\s+");
        sb.append(line).append('\n');

        long index = 1, indexR = 0;

        while ((line = br.readLine()) != null) {

            String[] cur = line.split("\\s+");

            if (cur[4].equals(prev[4]) &&
                cur[5].equals(prev[5]) &&
                cur[6].equals(prev[6]) &&
                cur[7].equals(prev[7]) &&
                cur[8].equals(prev[8]) &&
                cur[9].equals(prev[9]) &&
                cur[10].equals(prev[10]) &&
                cur[11].equals(prev[11]) &&
                cur[13].equals(prev[13])) {

                indexR++;
            } else {
                sb.append(line).append('\n');
                prev = cur;
            }

            index++;

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
