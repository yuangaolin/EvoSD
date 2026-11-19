import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

public class RecoverFast {

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

        while ((line = br.readLine()) != null) {

            String[] p = line.split("\\s+");

            if (p[16].charAt(0) == '1') {
                sb.append(p[0]).append(SEP).append(p[1]).append(SEP)
                  .append(p[2]).append(SEP).append(p[3]).append(SEP)
                  .append(p[7]).append(SEP).append(p[8]).append(SEP)
                  .append(p[9]).append(SEP)
                  .append(p[4]).append(SEP).append(p[5]).append(SEP)
                  .append(p[6]).append(SEP)
                  .append(p[10]).append(SEP).append(p[11]).append(SEP)
                  .append(p[12]).append(SEP).append(p[13]).append(SEP)
                  .append(p[14]).append(SEP).append(p[15]);
            } else {
                for (int i = 0; i < 16; i++) {
                    sb.append(p[i]).append(SEP);
                }
                sb.setLength(sb.length() - SEP.length()); // 去掉末尾多余空格
            }

            sb.append('\n');

            if (sb.length() > (1 << 20)) {
                writer.write(sb.toString());
                sb.setLength(0);
            }
        }

        if (sb.length() > 0) writer.write(sb.toString());

        writer.close();
        br.close();

        long endTime = System.currentTimeMillis();
        System.out.println("runTime:" + (endTime - startTime));
    }
}
