import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

public class ChangeFormatMulti {

    // 用于跨线程共享 genomeSize
    static class GenomeStats {
        double genomeSize = 0;
        int CRlength = 0;
        int Nlength = 0;
    }

    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        if (args.length < 3) {
            System.out.println("Usage: java ChangeFormatMultithread genome.fa genome.masked.fa outdir");
            return;
        }

        String genomeFile = args[0];
        String maskedFile = args[1];
        String outdir = args[2];

        File outDirFile = new File(outdir);
        if (!outDirFile.exists()) outDirFile.mkdirs();

        

        GenomeStats stats = new GenomeStats();

        // 线程1：处理 genome.fa
        Thread t1 = new Thread(() -> {
            processGenome(genomeFile, outdir, stats);
        });

        // 线程2：处理 genome.masked.fa
        Thread t2 = new Thread(() -> {
            processMasked(maskedFile, outdir);
        });

        t1.start();
        t2.start();

        try {
            t1.join();
            t2.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        long endTime = System.currentTimeMillis();
        System.out.println("The Input Genome Size is: " + stats.genomeSize);
        System.out.println("runTime: " + (endTime - startTime) + " ms");
    }

    private static void processGenome(String genomeFile, String outdir, GenomeStats stats) {
        try (Scanner in = new Scanner(new File(genomeFile));
             BufferedWriter writer = new BufferedWriter(new FileWriter(new File(outdir, "genome.txt")));
             BufferedWriter writerGenome = new BufferedWriter(new FileWriter(new File(outdir, "Genome_size_Indexes.txt")))) {

            double genomeSize = 0;
            int segNum = 0;
            int CRlength = 0;
            int Nlength = 0;

            String temp = in.nextLine();
            if (temp.charAt(0) != '>') {
                System.out.println("The first line of fa genome data is not start with character > ");
            } else {
                writer.write(">" + segNum);
                writer.newLine();
                writerGenome.write(temp);
                writerGenome.newLine();
                segNum++;
            }

            while (in.hasNextLine()) {
                temp = in.nextLine();
                if (temp.charAt(0) != '>') {
                    writer.write(temp);
                    genomeSize += temp.length();
                    for (int i = 0; i < temp.length(); i++) {
                        char c = temp.charAt(i);
                        if (c == 'a' || c == 'c' || c == 'g' || c == 't' ||
                            c == 'A' || c == 'C' || c == 'G' || c == 'T') {
                            CRlength++;
                        } else if (c == 'n' || c == 'N') {
                            Nlength++;
                        }
                    }
                } else {
                    writer.newLine();
                    writer.write(">" + segNum);
                    writer.newLine();
                    writerGenome.write(temp);
                    writerGenome.newLine();
                    segNum++;
                }
            }

            writer.newLine();
            writerGenome.write("GenomeSize:  " + genomeSize + "    segNum: " + segNum);
            writerGenome.newLine();
            writerGenome.write("CRlength:  " + CRlength + "    CR/Genome: " + (CRlength / genomeSize));
            writerGenome.newLine();
            writerGenome.write("Nlength:  " + Nlength + "    N/Genome: " + (Nlength / genomeSize));
            writerGenome.newLine();
            writerGenome.write("(CR+N)/Genome: " + ((CRlength + Nlength) / genomeSize));
            writerGenome.newLine();

            // 将统计结果写入共享对象
            stats.genomeSize = genomeSize;
            stats.CRlength = CRlength;
            stats.Nlength = Nlength;

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void processMasked(String maskedFile, String outdir) {
        try (Scanner in = new Scanner(new File(maskedFile));
             BufferedWriter writer = new BufferedWriter(new FileWriter(new File(outdir, "genome.masked.txt")))) {

            int segNum = 0;
            String temp = in.nextLine();
            if (temp.charAt(0) != '>') {
                System.out.println("The first line of masked genome data is not start with character > ");
            } else {
                writer.write(">" + segNum);
                writer.newLine();
                segNum++;
            }

            while (in.hasNextLine()) {
                temp = in.nextLine();
                if (temp.charAt(0) == '>') {
                    writer.newLine();
                    writer.write(">" + segNum);
                    writer.newLine();
                    segNum++;
                } else {
                    writer.write(temp);
                }
            }

            writer.newLine();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
