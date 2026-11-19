import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class CombineLineMulti {

    static class TaskResult {
        int index;
        String result;
        TaskResult(int i, String r) {
            index = i;
            result = r;
        }
    }

    public static void main(String[] args) {
        try {
            long startTime = System.currentTimeMillis();

            double denseBound = 0.01;
            int combineLength = 1000;
            double lengthBound = 1000;

            String freqFile = args[0];
            String genomeFile = args[1];
            String outFile = args[2];
            int threads = args.length >= 4 ? Integer.parseInt(args[3])
                    : Runtime.getRuntime().availableProcessors();

            Scanner in = new Scanner(new File(freqFile));
            Scanner inGenome = new Scanner(new File(genomeFile));
            BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

            ExecutorService pool = Executors.newFixedThreadPool(threads);
            List<Future<TaskResult>> futures = new ArrayList<>();

            int lineIndex = 0;

            while (in.hasNextLine()) {
                String header = in.nextLine();     // >id
                String freq = in.nextLine();       // freq string
                inGenome.nextLine();               // skip genome header
                String genome = inGenome.nextLine();

                final int taskId = lineIndex++;
                futures.add(pool.submit(() ->
                        new TaskResult(taskId,
                                processOneChromosome(header, freq, genome,
                                        denseBound, combineLength, lengthBound))
                ));
            }

            pool.shutdown();

            // 按输入顺序写出
            futures.sort(Comparator.comparingInt(a -> {
                try {
                    return a.get().index;
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }));

            for (Future<TaskResult> f : futures) {
                writer.write(f.get().result);
                writer.newLine();
            }

            writer.close();
            in.close();
            inGenome.close();

            long endTime = System.currentTimeMillis();
            System.out.println(" end Combine Line! ");
            System.out.println("runTime:" + (endTime - startTime));

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    // ================== 完整复制你原来的算法 ==================
    private static String processOneChromosome(String header, String freq, String genome,
                                               double denseBound, int combineLength,
                                               double lengthBound) {

        StringBuilder out = new StringBuilder();
        out.append(header).append("\n");

        int start, end, dis, kmerNum, length, index;
        double totalDis;
        boolean find, findN;
        char curFreq, curGenome;
        int lastIndex, nextIndex;

        length = freq.length();
        start = 0;
        end = 0;
        totalDis = 0;
        kmerNum = 0;
        index = 0;
        find = false;

        while (index < length && (!find)) {
            curFreq = freq.charAt(index);
            curGenome = genome.charAt(index);
            if ((curFreq != '0') && (curFreq != '1') &&
                    (curGenome != 'N') && (curGenome != 'n')) {
                find = true;
                start = index;
                end = start;
                totalDis = 1;
                kmerNum = 1;
            } else {
                out.append("0");
            }
            index++;
        }

        nextIndex = -1;
        lastIndex = -1;

        while (index < length) {

            find = false;
            findN = false;
            dis = 1;
            while (index < length && (!find) && (!findN)) {
                curFreq = freq.charAt(index);
                curGenome = genome.charAt(index);
                index++;
                if ((curFreq != '0') && (curFreq != '1')) {
                    find = true;
                } else if (curGenome == 'N' || curGenome == 'n') {
                    findN = true;
                } else {
                    if (curFreq == '1') {
                        dis++;
                    }
                }
            }

            if (findN == true) {
                if (totalDis >= lengthBound) {
                    if (lastIndex != -1) {
                        for (int i = start; i < lastIndex; i++) {
                            out.append("1");
                        }
                        index = lastIndex;
                        while ((index < length) &&
                                (freq.charAt(index) == '0' || freq.charAt(index) == '1')) {
                            index++;
                            out.append("0");
                        }
                        start = index;
                        end = start;
                        index++;
                    } else {
                        for (int i = start; i < nextIndex; i++) {
                            out.append("0");
                        }
                        start = nextIndex;
                        end = start;
                        index = nextIndex + 1;
                    }
                } else {
                    for (int i = start; i < index; i++) {
                        out.append("0");
                    }
                    while ((index < length) &&
                            (freq.charAt(index) == '0' || freq.charAt(index) == '1')) {
                        index++;
                        out.append("0");
                    }
                    start = index;
                    end = start;
                    index++;
                }
                totalDis = 1;
                kmerNum = 1;
                lastIndex = -1;
                nextIndex = -1;
            }

            if (find == true) {
                if (dis <= combineLength) {
                    end = index;
                    totalDis = totalDis + dis;
                    kmerNum++;
                    if (kmerNum == 2) {
                        nextIndex = index - 1;
                    }
                    if (totalDis >= lengthBound && kmerNum / totalDis >= denseBound) {
                        lastIndex = index;
                    }
                } else {
                    if (totalDis >= lengthBound) {
                        if (lastIndex != -1) {
                            for (int i = start; i < lastIndex; i++) {
                                out.append("1");
                            }
                            index = lastIndex;
                            while ((index < length) &&
                                    (freq.charAt(index) == '0' || freq.charAt(index) == '1')) {
                                index++;
                                out.append("0");
                            }
                            start = index;
                            end = start;
                            index++;
                        } else {
                            for (int i = start; i < nextIndex; i++) {
                                out.append("0");
                            }
                            start = nextIndex;
                            end = start;
                            index = nextIndex + 1;
                        }
                    } else {
                        for (int i = start; i < index - 1; i++) {
                            out.append("0");
                        }
                        start = index - 1;
                        end = start;
                    }
                    totalDis = 1;
                    kmerNum = 1;
                    lastIndex = -1;
                    nextIndex = -1;
                }
            }

            if (index >= length) {
                if (totalDis >= lengthBound) {
                    if (lastIndex > 0) {
                        for (int i = start; i < lastIndex; i++) {
                            out.append("1");
                        }
                        index = lastIndex;
                        while ((index < length) &&
                                (freq.charAt(index) == '0' || freq.charAt(index) == '1')) {
                            index++;
                            out.append("0");
                        }
                        start = index;
                        end = start;
                        index++;
                    } else {
                        for (int i = start; i < nextIndex; i++) {
                            out.append("0");
                        }
                        start = nextIndex;
                        end = start;
                        index = nextIndex + 1;
                    }
                    totalDis = 1;
                    kmerNum = 1;
                    lastIndex = -1;
                    nextIndex = -1;
                } else {
                    for (int i = start; i < length; i++) {
                        out.append("0");
                    }
                }
            }
        }

        return out.toString();
    }
}
