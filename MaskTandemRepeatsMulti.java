import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

public class MaskTandemRepeatsMulti {

    private static AtomicInteger TRLength = new AtomicInteger(0);

    public static void main(String[] args) {
        if(args.length < 3){
            System.out.println("Usage: java MaskTandemRepeatsMultiThreaded genome.fa trf.info output.fa");
            return;
        }

        String genomeFile = args[0];
        String trfFile = args[1];
        String outputFile = args[2];

        int numThreads = Runtime.getRuntime().availableProcessors(); // 可自定义线程数
        long startTime = System.currentTimeMillis();

        try {
            // 1. 读取 TRF 文件到 Map<染色体编号, List<区间>>
            Map<String, List<int[]>> trMap = new HashMap<>();
            try(Scanner trScanner = new Scanner(new File(trfFile))){
                String seqId = null;
                while(trScanner.hasNextLine()){
                    String line = trScanner.nextLine().trim();
                    if(line.isEmpty()) continue;

                    if(line.charAt(0) == '@'){
                        seqId = line.substring(1); // 去掉 @
                        trMap.putIfAbsent(seqId, new ArrayList<>());
                    } else {
                        String[] parts = line.split("\\s+");
                        if(parts.length >= 2){
                            int start = Integer.parseInt(parts[0]);
                            int end = Integer.parseInt(parts[1]);
                            trMap.get(seqId).add(new int[]{start, end});
                        }
                    }
                }
            }

            // 2. 线程池
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);
            List<Future<String>> futures = new ArrayList<>();

            // 3. 逐条染色体读取 genome 文件，每条染色体两行
            try(Scanner genomeScanner = new Scanner(new File(genomeFile))){
                while(genomeScanner.hasNextLine()){
                    String header = genomeScanner.nextLine();
                    if(!header.startsWith(">")) continue;

                    if(!genomeScanner.hasNextLine()) break;
                    String sequenceLine = genomeScanner.nextLine();
                    char[] bases = sequenceLine.toCharArray();

                    // 对应 TRF
                    String seqId = header.substring(1);
                    List<int[]> intervals = trMap.getOrDefault(seqId, Collections.emptyList());

                    Future<String> f = pool.submit(() -> processChromosome(header, bases, intervals));
                    futures.add(f);
                }
            }

            pool.shutdown();

            // 4. 按序写入输出文件
            try(BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))){
                for(Future<String> f : futures){
                    writer.write(f.get());
                    writer.newLine();
                }
            }

            long endTime = System.currentTimeMillis();
            System.out.println("Tandem Repeat Finder Detected: " + TRLength.get() + " bp tandem repeats in the input genome.");
            System.out.println("runTime: " + (endTime - startTime) + " ms");

        } catch(Exception e){
            e.printStackTrace();
        }
    }

    private static String processChromosome(String header, char[] bases, List<int[]> intervals){
        for(int[] interval : intervals){
            int start = interval[0];
            int end = interval[1];
            for(int i = start; i < end && i < bases.length; i++){
                if(bases[i] != 'N'){
                    bases[i] = 'N';
                    TRLength.incrementAndGet();
                }
            }
        }
        return header + "\n" + new String(bases);
    }
}
