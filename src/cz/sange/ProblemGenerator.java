package cz.sange;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;

/**
 * Created by sange on 01/05/16.
 */
public class ProblemGenerator {

    public static final int MIN = 1;
    public static final int MAX = 50;
    public static final int W = 4;
    public static final int C = 4;
    public static final String FILENAME = "src/cz/sange/04.in";

    private static Random random = new Random();

    public static void main(String[] args) {
        try {

            File file = new File(FILENAME);

            // if file doesnt exists, then create it
            if (!file.exists()) {
                file.createNewFile();
            }

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            // #W #C
            bw.write(W + " " + C + "\n");

            // supplies
            int [] supplies = new int[W];
            for (int i = 0; i < W; i++) {
                supplies[i] = getRandomInt(MIN, MAX);
                bw.write(supplies[i] + " ");
            }
            bw.write('\n');

            // demands
            for (int i = C - 1; i >= 0; i--) {
                bw.write(supplies[i] + " ");
            }
            bw.write('\n');

            // transport costs
            for (int i = 0; i < W; i++) {
                for (int j = 0; j < C; j++) {
                    bw.write(getRandomInt(MIN, MAX) + " ");
                }
                bw.write('\n');
            }


            // fixed costs
            for (int i = 0; i < W; i++) {
                for (int j = 0; j < C; j++) {
                    bw.write(getRandomInt(MIN, MAX) + " ");
                }
                bw.write('\n');
            }

            bw.close();



        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public static int getRandomInt(int min, int max){
        return random.nextInt(max - min + 1) + min;
    }
}
