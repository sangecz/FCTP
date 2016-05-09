package cz.sange;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by sange on 05/05/16.
 */
public class LogWriter {

    private String filePath;
    private BufferedWriter bw;
    private FileWriter fw;

    public LogWriter(String filePath) {
        this.filePath = filePath;

        open();
    }

    private void open() {
        File file = new File(filePath);

        fw = null;
        try {

            // if file doesnt exists, then create it
            if (!file.exists()) {
                file.createNewFile();
            }

            fw = new FileWriter(file.getAbsoluteFile(), true);
            bw = new BufferedWriter(fw);

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void write (String output) {
        try {

            bw.write(output + "\n");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void close () {
        try {

            if (bw != null) {
                bw.close();
            }

            if(fw != null) {
                fw.close();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
