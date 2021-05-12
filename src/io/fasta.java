package io;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * read or write the fasta file
 */
public class fasta {
    /**
     * 
     * @param path
     * @return String[lables[], Strings[]]
     * @throws IOException
     */
    public String[][] readFasta(String path) throws IOException {
        List<String> strs = new ArrayList<>();
        List<String> labels = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {     
            String temp;
            StringBuilder line = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                if (temp.charAt(0) == '>') {
                    if (line.length() != 0) {
                        strs.add(line.toString().toLowerCase());
                        line = new StringBuilder();
                    }
                    labels.add(temp.toLowerCase());
                }
                else {
                    line.append(temp);
                }
            }
            strs.add(line.toString().toLowerCase());
        }
        return new String[][]{labels.toArray(new String[0]), strs.toArray(new String[0])};
    }

    /**
     * 
     * @param strings
     * @param labels
     * @param path
     */
    public void writeFasta(String[] strings, String[] labels, String path) throws IOException{
        try (Writer write = new FileWriter(path); BufferedWriter bw = new BufferedWriter(write)) {
            for (int i = 0; i < strings.length; i++) {
                bw.write(labels[i]+"\n");
                bw.write(strings[i].toUpperCase()+"\n");
            }
        }
    }

    public void deleteGap(String inpath, String outpath) throws IOException {
        String[][] res = readFasta(inpath);
        for (int i = 0; i < res[1].length; i++) res[1][i] = res[1][i].replaceAll("-", "");
        writeFasta(res[1], res[0], outpath);
    }
}