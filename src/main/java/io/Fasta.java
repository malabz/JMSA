package io;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
//import java.util.regex.Pattern;

/**
 * read or write the fasta file
 */
public class Fasta {
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
                if (temp.length() > 0 && temp.charAt(0) == '>') {
                    if (line.length() != 0) {
                        strs.add(line.toString().toLowerCase());
                        line = new StringBuilder();
                    }
                    labels.add(temp.toLowerCase());
                }
                else if (temp.length() > 0) {
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
                int j = 1;
                for (; j * 60 < strings[i].length(); j++) {
                    bw.write(strings[i].substring(j*60-60, j*60)+"\n");
                }
                bw.write(strings[i].substring(j*60-60)+"\n");
            }
        }
    }

//    public String[][] NexToFasta(String path) throws IOException{
//        List<String> strs = new ArrayList<>();
//        List<String> labels = new ArrayList<>();
//        boolean state = false;
//        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
//            String temp;
//            while ((temp = br.readLine()) != null) {
//                if (temp.startsWith(";") && state) break;
//                if (state) {
//                    String[] labelStr = temp.split("\\s+", 3);
//                    labels.add(">" + labelStr[0]);
//                    strs.add(labelStr[1]);
//                }
//                if (temp.startsWith("matrix")) state = true;
//            }
//        }
//        return new String[][]{labels.toArray(new String[0]), strs.toArray(new String[0])};
//    }
//    public void deleteGap(String inpath, String outpath) throws IOException {
//        String[][] res = readFasta(inpath);
//        Pattern p = Pattern.compile("-");
//        for (int i = 0; i < res[1].length; i++) res[1][i] = p.matcher(res[1][i]).replaceAll("");
//        writeFasta(res[1], res[0], outpath);
//    }
}