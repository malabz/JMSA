package msa;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import measure.score;
import psa.STAlign;

public class reAlign {
    private String[] strsaligned;

    public reAlign(String[] strsed) {
        strsaligned = strsed;
        rAlign();
    }
    
    /**
     * 最基本的重新比对的策略，简单省时
     */
    public void rAlign() {
        System.out.println();
        System.out.println("Realign the results");
        System.out.println();
        System.out.println("find the mis area");
        int[][] record = findMis();
        if (record.length < 1) return;
        String[][] subStrsed = new String[record.length][strsaligned.length];
        int idxi = 0;
        System.out.println();
        System.out.println("get the seqs and align");
        System.out.println();
        for (int[] r : record) {
            String scren = "    " + (idxi + 1) + "/" + record.length;
            System.out.print(scren);
            String[] rdStrs = new String[strsaligned.length];
            for (int i = 0; i < rdStrs.length; i++) {
                rdStrs[i] = strsaligned[i].substring(r[0], r[1] + 1).replaceAll("-", "");
            }
            centerAlign cAlign = new centerAlign(rdStrs, 1);
            subStrsed[idxi++] = cAlign.getStrsAlign();
            System.out.print("\b".repeat(scren.length()));
        }
        System.out.println();
        System.out.println("\ncombine the results");
        System.out.println();
        for (int i = 0; i < strsaligned.length; i++) {
            String scren = "    " + (i + 1) + "/" + strsaligned.length;
            System.out.print(scren);
            idxi = 0;
            int idxj;
            StringBuilder strsed = new StringBuilder();
            for (int j = 0; j < subStrsed.length; j++) {
                idxj = record[j][0];
                strsed.append(strsaligned[i].subSequence(idxi, idxj));
                strsed.append(subStrsed[j][i]);
                idxi = record[j][1] + 1;
            }
            if (idxi < strsaligned[i].length()) {
                strsed.append(strsaligned[i].substring(idxi));
            }
            strsaligned[i] = strsed.toString();
            System.out.print("\b".repeat(scren.length()));
        }
        System.out.println("\n");
        score sc = new score();
        System.out.println(" sps: " + String.format("%.3f", sc.sps(strsaligned)));
        System.out.println("  tc: " + String.format("%.3f", sc.tc(strsaligned)));
    }

    public void rrAlign () {
        System.out.println();
        System.out.println("Realign the results");
        System.out.println();
        System.out.println("find the mis area");
        int[][] record = findMis();
        if (record.length < 1) return;
        String[][] subStrsed = new String[record.length][strsaligned.length];
        int idxi = 0;
        System.out.println();
        System.out.println("get the seqs and align");
        System.out.println();
        for (int[] r : record) {
            String scren = "    " + (idxi + 1) + "/" + record.length;
            System.out.print(scren);
            String[] rdStrs = new String[strsaligned.length];
            for (int i = 0; i < rdStrs.length; i++) {
                rdStrs[i] = strsaligned[i].substring(r[0], r[1] + 1).replaceAll("-", "");
            }
            centerAlign cAlign = new centerAlign(rdStrs, 1);
            subStrsed[idxi++] = cAlign.getStrsAlign();
            System.out.print("\b".repeat(scren.length()));
        }
        System.out.println();
        System.out.println("combine the results");
        System.out.println();
        for (int i = 0; i < strsaligned.length; i++) {
            String scren = "    " + (i + 1) + "/" + strsaligned.length;
            System.out.print(scren);
            idxi = 0;
            int idxj;
            StringBuilder strsed = new StringBuilder();
            for (int j = 0; j < subStrsed.length; j++) {
                idxj = record[j][0];
                strsed.append(strsaligned[i].subSequence(idxi, idxj));
                strsed.append(subStrsed[j][i]);
                idxi = record[j][1] + 1;
            }
            if (idxi < strsaligned[i].length()) {
                strsed.append(strsaligned[i].substring(idxi));
            }
            strsaligned[i] = strsed.toString();
            System.out.print("\b".repeat(scren.length()));
        }
        System.out.println();
        score sc = new score();
        System.out.println(" sps: " + String.format("%.3f", sc.sps(strsaligned)));
        System.out.println("  tc: " + String.format("%.3f", sc.tc(strsaligned)));
    }

    public void rreAlign(String[] strs) {
        int[][] record = findMis();
        // kmer km = new kmer(this.strs);
        // char[] alphabet = km.Counter();
        int round = 0;
        // int dxiA = -1, dxiB = -1;
        score sc = new score();
        double s1 = sc.tc(strsaligned);
        while ( record.length > 1 && round < record.length) {
            for (int[] r : record) {
                System.out.println(r[0] + " " + r[1]);
                List<String> Part1 = new ArrayList<>();
                List<String> Part2 = new ArrayList<>();
                for (int i = 0; i < strsaligned.length; i++) {
                    String temp = strsaligned[i].substring(r[0], r[1] + 1).replaceAll("-", "");
                    if (temp.equals("")) { 
                        Part1.add(strsaligned[i]);
                        // if (strs[i].length() > lenA) dxiA = i;
                    }
                    else {
                        Part2.add(strs[i]);
                        // if (strs[i].length() > lenB) dxiB = i;
                    }
                }
                String[] strPart1 = Part1.toArray(String[]::new);
                String[] strPart2 = Part2.toArray(String[]::new);

                delGapStrings(strPart1);
                // centerAlign cAlign1 = new centerAlign(strPart1);
                // strPart1 = cAlign1.getStrsAlign();
                centerAlign cAlign2 = new centerAlign(strPart2);
                strPart2 = cAlign2.getStrsAlign();
                String Astr = PickRealOne(strPart1);
                String Bstr = PickRealOne(strPart2);

                STAlign sAlign = new STAlign(Astr, Bstr);
                String[] ABstr = sAlign.getStrAlign();

                int[] markA = new int[Astr.length() + 1];
                int[] markB = new int[Bstr.length() + 1];

                int j = 0, counter = 0;
                for (char c : ABstr[0].toCharArray()) {
                    if (c == '-') counter++;
                    else { markA[j++] = counter; counter = 0; }
                }
                markA[j] = counter;
                j = 0; counter = 0;
                for (char c : ABstr[1].toCharArray()) { 
                    if (c == '-') counter++;
                    else { markB[j++] = counter; counter = 0; }
                }

                String[] nnnstrs = new String[strPart1.length + strPart2.length];
                for (int i = 0; i < strPart1.length; i++) {
                    nnnstrs[i] = insertGap(markA, strPart1[i]);
                }
                for (int i = 0; i < strPart2.length; i++) {
                    nnnstrs[i+strPart1.length] = insertGap(markB, strPart2[i]);
                }

                // multiKband mKband = new multiKband(strPart1, strPart2, alphabet);
                double s2 = sc.tc(nnnstrs);
                // double s2 = sc.tc(mKband.getStrsAlign());
                float l1 = strsaligned[0].length();
                // float l2 = mKband.getStrsAlign()[0].length();
                float l2 = nnnstrs[0].length();
                System.out.println(s1+":"+s2);
                round++;
                if (s2 > s1 || (s1 - s2) <= ((l1 - l2)/l2)) {
                    // strsaligned = mKband.getStrsAlign();
                    strsaligned = nnnstrs;
                    break;
                }
            }
            record = findMis();
            round++;
        }
    }

    private String insertGap(int[] mark, String seq) {
        StringBuilder seqGap = new StringBuilder();
        int len = mark.length;
        for (int i = 0; i < len; i++) {
            seqGap.append("-".repeat(mark[i]));
            if (i < len - 1) seqGap.append(seq.charAt(i));
        }
        return seqGap.toString();
    }

    private String PickRealOne (String[] strings) {
        char[] res = new char[strings[0].length()];
        for (int i = 0; i < res.length; i++) {
            HashMap<Character, Integer> charCounter = new HashMap<>();
            for (String string : strings) {
                char tempc = string.charAt(i);
                charCounter.put(tempc, charCounter.containsKey(tempc) ? charCounter.get(tempc) + 1 : 1);
            }
            charCounter.remove('-');
            char resc = 'n'; int num = 0;
            for (char c : charCounter.keySet()) {
                if (charCounter.get(c) > num) {
                    num = charCounter.get(c);
                    resc = c;
                }
            }
            res[i] = resc;
        }
        return new String(res);
    }


    private void delGapStrings (String[] strings) {
        int len = strings[0].length();
        int length = strings.length;
        int[] states = new int[len];
        for (int i = 0; i < len; i++) {
            int state = 0;
            for (String string : strings) {
                if (string.charAt(i) != '-') {
                    state = 1;
                    break;
                }
            }
            states[i] = state;
        }
        for (int i = 0; i < length; i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < len; j++) {
                if (states[j] == 1) sb.append(strings[i].charAt(j));
            }
            strings[i] = sb.toString();
        }
    }

    private int[][] findMis () {
        List<int[]> record = new ArrayList<>();
        int start = -1, end = -1;
        for (int i = 0; i < strsaligned[0].length(); i++) {
            int numgap = 0;
            for (String s : strsaligned) {
                if (s.charAt(i) == '-') numgap++;
            }
            if ((float) numgap/strsaligned.length > 0.6) { 
                start = start == -1 ? i : start; 
                end = i;
            }
            else { 
                if (start != end) {
                    record.add(new int[]{start, end});
                }
                start = -1; 
                end = -1; 
            }
        }
        if (start != end) {
            record.add(new int[]{start, end});
        }
        return record.toArray(int[][]::new);
    }
}
