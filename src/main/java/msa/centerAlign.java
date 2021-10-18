package msa;

import io.str;
import psa.dsa;
// import measure.score;
import psa.FMAlign;
import psa.STAlign;

public class centerAlign {
    private final String[] strsaligned;
    private final String mode;
    private final int num;
    private int maxcol;


    /**
     *
     * @param strs
     * @param mode "kband" "suffix" "fmindex"
     */
    public centerAlign(String[] strs, String mode) {
        this.num = strs.length;
        this.mode = mode;
        this.strsaligned = Align(strs);
    }

    /**
     * mode "fmindex"
     * @param strs
     */
    public centerAlign(String[] strs) {
        this.num = strs.length;
        this.mode = "fmindex";
        this.strsaligned = Align(strs);
    }

    /**
     * no output in silence mode "fmindex"
     * @param strs
     * @param silent
     */
    public centerAlign(String[] strs, int silent) {
        this.num = strs.length;
        this.mode = "fmindex";
        this.strsaligned = AlignSilent(strs);
    }

    public String[] getStrsAlign() { return this.strsaligned; }

    public int getLongestRow() { return maxcol; }

    private String[] Align(String[] strs) {

        long startTime, endTime;
        startTime = System.currentTimeMillis();
        System.out.println("\nfind the center seq");

        int maxcol = findMaxcol(strs);
        System.out.println("\n    No."+maxcol);
        System.out.println("\nAlign the center seq\n");

        String[][] resultsAlign = new String[num - 1][2];
        switch (this.mode) {
            case "suffix":
                STAlign stAlign = new STAlign(strs[maxcol]);
                for (int i = 0; i < num - 1; i++) {
                    String outToScreen = "    " + (i + 1) + " / " + (num - 1);
                    System.out.print(outToScreen);
                    int j = i >= maxcol ? i + 1 : i;
                    stAlign.AlignStrB(strs[j]);
                    resultsAlign[i] = stAlign.getStrAlign();
                    System.out.print(str.repeat("\b", outToScreen.length()));
                }
                System.out.println("    " + (num - 1) + " / " + (num - 1));
                break;
            case "kband":
                for (int i = 0; i < num - 1; i++) {
                    String outToScreen = "    " + (i + 1) + " / " + (num - 1);
                    System.out.print(outToScreen);
                    int j = i >= maxcol ? i + 1 : i;
                    dsa dalign = new dsa(strs[maxcol], strs[j], "kband");
                    resultsAlign[i] = dalign.getStrAlign();
                    System.out.print(str.repeat("\b", outToScreen.length()));
                }
                System.out.println("    " + (num - 1) + " / " + (num - 1));
                break;
            case "fmindex":
                FMAlign fmAlign = new FMAlign(strs[maxcol]);
                for (int i = 0; i < num - 1; i++) {
                    String outToScreen = "    " + (i + 1) + " / " + (num - 1);
                    System.out.print(outToScreen);
                    int j = i >= maxcol ? i + 1 : i;
                    fmAlign.AlignStrB(strs[j]);
                    resultsAlign[i] = fmAlign.getStrAlign();
                    System.out.print(str.repeat("\b", outToScreen.length()));
                }
                System.out.println("    " + (num - 1) + " / " + (num - 1));
                break;
            default:
                throw new IllegalArgumentException("unkown mode: " + mode);
        }

        String centerSeq = strs[maxcol];

        // mark the gaps in center seq
        int[] markInsertion = new int[centerSeq.length() + 1];

        for (String[] str2 : resultsAlign) {
            int i = 0, counter = 0;
            for (char c : str2[0].toCharArray()) {
                if (c == '-') counter++;
                else {
                    markInsertion[i] = Math.max(markInsertion[i], counter);
                    counter = 0;
                    i++;
                }
            }
            markInsertion[i] = Math.max(markInsertion[i], counter);
        }

        // insert the gap
        System.out.println("\nAlign all seqs\n");
        String[] newStrsaligned = new String[num];
        int idxinsert = 0;
        newStrsaligned[maxcol] = insertGap(markInsertion, centerSeq);

        for (String[] str2 : resultsAlign) {
            String outToScreen = "    " + (idxinsert + 1) + " / " + (num - 1);
            System.out.print(outToScreen);
            char[] tempA = str2[0].toCharArray();
            int[] mark = new int[tempA.length + 1];
            int pi = 0, pj = 0, total = 0;
            for (char c : tempA) {
                if (c == '-') total++;
                else {
                    mark[pi++] = markInsertion[pj++] - total;
                    while (total != 0) { pi++; total--; }
                }
            }
            mark[pi] = markInsertion[pj] - total;
            if (idxinsert >= maxcol) {
                newStrsaligned[++idxinsert] = insertGap(mark, str2[1]);
            }
            else {
                newStrsaligned[idxinsert++] = insertGap(mark, str2[1]);
            }
            System.out.print(str.repeat("\b", outToScreen.length()));
        }

        endTime = System.currentTimeMillis();
        System.out.println("    " + (num - 1) + " / " + (num - 1) + "\n");
        System.out.println("time: "+((endTime-startTime)/1000)+"s\n");
        // System.out.println(" sps: " + String.format("%.3f", score.sps(newStrsaligned)));
        // System.out.println("  tc: " + String.format("%.3f", score.tc(newStrsaligned)));
        // System.out.println("  tc: " + String.format("%.3f", score.tcStrict(newStrsaligned)));
        return newStrsaligned;
    }

    private String[] AlignSilent(String[] strs) {

        maxcol = findMaxcol(strs);
        String[][] resultsAlign = new String[num - 1][2];
        switch (this.mode) {
            case "suffix" :
                STAlign stAlign = new STAlign(strs[maxcol]);
                for (int i = 0; i < num - 1; i++) {
                    int j = i >= maxcol ? i + 1 : i;
                    stAlign.AlignStrB(strs[j]);
                    resultsAlign[i] = stAlign.getStrAlign();
                }
                break;
            case "kband" :
                for (int i = 0; i < num - 1; i++) {
                    int j = i >= maxcol ? i + 1 : i;
                    dsa dalign = new dsa(strs[maxcol], strs[j], "kband");
                    resultsAlign[i] = dalign.getStrAlign();
                }
                break;
            case "fmindex" :
                FMAlign fmAlign = new FMAlign(strs[maxcol]);
                for (int i = 0; i < num - 1; i++) {
                    int j = i >= maxcol ? i + 1 : i;
                    fmAlign.AlignStrB(strs[j]);
                    resultsAlign[i] = fmAlign.getStrAlign();
                }
                break;
            default :
                System.out.println("(centeralign) mode is wrong");
                System.exit(0);
        }
        String centerSeq = strs[maxcol];

        // mark the gaps in center seq
        int[] markInsertion = new int[centerSeq.length() + 1];
        for (String[] str2 : resultsAlign) {
            int i = 0;
            int counter = 0;
            char[] temp = str2[0].toCharArray();
            for (char c : temp) {
                if (c == '-') counter++;
                else {
                    markInsertion[i] = Math.max(markInsertion[i], counter);
                    counter = 0;
                    i++;
                }
            }
            markInsertion[i] = Math.max(markInsertion[i], counter);
        }

        // insert the gap
        String[] newstrsaligned = new String[num];
        int idxinsert = 0;
        newstrsaligned[maxcol] = insertGap(markInsertion, centerSeq);
        for (String[] str2 : resultsAlign) {
            char[] tempA = str2[0].toCharArray();
            int[] mark = new int[tempA.length + 1];
            int pi = 0, pj = 0;
            int total = 0;
            for (char c : tempA) {
                if (c == '-') total++;
                else {
                    mark[pi++] = markInsertion[pj++] - total;
                    while (total != 0) { pi++; total--; }
                }
            }
            mark[pi] = markInsertion[pj] - total;
            if (idxinsert >= maxcol) {
                newstrsaligned[++idxinsert] = insertGap(mark, str2[1]);
            }
            else {
                newstrsaligned[idxinsert++] = insertGap(mark, str2[1]);
            }
        }
        return newstrsaligned;
    }

    private String insertGap(int[] mark, String seq) {
        StringBuilder seqGap = new StringBuilder();
        int len = mark.length;
        for (int i = 0; i < len; i++) {
            seqGap.append(str.repeat("-", mark[i]));
            if (i < len - 1) seqGap.append(seq.charAt(i));
        }
        return seqGap.toString();
    }

    private int findMaxcol(String[] strs) {
        int max = 0;
        for (int i = 0; i < num; i++) {
            if (strs[i].length() > strs[max].length()) max = i;
        }
        return max;
    }
}
