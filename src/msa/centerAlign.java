package msa;
import psa.dsa;
import psa.STAlign;
import measure.score;

public class centerAlign {
    private String[] strs;
    private String[] strsaligned;
    private String mode;
    private int num;

    /**
     * 
     * @param strs
     * @param labels
     * @param mode "kband" "suffix"
     */
    public centerAlign(String[] strs, String[] labels, String mode) {
        this.strs = strs;
        this.num = strs.length;
        this.mode = mode;
        this.strsaligned = Align();
    }

    /**
     * 
     * @param strs
     * @param mode "kband" "suffix"
     */
    public centerAlign(String[] strs, String mode) {
        this.strs = strs;
        this.num = strs.length;
        this.mode = mode;
        this.strsaligned = Align();
    }

    /**
     * 
     * @param strs
     */
    public centerAlign(String[] strs) {
        this.strs = strs;
        this.num = strs.length;
        this.mode = "suffix";
        this.strsaligned = Align();
    }

    public String[] getStrsAlign() {
        return this.strsaligned;
    }

    private String[] Align() {
        
        long startTime, endTime;
        startTime = System.currentTimeMillis();
        System.out.println("\nfind the center seq");
        
        int maxcol = findMaxcol();
        System.out.println("\n    No."+maxcol);
        System.out.println("\nAlign the center seq\n");

        String[][] resultsAlign = new String[num - 1][2];
        if (this.mode.equals("suffix")) {    
            STAlign CTalign = new STAlign(strs[maxcol]);
            for (int i = 0; i < num - 1; i++) {
                String outToScreen = "    " + (i + 1) + " / " + (num - 1);
                System.out.print(outToScreen);
                int j = i >= maxcol ? i + 1 : i;
                CTalign.AlignStrB(strs[j]);
                resultsAlign[i] = CTalign.getStrAlign();
                System.out.print("\b".repeat(outToScreen.length()));
            }
            System.out.println("    " + (num - 1) + " / " + (num - 1));
        }
        else if (this.mode.equals("kband")) {
            for (int i = 0; i < num - 1; i++) {
                String outToScreen = "    " + (i + 1) + " / " + (num - 1);
                System.out.print(outToScreen);
                int j = i >= maxcol ? i + 1 : i;
                dsa dalign = new dsa(strs[maxcol], strs[j], "kband");
                resultsAlign[i] = dalign.getStrAlign();
                System.out.print("\b".repeat(outToScreen.length()));
            }
            System.out.println("    " + (num - 1) + " / " + (num - 1));
        }
        
        String centerSeq = this.strs[maxcol];
        this.strs = null;

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
        System.out.println("\nAlign all seqs\n");
        String[] newStrsaligned = new String[num];
        int idxinsert = 0;
        newStrsaligned[maxcol] = insertGap(markInsertion, centerSeq);

        for (String[] str2 : resultsAlign) {
            String outToScreen = "    " + (idxinsert + 1) + " / " + (num - 1);
            System.out.print(outToScreen);
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
                newStrsaligned[++idxinsert] = insertGap(mark, str2[1]);
            }
            else {
                newStrsaligned[idxinsert++] = insertGap(mark, str2[1]);
            }
            System.out.print("\b".repeat(outToScreen.length()));
        }
        
        resultsAlign = null;
        endTime = System.currentTimeMillis();
        System.out.println("    " + (num - 1) + " / " + (num - 1) + "\n");
        System.out.println("\ntime: "+((endTime-startTime)/1000)+"s\n");
        score sc = new score();
        System.out.println(" sps: " + String.format("%.3f", sc.sps(newStrsaligned)));
        System.out.println("  tc: " + String.format("%.3f", sc.tc(newStrsaligned)));
        return newStrsaligned;
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

    private int findMaxcol() {
        int max = 0;
        for (int i = 0; i < num; i++) {
            if (strs[i].length() > strs[max].length()) max = i;
        }
        return max;
    }
}
