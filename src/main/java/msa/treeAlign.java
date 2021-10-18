package msa;

import java.util.HashMap;

import hierCluster.guidetree;
import io.str;
import measure.*;
import psa.FastMultiAlign;
import psa.dsa;
import psa.multiDP;
import psa.multiKband;

public class treeAlign {
    private String[] straligned;
    private final String Treemode;
    private int[] orders;
    private final int num;
    private final int kk;


    /**
     * Choose one Gen tree mode "nj" or "upgma"
     */
    public treeAlign(String[] strs, String treemode) {
        this.num = strs.length;
        this.Treemode = treemode;
        kk = score.getK(strs, false);
        System.out.println("k:" + kk);
        Align(strs);
        reOrder();
    }

    /**
     * no output (silent)
     */
    public treeAlign(String[] strs, String treemode, int silent) {
        this.num = strs.length;
        this.Treemode = treemode;
        kk = score.getK(strs, false);
        AlignSlient(strs);
        reOrder();
    }

    /**
     * To get the alignment results.
     */
    public String[] getStrsAlign() { return this.straligned; }

    private int[] combineLabels(HashMap<Integer, int[]> labelsList, int l1, int l2) {
        int[] listL1 = l1 < num ? new int[]{l1} : labelsList.get(l1);
        int[] listL2 = l2 < num ? new int[]{l2} : labelsList.get(l2);
        int[] res = new int[listL1.length + listL2.length];
        System.arraycopy(listL1, 0, res, 0, listL1.length);
        System.arraycopy(listL2, 0, res, listL1.length, listL2.length);
        return res;
    }

    private void reOrder() {
        int i = 0;
        String[] temp = new String[num];
        for (int j : orders) {
            temp[j] = straligned[i++];
        }
        straligned = temp;
    }

    private void Align(String[] strs) {
        
        System.out.println();
        System.out.println("build the " + Treemode + " Tree");
        System.out.println();

        long startTime = System.currentTimeMillis();
        guidetree gTree = new guidetree(strs, this.Treemode);
        int[][] treeList = gTree.genTreeList(0);
        HashMap<Integer, String[]> strsList = new HashMap<>();
        HashMap<Integer, int[]> labelsList = new HashMap<>();
        HashMap<Integer, Integer> lengthList = new HashMap<>();
        long endTime = System.currentTimeMillis();

        System.out.println("time: "+((endTime-startTime)/1000)+"s\n");
        System.out.println("Align the seqs along the Tree\n");
        
        startTime = System.currentTimeMillis();
        kmer km = new kmer(strs, 4);
        char[] alphabet = km.Counter();
        int len = treeList.length, i = 0;
        for (int[] readyAlign : treeList) {
            String outToScreen = "    " + (i + 1) + " / " + len;
            System.out.print(outToScreen);
            String[] strsA, strsB;
            if (readyAlign[0] < this.num && readyAlign[1] < this.num) {
                dsa pa = new dsa(strs[readyAlign[0]], strs[readyAlign[1]], "fmindex");
                String[] strsC = pa.getStrAlign();
                
                strsList.put(readyAlign[2], strsC);
                labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
                lengthList.put(readyAlign[2], strs[readyAlign[0]].length() >= strs[readyAlign[1]].length() ? 0 : 1);

                i++;
                System.out.print(str.repeat("\b", outToScreen.length()));
                continue;
            }
            int idxcA = 0, idxcB = 0;
            if (readyAlign[0] < this.num) {
                strsA = new String[1];
                strsA[0] = strs[readyAlign[0]];
            }
            else {
                strsA = strsList.remove(readyAlign[0]);
                idxcA = lengthList.remove(readyAlign[0]);
            }
            if (readyAlign[1] < this.num) {
                strsB = new String[1];
                strsB[0] = strs[readyAlign[1]];
            }
            else {
                strsB = strsList.remove(readyAlign[1]);
                idxcB = lengthList.remove(readyAlign[1]);
            }
            int ln = Math.abs(strsB[0].length() - strsA[0].length());
            double diff = (double) ln / Math.max(strsB[0].length(), strsA[0].length());
            String[] strsC;
            if (diff > 0.4) { 
                multiDP mkband = new multiDP(strsA, strsB, alphabet);
                strsC = mkband.getStrsAlign();
            }
            else if (diff > 0.1) {
                multiKband mkband = new multiKband(strsA, strsB, alphabet, kk);
                strsC = mkband.getStrsAlign();
            }
            else {
                FastMultiAlign mkband = new FastMultiAlign(strsA, strsB, alphabet, kk, idxcA, idxcB);
                strsC = mkband.getStrsAlign();
            }
            strsList.put(readyAlign[2], strsC);
            
            labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
            labelsList.remove(readyAlign[0]);
            labelsList.remove(readyAlign[1]);

            if (strsA[idxcA].length() >= strsB[idxcB].length()) {
                lengthList.put(readyAlign[2], idxcA);
            }
            else {
                lengthList.put(readyAlign[2], idxcB + strsA.length);
            }

            i++;
            System.out.print(str.repeat("\b", outToScreen.length()));
        }
        endTime = System.currentTimeMillis();
        System.out.println("    " + len + " / " + len + "\n");
        System.out.println("time: "+((endTime-startTime)/1000)+"s\n");
        this.orders = labelsList.get(treeList[treeList.length-1][2]);
        this.straligned = strsList.get(treeList[treeList.length-1][2]);
        // System.out.println(" sps: " + String.format("%.3f", score.sps(straligned)));
        // System.out.println("  tc: " + String.format("%.3f", score.tc(straligned)));
        // System.out.println("  tc: " + String.format("%.3f", score.tcStrict(straligned)));
        
    }


    private void AlignSlient(String[] strs) {
        guidetree gTree = new guidetree(strs, this.Treemode);
        int[][] treeList = gTree.genTreeList(1);
        HashMap<Integer, String[]> strsList = new HashMap<>();
        HashMap<Integer, int[]> labelsList = new HashMap<>();
        HashMap<Integer, Integer> lengthList = new HashMap<>();
        kmer km = new kmer(strs, 4);
        char[] alphabet = km.Counter();
        for (int[] readyAlign : treeList) {
            String[] strsA, strsB;
            if (readyAlign[0] < this.num && readyAlign[1] < this.num) {
                dsa pa = new dsa(strs[readyAlign[0]], strs[readyAlign[1]], "fmindex");
                String[] strsC = pa.getStrAlign();
                strsList.put(readyAlign[2], strsC);
                labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
                lengthList.put(readyAlign[2], strs[readyAlign[0]].length() >= strs[readyAlign[1]].length() ? 0 : 1);
                continue;
            }
            int idxcA = 0, idxcB = 0;
            if (readyAlign[0] < this.num) {
                strsA = new String[1];
                strsA[0] = strs[readyAlign[0]];
            }
            else {
                strsA = strsList.remove(readyAlign[0]);
                idxcA = lengthList.remove(readyAlign[0]);
            }
            if (readyAlign[1] < this.num) {
                strsB = new String[1];
                strsB[0] = strs[readyAlign[1]];
            }
            else {
                strsB = strsList.remove(readyAlign[1]);
                idxcB = lengthList.remove(readyAlign[1]);
            }
            int ln = Math.abs(strsB[0].length() - strsA[0].length());
            double diff = (double) ln / Math.max(strsB[0].length(), strsA[0].length());
            String[] strsC;
            if (diff > 0.4) { 
                multiDP mkband = new multiDP(strsA, strsB, alphabet);
                strsC = mkband.getStrsAlign();
            }
            else if (diff > 0.1) {
                multiKband mkband = new multiKband(strsA, strsB, alphabet, kk);
                strsC = mkband.getStrsAlign();
            }
            else {
                FastMultiAlign mkband = new FastMultiAlign(strsA, strsB, alphabet, kk, idxcA, idxcB);
                strsC = mkband.getStrsAlign();
            }
            strsList.put(readyAlign[2], strsC);
            labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
            labelsList.remove(readyAlign[0]);
            labelsList.remove(readyAlign[1]);
            if (strsA[idxcA].length() >= strsB[idxcB].length()) {
                lengthList.put(readyAlign[2], idxcA);
            }
            else {
                lengthList.put(readyAlign[2], idxcB+strsA.length);
            }
        }
        this.orders = labelsList.get(treeList[treeList.length-1][2]);
        this.straligned = strsList.get(treeList[treeList.length-1][2]);
    }
}