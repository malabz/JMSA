package msa;


import hierCluster.guidetree;
import java.util.HashMap;
import psa.multiKband;
import measure.*;
import psa.dsa;

public class treeAlign {
    private String[] straligned;
    private final String Treemode;
    private int[] orders;
    private final int num;
    private int kk;


    /**
     * used for cenTree
     * @param strs
     * @param straligned
     */
    public treeAlign (String[] strs, String[] straligned, String treemode) {
        this.num = strs.length;
        this.Treemode = treemode;
        this.straligned = straligned;
        ReAlign(strs);
    }


    /**
     * Choose one Gen tree mode "nj" or "upgma"
     * @param strs
     */
    public treeAlign(String[] strs, String treemode) {
        this.num = strs.length;
        this.Treemode = treemode;
        score sc = new score();
        kk = sc.getK(strs, false);
        System.out.println("k:" + kk);
        Align(strs);
        reOrder();
    }

    /**
     * no output (silent)
     * @param strs
     * @param silent
     */
    public treeAlign(String[] strs, String treemode, int silent) {
        this.num = strs.length;
        this.Treemode = treemode;
        score sc = new score();
        kk = sc.getK(strs, false);
        AlignSlient(strs);
        reOrder();
    }

    /**
     * To get the alignment results.
     */
    public String[] getStrsAlign() {
        return this.straligned;
    }

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
        long endTime = System.currentTimeMillis();

        System.out.println("time: "+((endTime-startTime)/1000)+"s\n");
        System.out.println("Align the seqs along the Tree\n");
        
        startTime = System.currentTimeMillis();
        kmer km = new kmer(strs);
        char[] alphabet = km.Counter();
        int len = treeList.length, i = 0;
        for (int[] readyAlign : treeList) {
            String outToScreen = "    " + (i + 1) + " / " + len;
            System.out.print(outToScreen);
            String[] strsA, strsB;
            if (readyAlign[0] < this.num && readyAlign[1] < this.num) {
                System.currentTimeMillis();
                // TODO to use fmindex
                dsa pa = new dsa(strs[readyAlign[0]], strs[readyAlign[1]], "suffix");
                String[] strsC = pa.getStrAlign();
                
                strsList.put(readyAlign[2], strsC);
                strsList.remove(readyAlign[0]);
                strsList.remove(readyAlign[1]);
                
                labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
                labelsList.remove(readyAlign[0]);
                labelsList.remove(readyAlign[1]);
                System.currentTimeMillis();

                i++;
                System.out.print("\b".repeat(outToScreen.length()));
                continue;
            }
            if (readyAlign[0] < this.num) {
                strsA = new String[1];
                strsA[0] = strs[readyAlign[0]];
            }
            else {
                strsA = strsList.remove(readyAlign[0]);
            }
            if (readyAlign[1] < this.num) {
                strsB = new String[1];
                strsB[0] = strs[readyAlign[1]];
            }
            else {
                strsB = strsList.remove(readyAlign[1]);
            }

            multiKband mkband = new multiKband(strsA, strsB, alphabet, kk);
            String[] strsC = mkband.getStrsAlign();

            strsList.put(readyAlign[2], strsC);
            
            labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
            labelsList.remove(readyAlign[0]);
            labelsList.remove(readyAlign[1]);

            i++;
            System.out.print("\b".repeat(outToScreen.length()));
        }
        endTime = System.currentTimeMillis();
        System.out.println("    " + len + " / " + len + "\n");
        System.out.println("time: "+((endTime-startTime)/1000)+"s\n");
        this.orders = labelsList.get(treeList[treeList.length-1][2]);
        this.straligned = strsList.get(treeList[treeList.length-1][2]);
        score sc = new score();
        System.out.println(" sps: " + String.format("%.3f", sc.sps(straligned)));
        System.out.println("  tc: " + String.format("%.3f", sc.tc(straligned)));
        
    }


    /**
     * TODO
     * 这个函数是为了之后校正用的，对齐之后还是需要转换顺序的
     */
    private void ReAlign(String[] strs) {
        
        System.out.println();
        System.out.println("build the " + Treemode + " Tree");
        System.out.println();

        long startTime = System.currentTimeMillis();
        guidetree gTree = new guidetree(this.straligned, 2, Treemode);
        int[][] treeList = gTree.genTreeList(0);
        HashMap<Integer, String[]> strsList = new HashMap<>();
        HashMap<Integer, int[]> labelsList = new HashMap<>();
        long endTime = System.currentTimeMillis();

        System.out.println("time: "+((endTime-startTime)/1000)+"s\n");
        System.out.println("Align the seqs along the Tree\n");
        
        startTime = System.currentTimeMillis();
        kmer km = new kmer(strs);
        char[] alphabet = km.Counter();
        int len = treeList.length, i = 0;
        for (int[] readyAlign : treeList) {
            String outToScreen = "    " + (i + 1) + " / " + len;
            System.out.print(outToScreen);
            String[] strsA, strsB;
            if (readyAlign[0] < this.num && readyAlign[1] < this.num) {
                System.currentTimeMillis();
                dsa pa = new dsa(strs[readyAlign[0]], strs[readyAlign[1]], "suffix");
                String[] strsC = pa.getStrAlign();
                
                strsList.put(readyAlign[2], strsC);
                strsList.remove(readyAlign[0]);
                strsList.remove(readyAlign[1]);
                
                labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
                labelsList.remove(readyAlign[0]);
                labelsList.remove(readyAlign[1]);
                System.currentTimeMillis();

                i++;
                System.out.print("\b".repeat(outToScreen.length()));
                continue;
            }
            if (readyAlign[0] < this.num) {
                strsA = new String[1];
                strsA[0] = strs[readyAlign[0]];
            }
            else {
                strsA = strsList.remove(readyAlign[0]);
            }
            if (readyAlign[1] < this.num) {
                strsB = new String[1];
                strsB[0] = strs[readyAlign[1]];
            }
            else {
                strsB = strsList.remove(readyAlign[1]);
            }

            multiKband mkband = new multiKband(strsA, strsB, alphabet, kk);
            String[] strsC = mkband.getStrsAlign();

            strsList.put(readyAlign[2], strsC);
            
            labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
            labelsList.remove(readyAlign[0]);
            labelsList.remove(readyAlign[1]);

            i++;
            System.out.print("\b".repeat(outToScreen.length()));
        }
        endTime = System.currentTimeMillis();
        System.out.println("    " + len + " / " + len + "\n");
        System.out.println("\ntime: "+((endTime-startTime)/1000)+"s\n");
        this.orders = labelsList.get(treeList[treeList.length-1][2]);
        this.straligned = strsList.get(treeList[treeList.length-1][2]);
        score sc = new score();
        System.out.println(" sps: " + String.format("%.3f", sc.sps(straligned)));
        System.out.println("  tc: " + String.format("%.3f", sc.tc(straligned)));
    }

    private void AlignSlient(String[] strs) {
        guidetree gTree = new guidetree(strs, this.Treemode);
        int[][] treeList = gTree.genTreeList(1);
        HashMap<Integer, String[]> strsList = new HashMap<>();
        HashMap<Integer, int[]> labelsList = new HashMap<>();
        kmer km = new kmer(strs);
        char[] alphabet = km.Counter();
        for (int[] readyAlign : treeList) {
            String[] strsA, strsB;
            if (readyAlign[0] < this.num && readyAlign[1] < this.num) {
                System.currentTimeMillis();
                dsa pa = new dsa(strs[readyAlign[0]], strs[readyAlign[1]], "suffix");
                String[] strsC = pa.getStrAlign();
                strsList.put(readyAlign[2], strsC);
                strsList.remove(readyAlign[0]);
                strsList.remove(readyAlign[1]);
                labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
                labelsList.remove(readyAlign[0]);
                labelsList.remove(readyAlign[1]);
                System.currentTimeMillis();
                continue;
            }
            if (readyAlign[0] < this.num) {
                strsA = new String[1];
                strsA[0] = strs[readyAlign[0]];
            }
            else {
                strsA = strsList.get(readyAlign[0]);
            }
            if (readyAlign[1] < this.num) {
                strsB = new String[1];
                strsB[0] = strs[readyAlign[1]];
            }
            else {
                strsB = strsList.get(readyAlign[1]);
            }
            multiKband mkband = new multiKband(strsA, strsB, alphabet, kk);
            String[] strsC = mkband.getStrsAlign();
            strsList.put(readyAlign[2], strsC);
            strsList.remove(readyAlign[0]);
            strsList.remove(readyAlign[1]);
            labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
            labelsList.remove(readyAlign[0]);
            labelsList.remove(readyAlign[1]);
        }
        this.orders = labelsList.get(treeList[treeList.length-1][2]);
        this.straligned = strsList.get(treeList[treeList.length-1][2]);
    }
}
