package msa;

import java.util.HashMap;
import cluster.guidetree;
import cluster.kmer;
import measure.score;
import psa.dsa;
import psa.multiKband;
public class treeAlign {
    private String[] strs;
    private String[] straligned;
    private String Treemode;
    private int[] orders;
    private int num;


    /**
     * 
     * @param strs
     */
    public treeAlign(String[] strs, String treemode) {
        this.strs = strs;
        this.num = strs.length;
        this.Treemode = treemode;
        Align();
        reOrder();
    }

    /**
     * 
     * @param strs
     */
    public treeAlign(String[] strs) {
        this.strs = strs;
        this.num = strs.length;
        this.Treemode = "nj";
        Align();
        reOrder();
    }

    /**
     * To get the alignment results.
     * @return
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

    private void Align() {
        
        System.out.println();
        System.out.println("build the " + Treemode + " Tree");
        System.out.println();

        long startTime = System.currentTimeMillis();
        guidetree gTree = new guidetree(this.strs, this.Treemode);
        int[][] treeList = gTree.genTreeList();
        var strsList = new HashMap<Integer, String[]>();
        var labelsList = new HashMap<Integer, int[]>();
        long endTime = System.currentTimeMillis();

        System.out.println("time: "+((endTime-startTime)/1000)+"s\n");
        System.out.println("Align the seqs along the Tree\n");
        
        startTime = System.currentTimeMillis();
        kmer km = new kmer(this.strs);
        char[] alphabet = km.Counter();
        int len = treeList.length, i = 0;
        for (int[] readyAlign : treeList) {
            String outToScreen = "    " + (i + 1) + " / " + len;
            System.out.print(outToScreen);
            String[] strsA;
            String[] strsB;
            if (readyAlign[0] < this.num && readyAlign[1] < this.num) {
                startTime = System.currentTimeMillis();
                dsa pa = new dsa(strs[readyAlign[0]], strs[readyAlign[1]]);
                String[] strsC = pa.getStrAlign();
                
                strsList.put(readyAlign[2], strsC);
                strsList.remove(readyAlign[0]);
                strsList.remove(readyAlign[1]);
                
                labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
                labelsList.remove(readyAlign[0]);
                labelsList.remove(readyAlign[1]);
                endTime = System.currentTimeMillis();

                i++;
                System.out.print("\b".repeat(outToScreen.length()));
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

            multiKband mkband = new multiKband(strsA, strsB, alphabet);
            String[] strsC = new String[strsA.length + strsB.length];
            System.arraycopy(mkband.alignA, 0, strsC, 0, strsA.length);
            System.arraycopy(mkband.alignB, 0, strsC, strsA.length, strsB.length);

            strsList.put(readyAlign[2], strsC);
            strsList.remove(readyAlign[0]);
            strsList.remove(readyAlign[1]);
            
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
}
