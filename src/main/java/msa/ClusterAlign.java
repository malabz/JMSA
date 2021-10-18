package msa;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;

import hierCluster.guidetree;
import io.str;
import psa.FastMultiAlign;
import psa.dsa;
import psa.multiDP;
import psa.multiKband;
import strsCluster.*;
import measure.*;

public class ClusterAlign {
    private final int kk, numsStrs;
    private int[] order;
    private final char[] alphabet;
    private String[] labels;
    private String[] strsAligned;
    private HashMap<Integer, int[]> cluIdx;
    // mode1 "c" or "t" mode2 "c" or "t"
    private String mode1 = "t";
    private String mode2 = "t2";

    public String[] getStrsAlign() {
        reOrder();
        return strsAligned;
    }

    public ClusterAlign (String[] strs, boolean iter){
        numsStrs = strs.length;
        kk = score.getK(strs, false);
        kmer km = new kmer(strs, 4);
        alphabet = km.Counter();
        strsAligned = new String[strs.length];
        multiAlign(cluAlign(getCluster(strs, iter)));
        // new reAlign(strsAligned);
    }

    /**
     * get cluster from cdhitfilew
     */
    public ClusterAlign (String[] strs, String[] labels, String cdhitfile) throws IOException {
        numsStrs = strs.length;
        kk = score.getK(strs, false);
        kmer km = new kmer(strs, 4);
        alphabet = km.Counter();
        strsAligned = new String[strs.length];
        this.labels = labels;
        multiAlign(cluAlign(getCluCdhit(strs, cdhitfile)));
    }

    /**
     * reorder the strings
     */
    private void reOrder() {
        int j = 0;
        String[] res = new String[strsAligned.length];
        for (int i : order) {
            res[i] = strsAligned[j++];
        }
        strsAligned = res;
    }

    /**
     * gen the cluster by own
     */
    private HashMap<Integer, String[]> getCluster (String[] strs, boolean iter) {
        long startTime = System.currentTimeMillis();
        System.out.print("\nTo cluster ... ");
        // CenCluster cster = new CenCluster(strs, 0.85, iter);
        // cluIdx = cster.getClusters();
        FastCluster fastCluster = new FastCluster(strs, 0.85);
        cluIdx = fastCluster.getClusters();
        HashMap<Integer, String[]> multiStrs = new HashMap<>();
        for (int c : cluIdx.keySet()) {
            int[] newOne = new int[cluIdx.get(c).length + 1];
            newOne[0] = c;
            System.arraycopy(cluIdx.get(c), 0, newOne, 1, newOne.length - 1);
            cluIdx.put(c, newOne);

            String[] newStrs = new String[cluIdx.get(c).length];
            for (int i = 0; i < newStrs.length; i++) {
                newStrs[i] = strs[cluIdx.get(c)[i]];
            }
            multiStrs.put(c, newStrs);
        }
        long endTime = System.currentTimeMillis();
        System.out.println("\ntime : " + (endTime - startTime)/1000 + "s\n");
        return multiStrs;
    }

    /**
     * get cluster strings [][] by cdhitdfile
     */
    private HashMap<Integer, String[]> getCluCdhit (String[] strs, String cdhitfile) throws IOException {
        exCDhit cdhit = new exCDhit();
        String[][] names = cdhit.readClstr(cdhitfile);
        HashMap<String, Integer> idxmap = new HashMap<>();
        HashMap<Integer, String[]> multiStrs = new HashMap<>();
        cluIdx = new HashMap<>();
        Pattern p = Pattern.compile("(>\\s*\\w+).*");
        for (int i = 0; i < labels.length; i++) {
            Matcher m = p.matcher(labels[i]);
            if (m.matches()) idxmap.put(m.group(1), i);
            else throw new NumberFormatException();
        }
        for (String[] name : names) {
            String[] temp = new String[name.length];
            int[] tempidx = new int[name.length];
            for (int j = 0; j < name.length; j++) {
                tempidx[j] = idxmap.get(name[j]);
                temp[j] = strs[tempidx[j]];
            }
            multiStrs.put(tempidx[0], temp);
            cluIdx.put(tempidx[0], tempidx);
        }
        return multiStrs;
    }


    private HashMap<Integer, String[]> cluAlign (HashMap<Integer, String[]> multiStrs) {
        if (multiStrs.size() <= 1) {
            int length = 0;
            for (int c : multiStrs.keySet()) {
                length = multiStrs.get(c).length;
            }
            if (length <= 3000) { mode1 = "t"; }
            else { mode1 = "c"; }
        }
        int i = 1;
        System.out.println("To align the cluster\n");
        long startTime = System.currentTimeMillis();
        for (int c : multiStrs.keySet()) {
            String outToScreen = "    " + (i++) + " / " + multiStrs.size();
            System.out.print(outToScreen);
            if (multiStrs.get(c).length == 2) {
                dsa da = new dsa(multiStrs.get(c)[0], multiStrs.get(c)[1], "fmindex");
                multiStrs.put(c, da.getStrAlign());
            }
            else if (multiStrs.get(c).length < 2) {
                System.out.print(str.repeat("\b", outToScreen.length()));
                continue;
            }
            else if (mode1.equals("c") || multiStrs.get(c).length >= 3000) {
                centerAlign cAlign = new centerAlign(multiStrs.get(c), 1);
                multiStrs.put(c, cAlign.getStrsAlign());
            }
            else if (mode1.equals("t")) {
                treeAlign tAlign = new treeAlign(multiStrs.get(c), "upgma", 1);
                multiStrs.put(c, tAlign.getStrsAlign());
            }
            else {
                throw new IllegalArgumentException("Unknown mode " + mode1);
            }
            // Fasta ft = new Fasta();
            // String[] lables = new String[multiStrs.get(c).length];
            // try {
            //     ft.writeFasta(multiStrs.get(c), lables, "/home/kun/Documents/mywork/msa/" + c + ".cluster");
            // } catch (IOException e) {
            //     e.printStackTrace();
            // }
            System.out.print(str.repeat("\b", outToScreen.length()));
        }
        long endTime = System.currentTimeMillis();
        System.out.println("time : " + (endTime - startTime)/1000 + "s             \n");
        return multiStrs;
    }


    private void multiAlign (HashMap<Integer, String[]> multiStrsed) {
        if (numsStrs > 30000) { mode2 = "t1"; }
        switch (mode2) {
            case "t1" :
                TreeAlign1(multiStrsed);
                break;
            case "t2" :
                TreeAlign2(multiStrsed);
                break;
            default :
                throw new IllegalArgumentException("Unknown mode " + mode2);
        }
    }

    /**
     * 类中取1，再来比对
     */
    public void TreeAlign1(HashMap<Integer, String[]> multiStrs) {
        System.out.println("To combine the cluster");
        if (multiStrs.size() < 2) {
            for (int i : multiStrs.keySet()) {
                strsAligned = multiStrs.get(i);
                order = cluIdx.get(i);
            }
        }
        else {
            // cluIdx 索引映射 cluIdx [centerIdxc, ...] N + 1
            // multiStrs 类的字符串映射
            String[] centerStrs = new String[multiStrs.size()];
            int[] keys = new int[multiStrs.size()];
            int i = 0;
            for (int key : multiStrs.keySet()) {
                centerStrs[i] = PickRealOne(multiStrs.get(key));
                keys[i++] = key;
            }
            treeAlign tAlign = new treeAlign(centerStrs, "upgma", 1);
            centerStrs = tAlign.getStrsAlign();
            for (i = 0; i < centerStrs.length; i++) {
                insertGapStrings(centerStrs[i], multiStrs.get(keys[i]));
            }
            order = new int[numsStrs];
            i = 0;
            for (int key : keys) {
                int len = cluIdx.get(key).length;
                System.arraycopy(cluIdx.remove(key), 0, order, i, len);
                System.arraycopy(multiStrs.remove(key), 0, strsAligned, i, len);
                i += len;
            }
        }
        // System.out.println(" sps: " + String.format("%.3f", score.sps(strsAligned)));
        // System.out.println("  tc: " + String.format("%.3f", score.tc(strsAligned)));
        // System.out.println("  tc: " + String.format("%.3f", score.tcStrict(strsAligned)));
    }

    /**
     * 生成一颗中心序列的树，再来树比对
     */
    private void TreeAlign2(HashMap<Integer, String[]> multiStrs) {
        System.out.println("To combine the cluster\n");
        long startTime = System.currentTimeMillis();
        long time = 0;
        if (multiStrs.size() < 2) {
            for (int i : multiStrs.keySet()) {
                strsAligned = multiStrs.get(i);
                order = cluIdx.get(i);
            }
        }
        else {
            HashMap<Integer, int[]> newcluIdx = new HashMap<>();
            HashMap<Integer, Integer> idxmap = new HashMap<>();
            HashMap<Integer, String[]> NewStrsed = new HashMap<>();
            int[][] treeList = GenList(multiStrs, idxmap);
            int i = 0;
            for (int[] readyAlign : treeList) {
                String outToScreen = "    " + (i++) + " / " + treeList.length;
                System.out.print(outToScreen);
                String[] strsA, strsB;
                int[] IdxA, IdxB, IdxC;
                if (readyAlign[0] < idxmap.size()) {
                    strsA = multiStrs.remove(idxmap.get(readyAlign[0]));
                    IdxA = cluIdx.remove(idxmap.get(readyAlign[0]));
                }
                else {
                    strsA = NewStrsed.remove(readyAlign[0]);
                    IdxA = newcluIdx.remove(readyAlign[0]);
                }
                if (readyAlign[1] < idxmap.size()) {
                    strsB = multiStrs.remove(idxmap.get(readyAlign[1]));
                    IdxB = cluIdx.remove(idxmap.get(readyAlign[1]));
                }
                else {
                    strsB = NewStrsed.remove(readyAlign[1]);
                    IdxB = newcluIdx.remove(readyAlign[1]);
                }
                int ln = Math.abs(strsB[0].length() - strsA[0].length());
                double diff = (double) ln / Math.max(strsB[0].length(), strsA[0].length());
                if (diff > 0.4) {
                    multiDP mkband = new multiDP(strsA, strsB, alphabet);
                    NewStrsed.put(readyAlign[2], mkband.getStrsAlign());
                }
                else if (diff > 0.1) {
                    long startTime1 = System.currentTimeMillis();
                    multiKband mkband = new multiKband(strsA, strsB, alphabet, kk);
                    NewStrsed.put(readyAlign[2], mkband.getStrsAlign());
                    long endTime1 = System.currentTimeMillis();
                    time += endTime1 - startTime1;
                }
                else {
                    FastMultiAlign mkband = new FastMultiAlign(strsA, strsB, alphabet, kk, 0, 0);
                    NewStrsed.put(readyAlign[2], mkband.getStrsAlign());
                }
                IdxC = new int[IdxA.length + IdxB.length];
                System.arraycopy(IdxA, 0, IdxC, 0, IdxA.length);
                System.arraycopy(IdxB, 0, IdxC, IdxA.length, IdxB.length);
                newcluIdx.put(readyAlign[2], IdxC);

                System.out.print(str.repeat("\b", outToScreen.length()));
            }
            strsAligned = NewStrsed.get(treeList[treeList.length-1][2]);
            order = newcluIdx.get(treeList[treeList.length-1][2]);
        }
        long endTime = System.currentTimeMillis();
        System.out.println("time : " + (endTime - startTime)/1000 + "s             \n");
        System.out.println("time : " + (time)/1000 + "s             \n");
        // System.out.println(" sps: " + String.format("%.3f", score.sps(strsAligned)));
        // System.out.println("  tc: " + String.format("%.3f", score.tc(strsAligned)));
        // System.out.println("  tc: " + String.format("%.3f", score.tcStrict(strsAligned)));
    }


    private int[][] GenList(HashMap<Integer, String[]> multiStrs, HashMap<Integer, Integer> idxMap) {
        String[] strs = new String[multiStrs.size()];
        int tempIdx = 0;
        for (int key : multiStrs.keySet()) {
            strs[tempIdx] = multiStrs.get(key)[0];
            idxMap.put(tempIdx++, key);
        }
        guidetree gTree = new guidetree(strs, "upgma");
        return gTree.genTreeList(1);
    }


    private void insertGapStrings(String str, String[] strings) {
        int[] mark = new int[strings[0].length() + 1];
        int i = 0, nums = 0;
        for (char c : str.toCharArray()) {
            if (c == '-') nums++;
            else {
                mark[i++] = nums;
                nums = 0;
            }
        }
        mark[i] = nums;
        for (i = 0; i < strings.length; i++) {
            strings[i] = insertGap(mark, strings[i]);
        }
    }

    private String insertGap(int[] mark, String seq) {
        assert mark.length == seq.length() + 1;
        StringBuilder seqGap = new StringBuilder();
        int len = mark.length;
        for (int i = 0; i < len; i++) {
            seqGap.append(str.repeat("-", mark[i]));
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
}
