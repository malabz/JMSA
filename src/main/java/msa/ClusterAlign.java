package msa;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;

import hierCluster.guidetree;
// import io.fasta;
import psa.dsa;
import psa.multiKband;
import strsCluster.*;
import measure.*;

public class ClusterAlign {
    private int kk;
    private int[] order;
    private final char[] alphabet;
    private String[] labels;
    private String[] strsAligned;
    private HashMap<Integer, int[]> cluIdx;
    // mode1 "c" or "t" mode2 "c" or "t"
    private String mode1 = "t";
    private String mode2 = "t";

    public String[] getStrsAlign() {
        reOrder();
        return strsAligned;
    }

    public ClusterAlign (String[] strs){
        score sc = new score();
        kk = sc.getK(strs, false);
        kmer km = new kmer(strs);
        alphabet = km.Counter();
        strsAligned = new String[strs.length];
        multiAlign(cluAlign(getCluster(strs)));
        // new reAlign(strsAligned);
    }

    /**
     * get cluster from cdhitfile
     * @param strs
     * @param labels
     * @param cdhitfile
     */
    public ClusterAlign (String[] strs, String[] labels, String cdhitfile) throws IOException {
        kmer km = new kmer(strs);
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
     * @return
     */
    private HashMap<Integer, String[]> getCluster (String[] strs) {
        long startTime = System.currentTimeMillis();
        System.out.print("\nTo cluster ... ");
        CenCluster cster = new CenCluster(strs, 0.9);
        cluIdx = cster.getClusters();
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
     * @param cdhitfile
     * @return
     * @throws IOException
     */
    private HashMap<Integer, String[]> getCluCdhit (String[] strs, String cdhitfile) throws IOException {
        exCDhit cdhit = new exCDhit();
        String[][] names = cdhit.readClstr(cdhitfile);
        HashMap<String, Integer> idxmap = new HashMap<>();
        HashMap<Integer, String[]> multiStrs = new HashMap<>();
        cluIdx = new HashMap<>();
        Pattern p = Pattern.compile("(\\D\\w+).*");
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
        if (multiStrs.size() == 1) { mode1 = "c"; }
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
                System.out.print("\b".repeat(outToScreen.length()));
                continue;
            }
            else if (mode1.equals("c")) {
                centerAlign cAlign = new centerAlign(multiStrs.get(c), 1);
                multiStrs.put(c, cAlign.getStrsAlign());
            }
            else if (mode1.equals("t")) {
                treeAlign tAlign = new treeAlign(multiStrs.get(c), "upgma", 1);
                multiStrs.put(c, tAlign.getStrsAlign());
            }
            // fasta ft = new fasta();
            // String[] lables = new String[multiStrs.get(c).length];
            // try {
            //     ft.writeFasta(multiStrs.get(c), lables, "/home/kun/Documents/mywork/msa/code/" + c + ".cluster");
            // } catch (IOException e) {
            //     // TODO Auto-generated catch block
            //     e.printStackTrace();
            // }
            System.out.print("\b".repeat(outToScreen.length()));
        }
        long endTime = System.currentTimeMillis();
        System.out.println("time : " + (endTime - startTime)/1000 + "s             \n");
        return multiStrs;
    }


    private void multiAlign (HashMap<Integer, String[]> multiStrsed) {
        if (mode2.equals("c")) {
            // cenAlign1(multiStrsed);
            cenAlign2(multiStrsed);
        }
        else if (mode2.equals("t")) {
            TreeAlign2(multiStrsed);
        }
    }

    /**
     * 顺序的树比对
     */
    public void TreeAlign1(HashMap<Integer, String[]> multiStrs) {

    }

    /**
     * 生成一颗中心序列的树，再来树比对
     * 1. 既然是多序列的那就生成多序列的树来做！
     */
    private void TreeAlign2(HashMap<Integer, String[]> multiStrs) {
        System.out.println("To combine the cluster\n");
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
            for (int[] readyAlign : treeList) {
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

                multiKband mkband = new multiKband(strsA, strsB, alphabet, kk);
                NewStrsed.put(readyAlign[2], mkband.getStrsAlign());
                IdxC = new int[IdxA.length + IdxB.length];
                System.arraycopy(IdxA, 0, IdxC, 0, IdxA.length);
                System.arraycopy(IdxB, 0, IdxC, IdxA.length, IdxB.length);
                newcluIdx.put(readyAlign[2], IdxC);
                
                // fasta ft = new fasta();
                // String[] lables = new String[NewStrsed.get(readyAlign[2]).length];
                // try {
                //     ft.writeFasta(NewStrsed.get(readyAlign[2]), lables, "/home/kun/Documents/mywork/msa/code/" + readyAlign[2] + ".group");
                // } catch (IOException e) {
                //     // TODO Auto-generated catch block
                //     e.printStackTrace();
                // }

            }
            strsAligned = NewStrsed.get(treeList[treeList.length-1][2]);
            order = newcluIdx.get(treeList[treeList.length-1][2]);
        }
        score sc = new score();
        System.out.println(" sps: " + String.format("%.3f", sc.sps(strsAligned)));
        System.out.println("  tc: " + String.format("%.3f", sc.tc(strsAligned)));
    }


    private int[][] GenList(HashMap<Integer, String[]> multiStrs, HashMap<Integer, Integer> idxMap) {
        String[][] strings = new String[multiStrs.size()][];
        int tempIdx = 0;
        for (int key : multiStrs.keySet()) {
            strings[tempIdx] = multiStrs.get(key);
            idxMap.put(tempIdx++, key);
        }
        guidetree gTree = new guidetree(strings, "upgma");
        return gTree.genTreeListND();
    }

    /**
     * 中心比对1
     */
    public void cenAlign1(HashMap<Integer, String[]> multiStrs) {
        HashMap<Integer, Integer> idxName = new HashMap<>();
        String[] cenStrs = new String[multiStrs.size()];
        int tempi = 0;
        for (int key : multiStrs.keySet()) {
            cenStrs[tempi] = PickRealOne(multiStrs.get(key));
            idxName.put(tempi++, key);
        }
        
        String[] cenStrsed = new String[cenStrs.length];
        
        if (cenStrs.length > 2) {
            // TODO : choose 1
            centerAlign cAlign = new centerAlign(cenStrs, 1);
            cenStrsed = cAlign.getStrsAlign();
            // treeAlign tAlign = new treeAlign(cenStrs, 1);
            // cenStrsed = tAlign.getStrsAlign();
        }
        else if (cenStrs.length > 1) {
            multiKband mKband = new multiKband(multiStrs.get(0), multiStrs.get(1), alphabet, kk);
            strsAligned = mKband.getStrsAlign();
        }
        else {
            strsAligned = multiStrs.get(0);
        }
        
        for (int i = 0; i < cenStrs.length; i++) {
            int j = 0, counter = 0;
            int[] marks = new int[cenStrs[i].length() + 1];
            char[] tempc = cenStrsed[i].toCharArray();
            for (char c : tempc) {
                if (c == '-') counter++;
                else { marks[j++] = counter; counter = 0; }
            }
            marks[j] = counter;
            int idxD = 0;
            for (String str : multiStrs.get(idxName.get(i))) {
                strsAligned[cluIdx.get(idxName.get(i))[idxD++]] = insertGap(marks, str);
            }
        }
        score sc = new score();
        System.out.println(" sps: " + String.format("%.3f", sc.sps(strsAligned)));
        System.out.println("  tc: " + String.format("%.3f", sc.tc(strsAligned)));
    }

    /**
     * 中心比对2
     */
    public void cenAlign2 (HashMap<Integer, String[]> multiStrs) {
        int idxc = -1, length = 0;
        // 给每个组加上了一个guawazi,所以之后要去除掉瓜娃子，crazy
        for (int key : multiStrs.keySet()) {
            idxc = multiStrs.get(key)[0].length() > length ? key : idxc;
            length = Math.max(multiStrs.get(key)[0].length(), length);
            String[] newGUAWAZI = new String[multiStrs.get(key).length + 1];
            newGUAWAZI[0] = PickRealOne(multiStrs.get(key));
            System.arraycopy(multiStrs.get(key), 0, newGUAWAZI, 1, newGUAWAZI.length - 1);
            multiStrs.put(key, newGUAWAZI);
        }
        String[][] newStrsed = new String[multiStrs.size() - 1][2];
        HashMap<Integer, Integer> idxName = new HashMap<>();
        int tempi = 0;
        for (int key : multiStrs.keySet()) {
            if (key == idxc) continue; 
            multiKband mKband = new multiKband(multiStrs.get(idxc), multiStrs.get(key), alphabet, true);
            newStrsed[tempi][0] = mKband.getStrAlign()[0][0];
            newStrsed[tempi][1] = mKband.getStrAlign()[1][0];
            idxName.put(key, tempi++);
            
        }
        String centerSeqs = multiStrs.get(idxc)[0];
        int[] markInsertion = new int[centerSeqs.length() + 1];
        for (String[] strs2 : newStrsed) {
            int i = 0, counter = 0;
            for (char c : strs2[0].toCharArray()) {
                if (c == '-') counter++;
                else {
                    markInsertion[i] = Math.max(markInsertion[i], counter);
                    counter = 0;
                    i++;
                }
            }
            markInsertion[i] = Math.max(markInsertion[i], counter);
        }
        for (int i = 1; i < multiStrs.get(idxc).length; i++) {
            strsAligned[cluIdx.get(idxc)[i - 1]] = insertGap(markInsertion, multiStrs.get(idxc)[i]);
        }
        for (int key : multiStrs.keySet()) {
            if (key == idxc) continue;

            char[] tempA = multiStrs.get(key)[0].toCharArray();
            char[] tempB = newStrsed[idxName.get(key)][0].toCharArray();
            int[] mark1 = new int[tempA.length + 1];
            int[] mark2 = new int[tempB.length + 1];

            int i = 0, counter = 0;
            for (char c : newStrsed[idxName.get(key)][1].toCharArray()) {
                if (c == '-') counter++;
                else {
                    mark1[i] = Math.max(mark1[i], counter);
                    counter = 0;
                    i++;
                }
            }
            mark1[i] = Math.max(mark1[i], counter);
            
            int pi = 0, pj = 0, total = 0;
            for (char c : tempB) {
                if (c == '-') total++;
                else {
                    mark2[pi++] = markInsertion[pj++] - total;
                    while (total != 0) { pi++; total--; }
                }
            }
            mark2[pi] = markInsertion[pj] - total;
            for (i = 1; i < multiStrs.get(key).length; i++) {
                strsAligned[cluIdx.get(key)[i - 1]] =  insertGap(mark2, insertGap(mark1, multiStrs.get(key)[i]));
            }
        }
        score sc = new score();
        System.out.println(" sps: " + String.format("%.3f", sc.sps(strsAligned)));
        System.out.println("  tc: " + String.format("%.3f", sc.tc(strsAligned)));
    }


    private String insertGap(int[] mark, String seq) {
        assert mark.length == seq.length() + 1;
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
}
