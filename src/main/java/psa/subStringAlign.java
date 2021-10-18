package psa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This is the abstract class of the sub string align.
 */
public abstract class subStringAlign {

    protected String A, B, alignA, alignB;
    protected int thorehold;
    private float sim;
    protected boolean state = false;

    /**
     * To get the aligned results.
     */
    public String[] getStrAlign() {
        if (state) return new String[]{alignB, alignA};
        return new String[]{alignA, alignB};
    }

    /**
     * To get the sim of A and B
     */
    public float getSimstrB(String B) {
        this.B = B.toLowerCase();
        this.selectCommonStrings();
        return this.sim;
    }

    /**
     * string B to align string A
     */
    public void AlignStrB(String B) {
        this.B = B.toLowerCase();
        this.Align();
    }

    /**
     * get substrings between A and B
     * @param B
     */
    public int[][] getSubStrings(String B) {
        this.B = B.toLowerCase();
        return selectCommonStrings();
    }

    protected abstract List<Integer> selectprefix(String p);

    protected int MultiReg(List<Integer> l, int len, float rate) {
        int idx = 0;
        float ratex = Float.POSITIVE_INFINITY;
        for (int i = 0; i < l.size(); ++i) {
            float tmp = Math.abs((float) l.get(i) / len - rate);
            if (tmp < ratex) {
                idx = i;
                ratex = tmp;
            }
        }
        return idx;
    }

    protected abstract int[][] findCommonStrings();

    protected int[][] pickOneSubString(List<int[]> IdxLen, List<List<Integer>> locx) {
        List<int[]> res = new ArrayList<>();
        
        if (IdxLen.size() == 0) { return new int[0][]; }
        
        float start = Float.POSITIVE_INFINITY;
        float end = Float.NEGATIVE_INFINITY;
        for (List<Integer> l : locx) {
            if (Collections.max(l) > end) { end = Collections.max(l);}
            if (Collections.min(l) < start) { start = Collections.min(l);}
        }
        int length1 = IdxLen.get(IdxLen.size()-1)[0] + IdxLen.get(IdxLen.size()-1)[1] - IdxLen.get(0)[0];
        int length2 = (int) (end - start);
        for (int i = 0; i < locx.size(); ++i) {
            if (locx.get(i).size() > 1) {
                int idx = this.MultiReg(locx.get(i), length2, (float)IdxLen.get(i)[0]/length1);
                res.add(new int[]{IdxLen.get(i)[0], IdxLen.get(i)[1], locx.get(i).get(idx)});
            }
            else {
                res.add(new int[]{IdxLen.get(i)[0], IdxLen.get(i)[1], locx.get(i).get(0)});
            }
        }
        return res.toArray(new int[res.size()][]);
    }

    private int getIdxMax (int[] p) {
        int idx = 0, vmax = p[0];
        for (int i = 1; i < p.length; ++i) {
            if (p[i] > vmax) {
                idx = i;
                vmax = p[i];
            }
        }
        return idx;
    }

    private int[] traceback (int[][] results, int[] p) {
        int idx = this.getIdxMax(p), j;
        List<Integer> track = new ArrayList<>();
        track.add(idx);
        
        while (idx > 0) {
            j = idx - 1;
            if (p[idx] == results[idx][1]) { break; }
            while (j >= 0) {
                if ((results[idx][2] >= results[j][1] + results[j][2]) && p[idx] == p[j] + results[idx][1]) {
                    track.add(j);
                    idx = j;
                    break;
                }
                j--;
            }
        }
        int[] res = new int[track.size()];
        int i = res.length - 1;
        for (int t : track) { res[i--] = t; }
        return res;
    }

    protected int[][] selectCommonStrings() {
        if (A.length() >= B.length()) {
            this.thorehold = Math.max(this.A.length() / 100, 15);
        }
        else {
            this.thorehold = Math.max(this.B.length() / 100, 15);
        }
        
        int[][] results = findCommonStrings();
        if (results.length == 0) { return null; }

        int m = results.length;
        int[] p = new int[m];
        for (int i = 0; i < m; ++i) {
            p[i] = results[i][1];
        }
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < i; ++j) {
                if (results[i][2] >= results[j][2] + results[j][1]) {
                    p[i] = Math.max(p[i], p[j] + results[i][1]);
                }
            }
        }
        int[] seIdx = this.traceback(results, p);
        int[][] res = new int[seIdx.length][3];
        int j = 0;
        for (int idx : seIdx) { res[j++] = results[idx]; }
        this.sim = (float)p[this.getIdxMax(p)]/this.B.length();
        return res;
    }

    protected void Align() {
        int[][] results = this.selectCommonStrings();
        if (results == null) {
            Kband alignKb = new Kband(this.A, this.B);
            this.alignA = alignKb.getStrAlign()[0];
            this.alignB = alignKb.getStrAlign()[1];
        }
        else {
            List<String> readyStrsA = new ArrayList<>();
            List<String> readyStrsB = new ArrayList<>();
            List<String> alStrsA = new ArrayList<>();
            List<String> alStrsB = new ArrayList<>();

            int sa = 0, sb = 0;

            for (int[] r : results) {
                readyStrsA.add(this.A.substring(sa, r[2]));
                readyStrsB.add(this.B.substring(sb, r[0]));
                sa = r[1] + r[2];
                sb = r[1] + r[0];
            }

            readyStrsA.add(sa >= this.A.length() ? "" : this.A.substring(sa));
            readyStrsB.add(sb >= this.B.length() ? "" : this.B.substring(sb));

            for (int i = 0; i < readyStrsA.size(); ++i) {
                Kband alignKb = new Kband(readyStrsA.get(i), readyStrsB.get(i));
                alStrsA.add(alignKb.getStrAlign()[0]);
                alStrsB.add(alignKb.getStrAlign()[1]);
            }
            
            StringBuilder sA = new StringBuilder();
            StringBuilder sB = new StringBuilder();

            for (int j = 0; j < results.length; ++j) {
                sA.append(alStrsA.get(j));
                sB.append(alStrsB.get(j));
                sA.append(this.A, results[j][2], results[j][2] + results[j][1]);
                sB.append(this.B, results[j][0], results[j][0] + results[j][1]);
            }
            sA.append(alStrsA.get(alStrsA.size()-1));
            sB.append(alStrsB.get(alStrsB.size()-1));

            this.alignA = sA.toString();
            this.alignB = sB.toString();
        }
    }
    
}

