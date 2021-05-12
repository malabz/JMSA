package psa;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.ArrayList;

/** 
 * Suffix tree align
 * 
 * @author Juntao chen
*/
public class STAlign extends STree {

    private String B;
    private String alignA, alignB;
    private float sim;
    private int thorehold;

    /**
     * To align A and B
     * @param A
     * @param B
     */
    public STAlign (String A, String B) {
        super(A);
        this.B = B;
        this.buildTree();
        Align();
    }

    /**
     * To build the suffix tree for String A
     * @param A
     */
    public STAlign (String A) {
        super(A);
        this.buildTree();
    }

    /**
     * input string B to align string A
     * @param B
     */
    public void AlignStrB (String B) {
        this.B = B;
        Align();
    }

    /**
     * To get the align results.
     * @return
     */
    public String[] getStrAlign() {
        String[] res = new String[2];
        res[0] = alignA;
        res[1] = alignB;
        return res;
    }

    /**
     * To get the sim of A and B
     * @return
     */
    public float getSimstrB(String B) {
        this.B = B;
        thorehold = this.len/100 > 1 ? this.len/100 : 15;
        float[] sim = new float[]{0};
        this.slectCommonStrings(sim);
        this.sim = sim[0];
        return this.sim;
    }

    private boolean walkdownFcs(Node node, int step) {
        return super.edgelength(node) <= step ? true : false;
    }

    private void dfsleaves(Node node, List<Integer> results ,int length) {
        Collection<Node> sons = node.children.values();
        for (Node son : sons) {
            if (son.leaf && super.edgelength(son) > 1) {
                results.add(son.start - length);
            }
            else {
                int tmp = length + this.edgelength(son);
                this.dfsleaves(son, results, tmp);
            }
        }
    }

    private List<Integer> selectprefix(String s) {
        Node node = this.root;
        List<Integer>  starts = new ArrayList<>();
        int length = 0, step = 0, tag = 0;
        while (node.children.get(s.charAt(step + length)) != null) {
            node = node.children.get(s.charAt(step + length));
            step += 1;

            if (s.length() >= step + length) {
                if (this.walkdownFcs(node, step)) {
                    length += step;
                    step = 0;
                    if (s.length() == length) {break;}
                    else {continue;}
                }
                if (s.length() == step + length) {break;}
            }
            else {break;}

            if (this.T.charAt(node.start + step) != s.charAt(step + length)) { break; }

            while (this.T.charAt(node.start + step) == s.charAt(step + length)) {
                step += 1;
                if (s.length() >= step + length) {
                    if (this.walkdownFcs(node, step)) {
                        length += step;
                        step = 0;
                        if (s.length() == length) { tag = 1; }
                        break;
                    }
                    if (s.length() == step + length) {
                        tag = 1;
                        break;
                    }
                }
                else {
                    tag = 1;
                    break;
                }
            }

            if (tag == 1) { break; }
        }
        
        if (node.leaf) {starts.add(node.start - length);}
        else if (length + step >= thorehold) {this.dfsleaves(node, starts, length);}
        starts.add(length+step);
        return starts;
    }

    private int MultiReg(List<Integer> l, int len, float rate) {
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

    private List<int[]> findCommonStrings () {
        int index = 0;
        // int int
        List<int[]> loc2 = new ArrayList<>();
        List<List<Integer>> locx = new ArrayList<>();

        while (index <= this.B.length() - 1) {
            List<Integer> results = this.selectprefix(this.B.substring(index));
            int length = results.remove(results.size()-1);
            if (results.size() != 0 && length > thorehold ) {
                int[] tmp1 = {index, length};
                loc2.add(tmp1);
                locx.add(results);
                index += length;
            }
            else { index += 1; }
        }
        List<int[]> deresults = new ArrayList<>();
        if (loc2.size() != 0) {
            float start = Float.POSITIVE_INFINITY;
            float end = Float.NEGATIVE_INFINITY;
            int length1 = loc2.get(loc2.size()-1)[0] - loc2.get(loc2.size()-1)[1];
            int st = 0;
            for (int i = 0; i + st < locx.size(); ++i) {
                if (locx.get(i) == null || loc2.get(i)[1] < length1 / 100) {
                    locx.remove(i);
                    loc2.remove(i);
                    --i;
                    st++;
                }
            }
            for (List<Integer> l : locx) {
                if (Collections.max(l) > end) { end = Collections.max(l);}
                if (Collections.min(l) < start) { start = Collections.min(l);}
            }
            int length2 = (int)(end - start);
            for (int i = 0; i < locx.size(); ++i) {
                int [] tmp = {0,0,0};
                if (locx.get(i).size() > 1) {
                    int idx = this.MultiReg(locx.get(i), length2, (float)loc2.get(i)[0]/length1);
                    int [] tmp1 = {loc2.get(i)[0], loc2.get(i)[1], locx.get(i).get(idx)};
                    tmp = tmp1;
                }
                else {
                    int [] tmp1 = {loc2.get(i)[0], loc2.get(i)[1], locx.get(i).get(0)};
                    tmp = tmp1;
                }
                deresults.add(tmp);
            }
        }
        return deresults;
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

    private int[] traceback (List<int[]> results, int[] p) {
        int idx = this.getIdxMax(p);
        int j;
        List<Integer> track = new ArrayList<Integer>();
        track.add(idx);
        
        while (idx > 0) {
            j = idx - 1;
            if (p[idx] == results.get(idx)[1]) {
                break;
            }
            while (j >= 0) {
                if ((results.get(idx)[2] >= results.get(j)[1] + results.get(j)[2]) && p[idx] == p[j] + results.get(idx)[1]) {
                    track.add(j);
                    idx = j;
                    break;
                }
                j--;
            }
        }
        Collections.reverse(track);
        return track.stream().mapToInt(Integer::valueOf).toArray();
    }

    private int[][] slectCommonStrings(float[] sim) {
        List<int[]> results = this.findCommonStrings();
        if (results.size() != 0) {
            int m = results.size();
            int[] p = new int[m];
            for (int i = 0; i < m; ++i) {
                p[i] = results.get(i)[1];
            }
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < i; ++j) {
                    if (results.get(i)[2] >= results.get(j)[2] + results.get(j)[1]) {
                        p[i] = Math.max(p[i], p[j] + results.get(i)[1]);
                    }
                }
            }
            int[] seIdx = this.traceback(results, p);
            int[][] seresluts = new int[seIdx.length][3];
            int j = 0;
            for (int idx : seIdx) {
                seresluts[j++] = results.get(idx);
            }
            sim[0] = (float)p[this.getIdxMax(p)]/this.len;
            return seresluts;
        }
        return null;
    }

    public void Align() {
        thorehold = this.len/100 > 1 ? this.len/100 : 15;
        float[] sim = new float[]{0};
        int[][] results = this.slectCommonStrings(sim);
        this.sim = sim[0];
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
                sA.append(this.A.substring(results[j][2], results[j][2] + results[j][1]));
                sB.append(this.B.substring(results[j][0], results[j][0] + results[j][1]));
            }
            sA.append(alStrsA.get(alStrsA.size()-1));
            sB.append(alStrsB.get(alStrsB.size()-1));

            this.alignA = sA.toString();
            this.alignB = sB.toString();
        }
    }
}
