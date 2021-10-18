package psa;

import java.util.Collection;
import java.util.ArrayList;
import java.util.List;

import index.Node;
import index.STree;

/** 
 * Suffix tree align
 * 
 * @author Juntao chen
*/
public class STAlign extends subStringAlign{
    private final STree STA;

    /**
     * To align A and B
     * @param A
     * @param B
     */
    public STAlign (String A, String B) {
        if (A.length() >= B.length()) {
            this.A = A.toLowerCase();
            this.B = B.toLowerCase();
        }
        else { 
            this.A = B.toLowerCase();
            this.B = A.toLowerCase();
            this.state = true;
        }
        this.STA = new STree(this.A);
        Align();
    }

    /**
     * To build the suffix tree for String A
     * @param A
     */
    public STAlign (String A) {
        this.A = A.toLowerCase();
        this.STA = new STree(this.A);
    }

    private boolean walkdownFcs(Node node, int step) {
        return STA.edgelength(node) <= step;
    }

    private void dfsleaves(Node node, List<Integer> results ,int length) {
        Collection<Node> sons = node.children.values();
        for (Node son : sons) {
            if (son.leaf) {
                results.add(son.start - length - STA.edgelength(node));
            }
            else {
                int tmp = length + STA.edgelength(son);
                this.dfsleaves(son, results, tmp);
            }
        }
    }

    @Override
    protected List<Integer> selectprefix(String s) {
        Node node = STA.root;
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

            if (STA.T.charAt(node.start + step) != s.charAt(step + length)) { break; }

            while (STA.T.charAt(node.start + step) == s.charAt(step + length)) {
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
            if (STA.T.charAt(node.start + step) != s.charAt(step+length)) { break; }
        }
        
        if (node.leaf) {starts.add(node.start - length);}
        else if (length + step >= thorehold) {this.dfsleaves(node, starts, length);}
        starts.add(length+step);
        return starts;
    }

    @Override
    protected int[][] findCommonStrings () {
        int index = 0;
        // int int
        List<int[]> IdxLen = new ArrayList<>();
        List<List<Integer>> locx = new ArrayList<>();

        while (index <= this.B.length() - 1) {
            List<Integer> results = this.selectprefix(this.B.substring(index));
            int length = results.remove(results.size()-1);
            if (results.size() != 0 && length > thorehold ) {
                IdxLen.add(new int[]{index, length});
                locx.add(results);
                index += length;
            }
            else { index++; }
        }
        return pickOneSubString(IdxLen, locx);
    }
}
