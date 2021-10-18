package psa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import index.FMIndex;

public class FMAlign extends subStringAlign {
    private final FMIndex FMA;

    /**
     * To build bwt for string A
     * @param A
     */
    public FMAlign(String A) {
        this.A = A.toLowerCase();
        this.FMA = new FMIndex(this.A);
    }

    public FMAlign(String A, String B) {
        if (A.length() >= B.length()) {
            this.A = A.toLowerCase();
            this.B = B.toLowerCase();
        }
        else { 
            this.A = B.toLowerCase();
            this.B = A.toLowerCase();
            this.state = true;
        }
        this.FMA = new FMIndex(this.A);
        this.Align();
    }

    @Override
    protected List<Integer> selectprefix(String p) {
        int l = 0, r = FMA.getlen() - 1;
        int oldL = l, oldr = r + 1;
        int length = p.length();
        int i;
        for (i = length - 1; i >= 0; i--) {
            l = FMA.getgRank().rank(p.charAt(i), l-1) + FMA.count(p.charAt(i));
            r = FMA.getgRank().rank(p.charAt(i), r) + FMA.count(p.charAt(i)) - 1;
            if (r < l) { break;}
            oldL = l;
            oldr = r + 1;
        }
        List<Integer> res = new ArrayList<>();
        if (length - 1 - i < thorehold) { return res; }
        else {
            for (int o = oldL; o < oldr; o++) { res.add(FMA.resolve(o)); }
        }
        res.sort(Comparator.naturalOrder());

        res.add(length - 1 - i);
        return res;
    }

    @Override
    protected int[][] findCommonStrings() {
        int index = this.B.length();
        List<int[]> IdxLen = new ArrayList<>();
        List<List<Integer>> locx = new ArrayList<>();

        while (index > 0) {
            List<Integer> idxsPre = this.selectprefix(this.B.substring(0, index));
            if (idxsPre.size() > 0) {
                int length = idxsPre.remove(idxsPre.size() - 1);
                IdxLen.add(new int[]{index-length, length});
                locx.add(idxsPre);
                index -= length;
            }
            else { index--; }
        }
        List<int[]> res = new ArrayList<>();
        
        if (IdxLen.size() == 0) { return res.toArray(new int[0][]); }

        float start = Float.POSITIVE_INFINITY;
        float end = Float.NEGATIVE_INFINITY;
        for (List<Integer> l : locx) {
            if (Collections.max(l) > end) { end = Collections.max(l);}
            if (Collections.min(l) < start) { start = Collections.min(l);}
        }
        int length1 = IdxLen.get(0)[0] - IdxLen.get(IdxLen.size()-1)[0];
        int length2 = (int) (end - start);
        
        for (int i = locx.size() - 1; i >= 0; i--) {
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
}
