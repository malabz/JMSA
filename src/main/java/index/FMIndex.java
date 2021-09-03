package index;

import java.util.HashMap;

public class FMIndex {
    private final int sample = 8;
    private final int len;
    private final String L;
    private final HashMap<Character, Integer> F;
    private final GetRank gRank;
    private HashMap<Integer, Integer> sampleSA;


    public FMIndex (String str) {
        str = str.toLowerCase();
        if (str.charAt(str.length() - 1) != '$') str += '$';
        BWT bt = new BWT();
        int[] sa = bt.suffixArray(str);
        this.L = bt.bwt(str, sa);
        this.F = bt.getFirstCol(L);
        this.len = L.length();
        this.downSampleSA(sa);
        this.gRank = new GetRank(L);
    }

    public int getlen() {
        return this.len;
    }

    public GetRank getgRank() {
        return this.gRank;
    }

    private void downSampleSA(int[] sa) {
        sampleSA = new HashMap<>();
        int i = 0;
        for (int idx : sa) {
            if (idx % sample == 0) sampleSA.put(i, idx);
            i++;
        }
    }

    public int count (char c) {
        return F.getOrDefault(c, 0);
    }

    private int[] range(String p) {
        p = p.toLowerCase();
        int l = 0, r = len - 1;
        int length = p.length();
        for (int i = length - 1; i >= 0; i--) {
            l = gRank.rank(p.charAt(i), l-1) + count(p.charAt(i));
            r = gRank.rank(p.charAt(i), r) + count(p.charAt(i)) - 1;
            if (r < l) break;
        }
        return new int[]{l, r + 1};
    }

    /**
     * Step left according to character in given BWT row
     * @param row
     * @return row
     */
    private int stepLeft (int row) {
        char c = L.charAt(row);
        return gRank.rank(c, row - 1) + count(c);
    }

    /**
     * Given BWM row, return its offset
     * @param row
     * @return offset
     */
    public int resolve (int row) {
        int nsteps = 0;
        while (!sampleSA.containsKey(row)) {
            row = stepLeft(row);
            nsteps += 1;
        }
        return sampleSA.get(row) + nsteps;
    }

    public boolean hasSubstring(String p) {
        int[] res = range(p);
        return res[1] > res[0];
    }

    /**
     * Return offsets for all occurrences of p
     * @param p
     * @return
     */
    public int[] occurrences(String p) {
        p = p.toLowerCase();
        int[] rows = range(p);
        if (rows[1] <= rows[0]) return new int[0];
        int[] res = new int[rows[1] - rows[0]];
        int i = 0;
        for (int r = rows[0]; r < rows[1]; r++) {
            res[i++] = resolve(r);
        }
        System.out.println("L : " + L.toUpperCase());
        return res;
    }
}