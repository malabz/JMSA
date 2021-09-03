package index;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class GetRank {
    private static final int gap = 8;
    private final String T;
    private final HashMap<Character, int[]> ranks = new HashMap<>();
    private final HashMap<Character, Integer> count = new HashMap<>();

    public GetRank (String T) {
        this.T = T.toLowerCase();
        this.gen();
    }

    private void gen () {
        HashMap<Character, List<Integer>> rankList = new HashMap<>();
        for (char c : T.toCharArray()) {
            if (!count.containsKey(c)) {
                count.put(c, 0);
                rankList.put(c, new ArrayList<>());
            }
        }
        int i = 0;
        for (char c : T.toCharArray()) {
            count.put(c, count.get(c) + 1);
            if (i++ % gap == 0) {
                for (char cr : count.keySet()) {
                    rankList.get(cr).add(count.get(cr));
                }
            }
        }
        for (char c : rankList.keySet()) {
            List<Integer> tempList = rankList.get(c);
            int[] tempInt = new int[tempList.size()];
            i = 0;
            for (int j : tempList) tempInt[i++] = j;
            ranks.put(c, tempInt);
        }
    }

    /**
     * give the char and offset to get the rank
     * @param c
     * @param offset
     * @return rank
     */
    public int rank (char c, int offset) {
        if (offset < 0 || !count.containsKey(c)) return 0;
        int i = offset;
        int numc = 0;
        while (i % gap != 0) {
            if (T.charAt(i) == c) numc += 1;
            i -= 1;
        }
        return ranks.get(c)[i / gap] + numc;
    }
}
