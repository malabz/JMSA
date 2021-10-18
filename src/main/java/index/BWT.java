package index;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;


class IntIdx implements Comparable<IntIdx> {
    int index;
    int num1;
    int num2;

    public IntIdx(int idx, int num1, int num2){
        this.index = idx;
        this.num1 = num1;
        this.num2 = num2;
    }

    public int compareTo(IntIdx I){
        if (num1 != I.num1) return Integer.compare(num1, I.num1);
        return Integer.compare(num2, I.num2);
    }

    public boolean equals(IntIdx I) {
        return num1 == I.num1 && num2 == I.num2;
    }
}

public class BWT {

    /**
     * burrows-wheeler transform
     * @param str
     * @return L
     */
    public String bwt(String str) {
        if (str.charAt(str.length() - 1) != '$') {
            str = str.toLowerCase() + "$";
        }
        StringBuilder L = new StringBuilder();
        int[] sa = suffixArray(str);
        for (int i : sa) {
            if (i == 0) L.append("$");
            else L.append(str.charAt(i - 1));
        }
        return L.toString();
    }

    /**
     * burrows-wheeler transform
     * @param str
     * @param sa
     * @return L
     */
    public String bwt(String str, int[] sa) {
        assert str.length() == sa.length;
        StringBuilder L = new StringBuilder();
        for (int i :sa) {
            if (i == 0) L.append("$");
            else L.append(str.charAt(i - 1));
        }
        return L.toString().toLowerCase();
    }

    /**
     * Prefix doubling algorithms to get suffix array
     * @param str
     * @return sa
     */
    public int[] suffixArray(String str) {
        if (str.charAt(str.length() - 1) != '$') str = str.toLowerCase() + "$";
        int len = str.length();

        int[] indx = new int[len];
        IntIdx[] sIdxs = new IntIdx[len];
        
        for (int i = 0; i < len; i++) { 
            sIdxs[i] = new IntIdx(i, str.charAt(i) - '$', 0);
            indx[i] = i;
        }

        for (int k = 1; k < len * 2; k <<= 1) {
            
            for (int i = 0; i < len - k; i++) {
                sIdxs[indx[i]].num2 = sIdxs[indx[i+k]].num1;
            }
            for (int i = len - k; i < len; i++) {
                sIdxs[indx[i]].num2 = -1;
            }
            Arrays.sort(sIdxs);

            int idx = 0;
            for (int i = 0; i < len; i++) {
                indx[sIdxs[i].index] = i;
                if (i < len - 1 && !sIdxs[i].equals(sIdxs[i+1])) sIdxs[i].num1 = idx++;
                else sIdxs[i].num1 = idx;
            }
            if (idx + 1 >= len) break;

        }

        int[] sa = new int[len];
        for (int i = 0; i < len; i++) sa[indx[i]] = i;
        
        return sa;
    }

    /**
     * reverse burrows-wheeler transform
     * @param L
     * @return original string
     */
    public String reversebwt(String L) {
        L = L.toLowerCase();
        int rowi = 0;
        StringBuilder str = new StringBuilder("$");
        int[] ranks = rankBwt(L);
        HashMap<Character, Integer> F = getFirstCol(L);

        while (L.charAt(rowi) != '$') {
            char c = L.charAt(rowi);
            str.insert(0, c);
            rowi = F.get(c) + ranks[rowi];
        }

        return str.toString();
    }

    private int[] rankBwt(String L) {
        HashMap<Character, Integer> total = new HashMap<>();
        int[] ranks = new int[L.length()];
        int i = 0;
        for (char c : L.toCharArray()) {
            if (! total.containsKey(c)) total.put(c, 0);
            ranks[i++] = total.get(c);
            total.put(c, total.get(c) + 1);
        }
        return ranks;
    }

    public HashMap<Character, Integer> getFirstCol(String L) {
        HashMap<Character, Integer> F = new HashMap<>();
        for (char c : L.toCharArray()) {
            if (F.containsKey(c)) F.put(c, F.get(c) + 1);
            else F.put(c, 1);
        }
        List<Character> keyList = new ArrayList<>(F.keySet());
        keyList.sort(Comparator.naturalOrder());
        int total = 0;
        for (char c : keyList) {
            int count = F.get(c);
            F.put(c, total);
            total += count;
        }
        return F;
    }

}
