package sample;

import java.util.ArrayList;
import java.util.List;

public class sampleStrings {
    // sample 2/rate%
    private final int rate = 20;

    public String[] getSampleStrs(String[] strs) {
        int len = rate;
        if (strs.length <= rate / 2) return strs;
        else if (strs.length < rate * rate / 4) len = strs.length * 4 / rate;
        else if (strs.length >= 1000) len = strs.length / 50;
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i <= strs.length / len; i++) {
            int[] temp = LongShortest(strs, i * len, (i + 1) * len);
            if (temp == null) break;
            if (temp[0] == temp[1]) {
                res.add(temp[0]);
            }
            else {
                res.add(temp[0]);
                res.add(temp[1]);
            } 
        }
        String[] sampledStrs = new String[res.size()];
        int j = 0;
        for (int i : res) sampledStrs[j++] = strs[i];
        return sampledStrs;
    }

    private int[] LongShortest(String[] strs, int start, int end) {
        if (start >= strs.length) return null;
        end = Math.min(end, strs.length);
        int longOne = start, shortOne = end - 1;
        for (int j = start; j < end; j++) {
            if (strs[j].length() > strs[longOne].length()) { longOne = j; }
            else if (strs[j].length() <= strs[shortOne].length()) { shortOne = j; }
        }
        return new int[]{longOne, shortOne};
    } 

}
