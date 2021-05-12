package measure;

import java.util.HashMap;

// import cluster.kmer;

public class score {
    
    /**
     * compute the sps score in [0, 1]
     * @param strs
     * @return score
     */
    public double sps(String[] strs) {
        int nums = strs.length, len = strs[0].length();
        long match = 0, gap = 0, all;
        all = (long) nums * (nums - 1) * len / 2;
        for (int i = 0; i < nums; i++) {
            String A = strs[i];
            for (int j = i + 1; j < nums; j++) {
                String B = strs[j];
                for (int k = 0; k < len; k++) {
                    char a = A.charAt(k), b = B.charAt(k);
                    if (a == b) { 
                        if (a != '-') match++;
                        else gap++; 
                    }

                }
            }
        }
        return (double) match/(all-gap);
    }


    public double tc(String[] strs) {
        int nums = strs.length;
        int len = strs[0].length();
        int match = 0;
        for (int i = 0; i < len; i++) {
            String scren = (i + 1) + "/" + len;
            System.out.print(scren);
            HashMap<Character, Integer> charIdx = new HashMap<>();
            for (String str : strs) {
                char c = str.charAt(i);
                if (charIdx.containsKey(c)) charIdx.put(c, charIdx.get(c) + 1);
                else charIdx.put(c, 1);
            }
            if ( charIdx.containsKey('-') && charIdx.get('-') / nums > 0.6) {
                System.out.print("\b".repeat(scren.length()));
                continue;
            }
            charIdx.remove('-');
            if (charIdx.size() == 1) { 
                match++; 
                System.out.print("\b".repeat(scren.length()));
                continue; 
            }
            int all = 0, max = 0;
            for (int value : charIdx.values()) {
                max = value > max ? value : max;
                all += value;
            }
            if ((double) max / all > 0.8) match++;
            System.out.print("\b".repeat(scren.length()));
        }
        System.out.println();
        return (double) match / len;
    }
}
