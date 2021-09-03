package measure;

import java.util.HashMap;

import msa.centerAlign;
import sample.sampleStrings;

public class score {

    private final int ms = 7, mis = -3;
    private final int d = 13, e = 2;

    
    public double sp(String A, String B) {
        int gap = 0;
        boolean stat = false;
        long score = 0;
        assert A.length() == B.length();
        for (int i = 0; i < A.length(); i++) {
            if (A.charAt(i) == B.charAt(i)) {
                if (A.charAt(i) == '-') gap++;
                else { score += ms; stat = false;}
            }
            else if(A.charAt(i) == '-' || B.charAt(i) == '-') {
                score -= stat ? e : d;
                stat = true;
            }
            else {
                score += mis;
                stat = false;
            }
        }
        return (double) score/((A.length() - gap) * ms);
    }

    /**
     * used to get the 2D distance matrix
     * @param strsed
     */
    public double[][] getDist(String[] strsed) {
        double[][] dis = new double[strsed.length][strsed.length];
        for (int i = 0; i < strsed.length; i++) {
            for (int j = i + 1; j < strsed.length; j++) {
                dis[i][j] = dis[j][i] = sp(strsed[i], strsed[j]);
            }
        }
        return dis;
    }
    
    /**
     * compute the sps score in [0, 1]
     * @param strs
     * @return score
     */
    public double sps(String[] strs) {
        long nums = strs.length, len = strs[0].length();
        double match = 0, gap = 0;
        for (int i = 0; i < len; i++) {
            long tempmatch = 0, tempgap = 0;
            String scren = (i + 1) + "/" + len;
            System.out.print(scren);
            HashMap<Character, Long> al = new HashMap<>();
            for (String str : strs) {
                char c = str.charAt(i);
                al.put(c, al.containsKey(c) ? al.get(c) + 1 : 1);
            }
            long numg = al.containsKey('-') ? al.get('-') : 0; 
            long numn = al.containsKey('n') ? al.get('n') : 0;
            for (char c : al.keySet()) {
                if (c != '-') tempmatch += (al.get(c) * (al.get(c) - 1) / 2);
            }
            tempmatch += (numn * (nums - numg - numn));
            tempgap += (numg * (numg - 1) / 2);
            match += ((double) tempmatch / (double) (nums * (nums - 1) / 2));
            gap += ((double) tempgap / (double) (nums * (nums - 1) / 2));
            System.out.print("\b".repeat(scren.length()));
        }
        System.out.print("match:" + String.format("%.2f", match) + " ");
        System.out.print("gap:" + String.format("%.2f", gap) + " ");
        System.out.print("len:" + len + "\n");
        return match /(len - gap);
    }


    /**
     * compute the tc score in [0,1]
     * @param strs
     * @return
     */
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
            if ( charIdx.containsKey('-') && ((double) charIdx.get('-') / nums) > 0.6) {
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
                max = Math.max(value, max);
                all += value;
            }
            if ((double) max / all > 0.8) match++;
            System.out.print("\b".repeat(scren.length()));
        }
        return (double) match / len;
    }


    /**
     * try to get a appropriate k
     * @param strs
     * @param sampled
     * @return k
     */
    public int getK(String[] strs, boolean sampled) {
        if (!sampled) {
            sampleStrings spStrs = new sampleStrings();
            strs = spStrs.getSampleStrs(strs);    
        }
        centerAlign cAlign = new centerAlign(strs, 1);
        strs = cAlign.getStrsAlign();
        int nums = strs.length;
        // int len = strs[0].length();
        // int gap1 = 0, gap2 = 0;
        int gap2 = 0;
        // method 1: compute the whole gap
        // for (int i = 0; i < len; i++) {
        //     int numsGap = 0;
        //     for (String str : strs) {
        //         if (str.charAt(i) == '-') numsGap++;
        //     }
        //     if ((float) numsGap / nums > 0.5) gap1++;
        // }
        // method 2: compute the longest string's gap
        // int longRow = cAlign.getLongestRow();
        // gap += strs[longRow].length() -  strs[longRow].replaceAll("-", "").length();
        // method 3: 
        for (int i = 0; i < nums; i++) {
            for (int j = i; j < nums; j++) {
                // gap2 = Math.max(countGap(strs[i], strs[j]), gap2);
                gap2 += countGap(strs[i], strs[j]);
            }
        }
        gap2 /= (nums*nums/2);
        return Math.max(gap2, 1);
    }

    private int countGap(String A, String B) {
        assert A.length() == B.length();
        int nums = 0, len = A.length();
        for (int i = 0; i < len; i++) {
            if (A.charAt(i) == B.charAt(i) && A.charAt(i) == '-') nums++;
        }
        int lenA = A.replaceAll("-", "").length();
        int lenB = B.replaceAll("-", "").length();
        return Math.min(len - lenA - nums, len - lenB - nums);
    }
}
