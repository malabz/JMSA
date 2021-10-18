package psa;

import java.util.ArrayList;
import java.util.List;

public class FastMultiAlign {
    private final String[] A, B;
    private final int idxcA, idxcB;
    private String[] alignA, alignB;
    private final char[] alphabet;
    private boolean state = false;

    public FastMultiAlign(String[] A, String[] B, char[] alphabet, int kk, int idxcA, int idxcB) {
        if (A[0].length() >= B[0].length()) {
            this.A = A;
            this.B = B;
            this.idxcA = idxcA;
            this.idxcB = idxcB;
        }
        else {
            this.A = B;
            this.B = A;
            this.idxcA = idxcB;
            this.idxcB = idxcA;
            state = true;
        }
        this.alphabet = alphabet;
        this.Align(kk);
    }
    
    /**
     * To get the aligned results.
     */
    public String[][] getStrAlign() {
        if (state) { return new String[][]{alignB, alignA}; }
        return new String[][]{alignA, alignB};
    }

    public String[] getStrsAlign() {
        String[] res = new String[alignA.length + alignB.length];
        if (state) {
            System.arraycopy(alignB, 0, res, 0, alignB.length);
            System.arraycopy(alignA, 0, res, alignB.length, alignA.length);
        }
        else {
            System.arraycopy(alignA, 0, res, 0, alignA.length);
            System.arraycopy(alignB, 0, res, alignA.length, alignB.length);
        }
        return res;
    }
    
    private void Align(int kk) {
        // 1. get Homo-areas
        int[][] regions;
        regions = getNewRegions(idxcA, idxcB);
        
        // 2. evaluate the homo areas
        if (regions.length == 0) {
            multiKband mkband = new multiKband(A, B, alphabet, kk);
            alignA = mkband.getStrAlign()[0];
            alignB = mkband.getStrAlign()[1];
            return;
        }
        int oldA = 0, oldB = 0;
        StringBuilder[] sbA = new StringBuilder[A.length];
        StringBuilder[] sbB = new StringBuilder[B.length];
        for (int i = 0; i < sbA.length; i++) { sbA[i] = new StringBuilder(); }
        for (int i = 0; i < sbB.length; i++) { sbB[i] = new StringBuilder(); }
        for (int[] region : regions) {
            String[] tempA = getSubRegion(oldA, region[2], A);
            String[] tempB = getSubRegion(oldB, region[0], B);
            oldA = region[2] + region[1];
            oldB = region[0] + region[1];
            multiKband mkband = new multiKband(tempA, tempB, alphabet);
            for (int j = 0; j < sbA.length; j++) {
                sbA[j].append(mkband.getStrAlign()[0][j]);
            }
            for (int j = 0; j < sbB.length; j++) {
                sbB[j].append(mkband.getStrAlign()[1][j]);
            }
            for (int j = 0; j < sbB.length; j++) {
                sbB[j].append(B[j], region[0], oldB);
            }
            for (int j = 0; j < sbA.length; j++) {
                sbA[j].append(A[j], region[2], oldA);
            }
        }
        String[] tempA = getSubRegion(oldA, A[0].length(), A);
        String[] tempB = getSubRegion(oldB, B[0].length(), B);
        multiKband mkband = new multiKband(tempA, tempB, alphabet);
        for (int j = 0; j < sbA.length; j++) {
            sbA[j].append(mkband.getStrAlign()[0][j]);
        }
        for (int j = 0; j < sbB.length; j++) { 
            sbB[j].append(mkband.getStrAlign()[1][j]);
        }
        alignA = new String[A.length];
        alignB = new String[B.length];
        for (int j = 0; j < sbA.length; j++) {
            alignA[j] = sbA[j].toString();
        }
        for (int j = 0; j < sbB.length; j++) {
            alignB[j] = sbB[j].toString();
        }
    }

    private String[] getSubRegion(int s, int t, String[] strs) {
        String[] res = new String[strs.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = strs[i].substring(s, t);
        }
        return res;
    }

    private int[][] getNewRegions(int lgA, int lgB) {
        String ctA = A[lgA].replaceAll("-", "");
        String ctB = B[lgB].replaceAll("-", "");
        int[][] regions;
        FMAlign fmAlign = new FMAlign(ctA);    
        regions = fmAlign.getSubStrings(ctB);
        if (regions == null) return new int[0][0];
        // return regions;
        String orA = A[lgA];
        String orB = B[lgB];
        int[] mapA = new int[ctA.length()];
        int[] mapB = new int[ctB.length()];
        int j = 0;
        for (int i = 0; i < orA.length(); i++) {
            if (orA.charAt(i) == '-') continue;
            mapA[j++] = i;
            if (j == ctA.length()) break;
        }
        j = 0;
        for (int i = 0; i < orB.length(); i++) {
            if (orB.charAt(i) == '-') continue;
            mapB[j++] = i;
            if (j == ctB.length()) break;
        }
        List<int[]> regionsSplit = new ArrayList<>();
        for (int[] idx : regions) {
            boolean sta = true;
            boolean stb = true;
            if (mapA[idx[2]] + idx[1] == mapA[idx[2] + idx[1] - 1] + 1) sta = false;
            if (mapB[idx[0]] + idx[1] == mapB[idx[0] + idx[1] - 1] + 1) stb = false;

            if ((!sta) && (!stb)) regionsSplit.add(idx);
            else if (sta && stb) {
                int[] tempA = getSplits(idx[2], idx[2] + idx[1], mapA);
                int[] tempB = getSplits(idx[0], idx[0] + idx[1], mapB);
                int i = 0, oldValueA = tempA[0], oldValueB = tempB[0];
                j = 0;                
                while (i < tempA.length || j < tempB.length) {
                    if ((j < tempB.length) && (i == tempA.length || tempA[i] - tempA[0] > tempB[j] - tempB[0])) {
                        regionsSplit.add(new int[]{oldValueB, tempB[j] - oldValueB, oldValueA});
                        oldValueA += (tempB[j] - oldValueB);
                        oldValueB = tempB[j++];
                    }
                    else if (j == tempB.length || tempA[i] - tempA[0] < tempB[j] - tempB[0]) {
                        regionsSplit.add(new int[]{oldValueB, tempA[i] - oldValueA, oldValueA});
                        oldValueB += (tempA[i] - oldValueA);
                        oldValueA = tempA[i++];
                    }
                    else {
                        if (i != 0 && j != 0) { 
                            regionsSplit.add(new int[]{oldValueB, tempB[j] - oldValueB, oldValueA});
                        }
                        oldValueA = tempA[i++];
                        oldValueB = tempB[j++];
                    }
                }
                regionsSplit.add(new int[]{oldValueB, idx[1] - oldValueA + idx[2], oldValueA});
            }
            else if (sta) {
                int[] tempA = getSplits(idx[2], idx[2] + idx[1], mapA);
                int oldValueB = idx[0];
                for (int i = 1; i < tempA.length; i++) {
                    regionsSplit.add(new int[]{oldValueB, tempA[i] - tempA[i-1], tempA[i-1]});
                    oldValueB += (tempA[i] - tempA[i-1]);
                }
                regionsSplit.add(new int[]{oldValueB, idx[0] + idx[1] - oldValueB, tempA[tempA.length-1]});
            }
            else {
                int[] tempB = getSplits(idx[0], idx[0] + idx[1], mapB);
                int oldValueA = idx[2];
                for (int i = 1; i < tempB.length; i++) {
                    regionsSplit.add(new int[]{tempB[i-1], tempB[i] - tempB[i-1], oldValueA});
                    oldValueA += (tempB[i] - tempB[i-1]);
                }
                regionsSplit.add(new int[]{tempB[tempB.length-1], idx[2] + idx[1] - oldValueA, oldValueA});
            }
        }
        
        int[][] res = select(regionsSplit.toArray(new int[regionsSplit.size()][]), mapA, mapB, A[lgA].length(), B[lgB].length());
        res = eval(res, A[lgA].length(), B[lgB].length());

        for (int i = 0; i < res.length; i++) {
            res[i][0] = mapB[res[i][0]];
            res[i][2] = mapA[res[i][2]];
        }
        return res;
    }

    private int[][] select(int[][] regions, int[] mapA, int[] mapB, int lenA, int lenB) {
        List<int[]> newValues = new ArrayList<>();
        double tolerance = 1.05 * Math.abs(lenA - lenB) / Math.max(lenB, lenA);
        for (int[] ints : regions) {
            if (ints[1] >= 15) {
                double off1 = (double) mapB[ints[0]] / lenB;
                double off2 = (double) mapA[ints[2]] / lenA;
                double diff = Math.abs(off1 - off2);
                if (diff < tolerance) newValues.add(ints);
            }
        }
        return newValues.toArray(new int[newValues.size()][]);
    }

    private int[][] eval(int[][] regions, int lenA, int lenB) {
        if (regions.length == 0) { return regions; }
        else if (regions.length == 1) { 
            int forward = Math.abs(regions[0][2] - regions[0][0]);
            int backward = Math.abs((lenA - regions[0][2]) - (lenB - regions[0][0]));
            if (forward >= regions[0][1] && backward >= regions[0][1]) return new int[0][];
            return regions;
        }
        int idx = 0;
        // for first one
        while (idx < regions.length) {
            int ave = 0;
            for (int[] region : regions) ave += (region[2]-region[0]);
            ave /= regions.length;
            int forward = regions[idx][2] - regions[idx][0];
            if (forward > ave && Math.abs(forward) >= regions[0][1]) regions[idx++] = new int[3];
            else break;
        }
        if (idx >= regions.length) return new int[0][];

        // for other one
        List<int[]> res = new ArrayList<>();
        res.add(regions[idx]);
        for (int i = idx+1; i < regions.length; i++) {
            int mid = Math.abs(regions[i][2]-regions[idx][2]-regions[i][0]+regions[idx][0]);
            if (mid > regions[i][1]) regions[i][1] = -1;
            else {
                idx = i;
                res.add(regions[idx]);
            }
        }
        return res.toArray(new int[res.size()][]);
    }

    /**
     * start, end - 1
     * @param start
     * @param end
     * @param idx
     */
    private int[] getSplits(int start, int end, int[] idx) {
        List<Integer> res = new ArrayList<>();
        res.add(start);
        for (int i = start + 1; i < end; i++) {
            if (idx[i] != idx[i-1] + 1) res.add(i);
        }
        int[] resInt = new int[res.size()];
        int i = 0;
        for (int j : res) resInt[i++] = j;
        return resInt;
    }
}
