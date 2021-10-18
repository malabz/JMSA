package strsCluster;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

// import measure.strsdist;

public class FastCluster {
    private final double sim;
    private final String[] strs;
    private final int[] lens;
    private HashMap<Integer, int[]> clusters;
    
    public FastCluster(String[] strs, double sim) {
        this.sim = sim;
        this.strs = strs;
        this.lens = getLens(strs);
        this.genClusters();
    }

    /**
     * return the idxs of each cluster
     * @return clusters
     */
    public HashMap<Integer, int[]> getClusters () {
        return clusters;
    }

    private int[] getLens(String[] strs) {
        int[] res = new int[strs.length];
        for (int i = 0; i < strs.length; i++) { res[i] = strs[i].length(); }
        return res;
    }


    /**
     * find the longest one
     */
    private int pickLongest() {
        int res = 0;
        for (int i = 1; i < lens.length; i++) {
            res = lens[i] > lens[res] ? i : res;
        }
        return res;
    }

    private int[] LensToFind(int idxc) {
        int length = (int) (lens[idxc] * 0.90);
        lens[idxc] = -1;
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < lens.length; i++) {
            if (lens[i] >= length) {
                res.add(i);
                lens[i] = -1;
            }
        }
        int[] resInt = new int[res.size()+1];
        int i = 1;
        resInt[0] = idxc;
        for (int r : res) resInt[i++] = r;
        // res = new ArrayList<>();

        // String[] strings = new String[resInt.length];
        // for (i = 0; i < strings.length; i++) strings[i] = strs[resInt[i]];

        // strsdist sdist = new strsdist(strings, "kmer", 0);
        // double[] sims = sdist.getDismatrix1D();
        
        // for (i = 0; i < sims.length; i++) {
        //     if (sims[i] >= 0.6) {
        //         res.add(resInt[i+1]);
        //         lens[resInt[i+1]] = -1;
        //     }
        // }
        // resInt = new int[res.size()+1];
        // i = 1;
        // resInt[0] = idxc;
        // for (int r : res) resInt[i++] = r;
        
        return resInt;
    }

    private void genClusters() {
        int remainder = this.strs.length;
        clusters = new HashMap<>();
        HashMap<Integer, int[]> LensClusters = new HashMap<>();
        while (remainder > 0) {
            int idxc = pickLongest();
            LensClusters.put(idxc, LensToFind(idxc));
            remainder -= (LensClusters.get(idxc).length);
        }
        for (int key : LensClusters.keySet()) {
            int[] idxs = LensClusters.get(key);
            String[] strings = new String[idxs.length];
            for (int i = 0; i < strings.length; i++) { strings[i] = strs[idxs[i]]; }
            CenCluster ccluster = new CenCluster(strings, sim, false);
            HashMap<Integer, int[]> temp = ccluster.getClusters();
            for (int k : temp.keySet()) {
                int[] ints = temp.get(k);
                for (int i = 0; i < ints.length; i++) ints[i] = idxs[ints[i]];
                clusters.put(idxs[k], ints);
            }
        }
        System.out.println("  clusters : " + clusters.size() + "                    ");
    }
}
