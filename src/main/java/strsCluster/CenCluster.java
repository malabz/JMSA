package strsCluster;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import measure.score;
import measure.starDist;
import msa.centerAlign;
import psa.FMAlign;


public class CenCluster {
    private final double sim;
    private String[] strs;
    private HashMap<Integer, int[]> clusters;
    
    /**
     * gen clusters
     * @param strs
     * @param sim
     */
    public CenCluster (String[] strs, double sim) {
        this.sim = sim;
        this.strs = strs;
        genClusters();
    }

    /**
     * return the idxs of each cluster
     * @return clusters
     */
    public HashMap<Integer, int[]> getClusters () {
        return clusters;
    }

    /**
     * gen clusters
     */
    private void genClusters () {
        clusters = new HashMap<>();
        // a map of change idx --> origin idx
        HashMap<Integer, Integer> idxMap = new HashMap<>();
        for (int i = 0; i < strs.length; i++) idxMap.put(i, i);

        centerAlign cAlign = new centerAlign(strs, 1);
        strs = cAlign.getStrsAlign();

        while (strs.length > 1) {
            // pick the longest one
            int idxc = pickLongest();
            // compute the distance between the longest and the others
            starDist sDist = new starDist(strs, true);
            
            // pick the similars to be a class
            int[] simOnes = pickSimilars(sDist.getDismatrix1D(idxc));
            for (int i = 0; i < simOnes.length; i++) {
                simOnes[i] = simOnes[i] >= idxc ? simOnes[i] + 1 : simOnes[i];
                // simOnes[i] = idxMap.remove(simOnes[i]);
            }
            simOnes = delOutliers(idxMap, idxc, simOnes);
            clusters.put(idxMap.remove(idxc), simOnes);
            
            Integer[] keys = idxMap.keySet().toArray(Integer[]::new);
            if (idxMap.size() == 0) break;
            else if (idxMap.size() == 1) {
                clusters.put(idxMap.remove(keys[0]), new int[0]);
                break;
            }
            String[] newStrs = new String[idxMap.size()];
            for (int i = 0; i < keys.length; i++) {
                newStrs[i] = strs[keys[i]];
                int value = idxMap.remove(keys[i]);
                // idxMap.remove(keys[i]);
                idxMap.put(i, value);
            }
            this.strs = newStrs;
        }
        if (idxMap.size() == 1) {
            Integer[] keys = idxMap.keySet().toArray(Integer[]::new);
            clusters.put(idxMap.remove(keys[0]), new int[0]);
        }
        System.out.println("clusters : " + clusters.size());
    }

    private int[] delOutliers(HashMap<Integer, Integer> idxMap, int idxc, int[] clusters) {
        if (clusters.length <= 2) { return clusters; }
        FMAlign fmAlign = new FMAlign(strs[idxc].replaceAll("-", ""));
        score sc = new score();
        double[] scores = new double[clusters.length];
        double ave = 0.0, sigma = 0.0;
        for (int i = 0; i < scores.length; i++) { 
            fmAlign.AlignStrB(strs[clusters[i]].replaceAll("-", ""));
            scores[i] = sc.sp(fmAlign.getStrAlign()[0], fmAlign.getStrAlign()[1]);
            ave += scores[i];
        }
        ave /= scores.length;
        for (double score : scores) { sigma += Math.pow(score - ave, 2); }
        sigma = Math.sqrt(sigma/scores.length);
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < scores.length; i++) {
            if (scores[i] >= ave - sigma || scores[i] > 0.85) { res.add(clusters[i]); }
        }
        int[] resInt = new int[res.size()];
        int i = 0;
        for ( int value : res) { resInt[i++] = idxMap.remove(value); }
        return resInt;
    }

    /**
     * find the longest one
     */
    private int pickLongest () {
        int res = 0;
        for (int i = 1; i < strs.length; i++) {
            res = strs[i].replaceAll("-", "").length() > strs[res].replaceAll("-", "").length() ? i : res;
        }
        return res;
    }

    /**
     * pick similar ones
     * @param dismatrix
     * @return idxs
     */
    private int[] pickSimilars (double[] dismatrix) {
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < dismatrix.length; i++) {
            if (dismatrix[i] >= sim) {
                res.add(i);
            }
        }
        int i = 0;
        int[] resInt = new int[res.size()];
        for (Integer r : res) resInt[i++] = r;
        return resInt;
    }

}
