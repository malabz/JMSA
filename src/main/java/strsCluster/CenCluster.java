package strsCluster;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import io.str;
import measure.score;
import measure.starDist;
import msa.centerAlign;
import psa.FMAlign;


public class CenCluster {
    private final double sim;
    private String[] strs;
    private int[] lens;
    private HashMap<Integer, int[]> clusters;
    
    /**
     * gen clusters
     * @param strs
     * @param sim
     */
    public CenCluster (String[] strs, double sim, boolean iter) {
        this.sim = sim;
        this.strs = strs;
        this.lens = getLens(strs);
        genClusters(iter);
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
     * gen clusters
     */
    private void genClusters(boolean iter) {
        clusters = new HashMap<>();
        // a map of change idx --> origin idx
        int[] idxMap = new int[strs.length];
        for (int i = 0; i < strs.length; i++) idxMap[i] = i;

        // method 2
        if (!iter) {
            centerAlign cAlign = new centerAlign(strs, 1);
            strs = cAlign.getStrsAlign();
        }
        
        while (strs.length > 1) {
            String output = "  clu : " + clusters.size() + "  remain : " + strs.length + "      ";
            System.out.print(output);
            // pick the longest one
            int idxc = pickLongest();
            // compute the distance between the longest and the others
            starDist sDist = new starDist(strs, !iter);
            // pick the similars to be a class
            int[] simOnes = pickSimilars(sDist.getDismatrix1D(idxc));
            for (int i = 0; i < simOnes.length; i++) {
                simOnes[i] = simOnes[i] >= idxc ? simOnes[i] + 1 : simOnes[i];
                // int temp = simOnes[i];
                // simOnes[i] = idxMap[temp];
                // idxMap[temp] = -1;
            }
            simOnes = delOutliers(idxMap, idxc, simOnes);
            clusters.put(idxMap[idxc], simOnes);
            idxMap[idxc] = -1;

            if (idxMap.length - simOnes.length == 1) {
                System.out.print(str.repeat("\b", output.length()));
                break;
            }
            else if (idxMap.length - simOnes.length == 0) {
                clusters.put(idxMap[0], new int[0]);
                System.out.print(str.repeat("\b", output.length()));
                break;
            }
            String[] newStrs = new String[idxMap.length-simOnes.length-1];
            int[] newlens = new int[idxMap.length-simOnes.length-1];
            int[] newIdx = new int[idxMap.length-simOnes.length-1];
            int j = 0;
            for (int i = 0; i < idxMap.length; i++) {
                if (idxMap[i] != -1) {
                    newStrs[j] = strs[i];
                    newlens[j] = lens[i];
                    newIdx[j++] = idxMap[i];
                }
            }
            this.strs = newStrs;
            this.lens = newlens;
            idxMap = newIdx;
            System.out.print(str.repeat("\b", output.length()));
        }
        if (idxMap.length == 1) {
            clusters.put(idxMap[0], new int[0]);
        }
        // System.out.println("  clusters : " + clusters.size() + "                    ");
    }

    private int[] delOutliers(int[] idxMap, int idxc, int[] clusters) {
        if (clusters.length <= 2) { 
            int[] resInt = new int[clusters.length];
            for (int i = 0; i < clusters.length; i++) {
                resInt[i] = idxMap[clusters[i]];
                idxMap[clusters[i]] = -1;
            }
            return resInt; 
        }
        FMAlign fmAlign = new FMAlign(strs[idxc].replaceAll("-", ""));
        double[] scores = new double[clusters.length];
        double ave = 0.0;
        // double sigma = 0.0;
        for (int i = 0; i < scores.length; i++) { 
            fmAlign.AlignStrB(strs[clusters[i]].replaceAll("-", ""));
            scores[i] = score.sp(fmAlign.getStrAlign()[0], fmAlign.getStrAlign()[1]);
            ave += scores[i];
        }
        ave /= scores.length;
        // for (double score : scores) { sigma += Math.pow(score - ave, 2); }
        // sigma = Math.sqrt(sigma/scores.length);
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < scores.length; i++) {
            if (scores[i] >= ave || scores[i] >= 0.9) { res.add(clusters[i]); }
        }
        int[] resInt = new int[res.size()];
        int i = 0;
        for ( int value : res) { 
            resInt[i++] = idxMap[value]; 
            idxMap[value] = -1;
        }
        return resInt;
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
