package measure;

import msa.centerAlign;


/**
 * 计算序列间与最长子序列的距离
 */
public class starDist {
    private String[] strs;
    private final boolean aligned;

    /**
     * 输入序列，对齐或者未对齐
     * @param strs
     * @param aligned
     */
    public starDist(String[] strs, boolean aligned) {
        this.strs = strs;
        this.aligned = aligned;
    }

    /**
     * 得到中心序列与其他序列之间的距离
     * @param idxc
     * @return distance
     */
    public double[] getDismatrix1D (int idxc) {
        double[] dismatrix = new double[strs.length - 1];
        if (!aligned) {
            centerAlign cAlign = new centerAlign(strs, 1);
            strs = cAlign.getStrsAlign();
        }
        for (int i = 0; i < idxc; i++) {
            dismatrix[i] = getDist(strs[idxc], strs[i]);
        }
        for (int i = idxc + 1; i < strs.length; i++) {
            dismatrix[i - 1] = getDist(strs[idxc], strs[i]);
        }
        return dismatrix;
    }

    private double getDist (String str1, String str2) {
        int match = 0, gap = 0, len = str1.length();
        assert (str1.length() == str2.length());
        for (int i = 0; i < str1.length(); i++) {
            if (str1.charAt(i) == str2.charAt(i)) {
                if (str1.charAt(i) != '-') match++;
                else gap++;
            }
        }
        return (double) match / (len - gap);
    }

    public double[][] getDismatrix2D() {
        double[][] dismatrix = new double[strs.length][strs.length];
        if (!aligned) {
            centerAlign cAlign = new centerAlign(strs, 1);
            strs = cAlign.getStrsAlign();
        }
        for (int i = 0; i < strs.length; i++) {
            for (int j = 0; j < i; j++) {
                dismatrix[i][j] = dismatrix[j][i] = getDist(strs[i], strs[j]);
            }
        }
        return dismatrix;
    }
}
