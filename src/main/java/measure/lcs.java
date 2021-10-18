package measure;


import psa.STAlign;

/**
 * 计算序列间的最长公共子序列的距离
 */
public class lcs {
    String[] strs;

    public lcs (String[] strs) {
        this.strs = strs;
    }

    /**
     * 计算序列间两两之间的距离
     */
    public double[][] getDismatrix() {
        double[][] dismatrix = new double[strs.length][strs.length];
        for (int i = 0; i < strs.length; i++) {
            STAlign sAlign = new STAlign(strs[i]);
            for (int j = i + 1; j < strs.length; j++) {
                dismatrix[i][j] = dismatrix[j][i] = sAlign.getSimstrB(strs[j]);
            }
        }
        return dismatrix;
    }

    /**
     * 计算中心序列与其他序列之间的距离
     */
    public double[] getDismatrix1D(int idxc) {
        double[] dismatrix = new double[strs.length - 1];
        STAlign sAlign = new STAlign(strs[idxc]);
        for (int i = 0; i < idxc; i++) {
            dismatrix[i] = sAlign.getSimstrB(strs[i]);
        }
        for (int i = idxc + 1; i < strs.length; i++) {
            dismatrix[i - 1] = sAlign.getSimstrB(strs[i]);
        }
        return dismatrix;
    }
}