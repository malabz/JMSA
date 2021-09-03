package measure;


public class strsdist {
    private String[][] strings;
    private String[] strs;
    private String mode;
    public int idxc = -1;
    private double[][] dismatrix;
    private double[] dismatrix1D;

    /**
     * get 1D distance matrix
     * this used for cluster to choose the most similar
     * @param strs
     * @param mode "kmer" "lcs" "star"
     * @param idxc
     */
    public strsdist (String[] strs, String mode, int idxc) {
        this.strs = strs;
        this.mode = mode.toLowerCase();
        this.idxc = idxc;
        dismatrix1D = getDist1D();
    }

    /**
     * get 2D distance matrix
     * @param strs
     * @param mode "kmer" "lcs"
     */
    public strsdist (String[] strs, String mode) {
        this.strs = strs;
        this.mode = mode.toLowerCase();
        dismatrix = getDist2D();
    }

    /**
     * 用于计算profile之间的距离
     * @param strings
     */
    public strsdist (String[][] strings) {
        this.strings = strings;
        this.mode = "kmer";
        dismatrix = getDistND();
    }

    /**
     * get 1D dist matrix,
     * mode is kmer
     * @param strs
     * @param idxc
     */
    public strsdist (String[] strs, int idxc) {
        this.strs = strs;
        this.mode = "kmer";
        this.idxc = idxc;
        dismatrix1D = getDist1D();
    }

    /**
     * get the distance2D
     * @return dist[][]
     */
    public double[][] getDismatrix2D () {
        if (idxc == -1) return dismatrix;
        else return null;
    }

    /**
     * get the distance1D
     * @return dist[]
     */
    public double[] getDismatrix1D () {
        if (idxc != -1) return dismatrix1D;
        else return null;
    }

    private double[][] getDistND () {
        kmer km = new kmer(strings);
        return km.getDismatrixND();
    }

    /**
     * 计算序列间两两之间的距离
     */
    private double[][] getDist2D () {
        switch (mode) {
            case "kmer" -> {
                kmer ker = new kmer(strs);
                return ker.getDismatrix();
            }
            case "lcs" -> {
                lcs ls = new lcs(strs);
                return ls.getDismatrix();
            }
            default -> {
                System.out.println("(strsdist) mode is wrong");
                System.exit(0);
            }
        }
        return null;
    }

    /**
     * 计算中心序列与其他序列之间的距离
     */
    private double[] getDist1D () {
        switch (mode) {
            case "kmer" -> {
                kmer ker = new kmer(strs);
                return ker.getDismatrix1D(idxc);
            }
            case "lcs" -> {
                lcs ls = new lcs(strs);
                return ls.getDismatrix1D(idxc);
            }
            case "star" -> {
                starDist sdist = new starDist(strs, false);
                return sdist.getDismatrix1D(idxc);
            }
            default -> {
                System.out.println("mode is wrong");
                System.exit(0);
            }
        }
        return null;
    }
}
