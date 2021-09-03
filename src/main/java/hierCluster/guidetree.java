package hierCluster;

import measure.score;
// import measure.starDist;
import measure.strsdist;

public class guidetree {
    private String[][] strings;
    private String[] strs;
    private String[] strsed;
    private final String mode;

    /**
     * mode: "nj" or "upgma"
     * @param strs
     * @param mode
     */
    public guidetree(String[] strs, String mode) {
        this.strs = strs;
        this.mode = mode;
    }

    /**
     * 用于后续比对的生成树
     * @param strsed
     */
    public guidetree(String[] strsed, int nonem, String mode) {
        this.strsed = strsed;
        this.mode = mode;
    }

    /**
     * 用于聚类过后的profile生成树
     * @param strings
     * @param mode
     */
    public guidetree(String[][] strings, String mode) {
        this.strings = strings;
        this.mode = mode;
    }


    /**
     * 用于profile的生成树
     */
    public int[][] genTreeListND() {
        strsdist sdist = new strsdist(strings);
        double[][] simMatrix = sdist.getDismatrix2D();
        if (mode.equals("upgma")) {
            upgma htree = new upgma(simMatrix);
            htree.genTree();
            return htree.TreeList.toArray(new int[htree.TreeList.size()][3]);
        }
        else if (mode.equals("nj")) {
            NeighborJoining nJtree = new NeighborJoining(simMatrix);
            nJtree.genTree();
            return nJtree.TreeList.toArray(new int[nJtree.TreeList.size()][3]);
        }
        else {
            System.out.println("The mode that gens the guidetree is wrong! (two options: nj or upgma)");
            System.exit(0);
        }
        return null;
    }

    /**
     * gen a order list of the treealign 
     */
    public int[][] genTreeList(int silent) {
        double[][] simMatrix;
        // compute similarity matrix
        if (strsed == null) {
            strsdist sdist = new strsdist(strs, "kmer");
            // starDist sdist = new starDist(strs, false);
            simMatrix = sdist.getDismatrix2D();
        }
        else {
            score sc = new score();
            simMatrix = sc.getDist(strsed);
        }
        if (mode.equals("upgma")) {
            upgma htree = new upgma(simMatrix);
            htree.genTree();
            return htree.TreeList.toArray(new int[htree.TreeList.size()][3]);
        }
        else if (mode.equals("nj")) {
            NeighborJoining nJtree = new NeighborJoining(simMatrix);
            if (silent == 1) nJtree.genTreeSilent();
            else nJtree.genTree();
            return nJtree.TreeList.toArray(new int[nJtree.TreeList.size()][3]);
        }
        else {
            System.out.println("The mode that gens the guidetree is wrong! (two options: nj or upgma)");
            System.exit(0);
        }
        return null;
    }
}
