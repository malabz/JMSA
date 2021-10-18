package hierCluster;

import measure.score;
import measure.strsdist;

public class guidetree {
    private String[] strs;
    private String[] strsed;
    private final String mode;

    /**
     * mode: "nj" or "upgma"
     */
    public guidetree(String[] strs, String mode) {
        this.strs = strs;
        this.mode = mode;
    }

    /**
     * 用于后续比对的生成树
     */
    public guidetree(String[] strsed, int silent, String mode) {
        this.strsed = strsed;
        this.mode = mode;
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
            simMatrix = score.getDist(strsed);
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
            throw new IllegalArgumentException("unkown mode: " + mode);
        }
    }
}
