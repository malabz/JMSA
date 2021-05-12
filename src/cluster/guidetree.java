package cluster;

import java.io.IOException;

import io.fasta;

public class guidetree {
    private String[] strs;
    private String[] labels;
    private String mode;
    private int k;
    
    public static void main(String[] args) throws IOException {
        String path = "/home/siyang/Documents/MyWork/msa/code/data/origin/"; 
        String data = "steps-1024_0.fasta";

        fasta readfasta = new fasta();
        String[][] res = readfasta.readFasta(path+data);
        String[] strs = res[1];

        kmer km = new kmer(strs);
        double[][] ma = km.getDismatrix();
        long startTime = System.currentTimeMillis();
        upgma htree = new upgma(ma);
        System.out.println(htree.genTree());
        long endTime = System.currentTimeMillis();
        System.out.println("upgma cost time: "+((endTime-startTime)/1000)+"s");
        startTime = System.currentTimeMillis();
        NeighborJoining nj = new NeighborJoining(ma);
        System.out.println(nj.genTree());
        endTime = System.currentTimeMillis();
        System.out.println("nj cost time: "+((endTime-startTime)/1000)+"s");
    }
    /**
     * mode: "nj" or "upgma"
     * @param strs
     * @param labels
     * @param mode
     */
    public guidetree(String[] strs, String[] labels, String mode) {
        this.strs = strs;
        this.labels = labels;
        this.mode = mode;
        this.k = 4;
    }

    /**
     * mode: "nj" or "upgma"
     * @param strs
     * @param mode
     */
    public guidetree(String[] strs, String mode) {
        this.strs = strs;
        this.mode = mode;
        this.k = 4;
    }

    public guidetree(String[] strs) {
        this.strs = strs;
        this.mode = "upgma";
        this.k = 4;
    }

    /**
     * 
     * @return root node
     */
    public node genTree() {
        kmer km = new kmer(strs, k);
        double[][] simmatrix = km.getDismatrix();
        upgma htree = new upgma(simmatrix, labels);
        return htree.genTree();
    }

    public int[][] genTreeList() {
        kmer km = new kmer(strs, k);
        double[][] simMatrix = km.getDismatrix();
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
}
