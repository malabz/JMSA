package cluster;

import java.util.ArrayList;
import java.util.List;

public class NeighborJoining {
    private double[][] dmatrix, qmatrix;
    private double[] numsMatrix;
    public List<int[]> TreeList;
    private node[] nodes;
    private int n;
    private int global_n;

    /**
     * 
     * @param matrix
     * @param names
     */
    public NeighborJoining(double[][] matrix, String[] names) {
        this.dmatrix = matrix;
        this.n = names.length;
        this.nodes = new node[n];
        for (int i = 0; i < n; i++) {
            node temp = new leafnode(names[i], i);
            this.nodes[i] = temp;
        }
        this.global_n = this.n;
        this.TreeList = new ArrayList<>();
    }

    public NeighborJoining(double[][] matrix) {
        this.dmatrix = matrix;
        this.n = matrix.length;
        this.nodes = new node[this.n];
        String[] names = new String[this.n];
        for (int i = 0; i < n; i++) {
            names[i] = ""+i;
            node temp = new leafnode(names[i], i);
            this.nodes[i] = temp;
        }
        this.global_n = this.n;
        this.TreeList = new ArrayList<>();
    }

    // gen the Q-matrix
    private void genQcriterion() {
        numsMatrix = new double[this.n];
        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                numsMatrix[i] += dmatrix[i][j];
            }
        }
        qmatrix = new double[this.n][this.n];
        for (int i = 0; i < this.n; i++) {
            for (int j = i+1; j < this.n; j++) {
                qmatrix[i][j] = dmatrix[i][j] - (numsMatrix[i] + numsMatrix[j]) / (this.n - 2);
                qmatrix[j][i] = qmatrix[i][j];
            }
        }
    }

    /**
     * 
     * @return unrooted tree
     */
    public midnode genTree() {
        if (this.n < 2) {
            System.out.println("The number of types is smaller than 2!");
            return null;
        }
        int tempN = this.n;
        while (this.n > 2) {
            String outToScreen = "    " + (tempN - n + 1) + " / " + (tempN - 1);
            System.out.print(outToScreen);
            // 1. gen the Q-matrix
            genQcriterion();
            // 2.find the minimum value and its idx
            // idxj > idxi
            int idxi = 0, idxj = 1;
            double minimum = this.qmatrix[idxi][idxj];
            for (int i = 0; i < this.n; i++) {
                for (int j = i + 1; j < this.n; j++) {
                    if (this.qmatrix[i][j] < minimum) {
                        idxi = i;
                        idxj = j;
                        minimum = qmatrix[idxi][idxj];
                    }
                }
            }
            // 3. combine the two nodes
            node[] newNodes = new node[this.n-1];
            node newnode = new midnode(this.nodes[idxi], this.nodes[idxj], this.global_n++);
            // compute the length of branch
            double vx = 0.5 * dmatrix[idxi][idxj] + (numsMatrix[idxi] - numsMatrix[idxj]) / (2 * n - 4);
            double vy = dmatrix[idxi][idxj] - vx;
            this.nodes[idxi].setLen(vx);
            this.nodes[idxj].setLen(vy);
            int[] treelist = {this.nodes[idxi].getNum(), this.nodes[idxj].getNum(), this.global_n - 1};
            this.TreeList.add(treelist.clone());

            // renew the nodes arrays
            for (int i = 0; i < this.n; i++) {
                if (i == idxi) { newNodes[i] = newnode; }
                else if (i < idxj) { newNodes[i] = this.nodes[i]; }
                else if (i > idxj) { newNodes[i-1] = this.nodes[i]; }
            }
            this.nodes = newNodes;

            // 4. renew upright area of the distance matrix
            double[][] newMatrixs = new double[this.n-1][this.n-1];
            for (int i = 0; i < this.n - 1; i++) {
                for (int j = i + 1; j < this.n -1; j++) {
                    if (i != idxi && j != idxi) {
                        int ii = i < idxj ? i : i + 1;
                        int jj = j < idxj ? j : j + 1;
                        newMatrixs[i][j] = newMatrixs[j][i] = this.dmatrix[ii][jj];
                    }
                }
            }
            for (int i = 0; i < idxi; i++) {
                newMatrixs[i][idxi] = newMatrixs[idxi][i] = 0.5 * (dmatrix[idxi][i] + dmatrix[idxj][i] - dmatrix[idxi][idxj]);
            }
            int counter = 0;
            for (int j = idxi + 1; j < this.n; j++) {
                if (j == idxj) {
                    counter++;
                    continue;
                }
                newMatrixs[idxi][j - counter] = newMatrixs[j - counter][idxi] = 0.5 * (dmatrix[idxi][j] + dmatrix[idxj][j] - dmatrix[idxi][idxj]);
            }
            this.dmatrix = newMatrixs;

            // renew n
            this.n--;

            System.out.print("\b".repeat(outToScreen.length()));
        }
        System.out.println("    " + (tempN - n + 1) + " / " + (tempN - 1) + "\n");
        midnode temp = new midnode(this.nodes[0], this.nodes[1], this.global_n);
        double len = this.dmatrix[0][1]/2;
        this.nodes[0].setLen(len);
        this.nodes[1].setLen(len);
        int[] treelist = {this.nodes[0].getNum(), this.nodes[1].getNum(), this.global_n};
        this.TreeList.add(treelist.clone());
        return temp;    
    }
}
