package hierCluster;

import java.util.ArrayList;
import java.util.List;

public class upgma {
    private double[][] dmatrix;
    // private String[] names;
    public List<int[]> TreeList;
    private node[] nodes;
    private int[] nums;
    private int n, global_n;

    /**
     * 
     * @param matrix
     * @param names
     */
    public upgma(double[][] matrix, String[] names) {
        this.dmatrix = matrix;
        // this.names = names;
        this.n = names.length;
        this.nums = new int[this.n];
        for (int i = 0; i < n; i++) {
            nums[i] = 1;
        }
        this.nodes = new node[n];
        for (int i = 0; i < n; i++) {
            node temp = new leafnode(names[i], i);
            this.nodes[i] = temp;
        }
        this.global_n = this.n;
        this.TreeList = new ArrayList<>();
    }

    public upgma(double[][] matrix) {
        this.dmatrix = matrix;
        this.n = matrix.length;
        this.nums = new int[this.n];
        for (int i = 0; i < n; i++) {
            this.nums[i] = 1;
        }
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

    /**
     *
     */
    public void genTree() {
        if (this.n < 2) {
            System.out.println("The number of types is smaller than 2!");
            return;
        }
        while (this.n > 2) {    
            // find the minimum value and its idx
            // idxj > idxi
            int idxi = 0;
            int idxj = 1;
            double minimum = this.dmatrix[idxi][idxj];

            for (int i = 0; i < this.n; i++) {
                for (int j = i + 1; j < this.n; j++) {
                    if (this.dmatrix[i][j] < minimum) {
                        idxi = i;
                        idxj = j;
                        minimum = dmatrix[idxi][idxj];
                    }
                }
            }

            // combine the two nodes
            node[] newNodes = new node[this.n-1];
            node newnode = new midnode(this.nodes[idxi], this.nodes[idxj], this.global_n++);
            this.nodes[idxi].setLen(minimum/2 - this.nodes[idxi].getDistance());
            this.nodes[idxj].setLen(minimum/2 - this.nodes[idxj].getDistance());
            int[] treelist = {this.nodes[idxi].getNum(), this.nodes[idxj].getNum(), this.global_n - 1};
            // System.out.println(treelist[0]+", "+treelist[1]+", "+treelist[2]);

            this.TreeList.add(treelist.clone());

            // renew the nodes arrays
            for (int i = 0; i < this.n; i++) {
                if (i == idxi) { newNodes[i] = newnode; }
                else if (i < idxj) { newNodes[i] = this.nodes[i]; }
                else if (i > idxj) { newNodes[i-1] = this.nodes[i]; }
            }
            this.nodes = newNodes;

            // renew upright area of the distance matrix
            double[][] newMatrixs = new double[this.n-1][this.n-1];
            for (int i = 0; i < this.n - 1; i++) {
                for (int j = i + 1; j < this.n -1; j++) {
                    if (i != idxi && j != idxi) {
                        int ii = i < idxj ? i : i + 1;
                        int jj = j < idxj ? j : j + 1;
                        newMatrixs[i][j] = this.dmatrix[ii][jj];
                    }
                }
            }
            for (int i = 0; i < idxi; i++) {
                newMatrixs[i][idxi] = (dmatrix[i][idxi]*nums[idxi] + dmatrix[i][idxj]*nums[idxj])/(nums[idxi] + nums[idxj]);
            }
            int counter = 0;
            for (int j = idxi + 1; j < this.n; j ++) {
                if (j < idxj) {
                    dmatrix[idxj][j] = dmatrix[j][idxj];
                }
                else if (j == idxj) {
                    counter++;
                    continue;
                }
                newMatrixs[idxi][j-counter] = (dmatrix[idxi][j]*nums[idxi] + dmatrix[idxj][j]*nums[idxj])/(nums[idxi] + nums[idxj]);
            }
            this.dmatrix = newMatrixs;

            // renew n
            this.n--;

            // renew nums
            this.nums[idxi]++;
            int[] newnums = new int[this.n];
            for (int i = 0; i < this.n; i++) {
                if (i < idxj) {
                    newnums[i] = this.nums[i];
                }
                else {
                    newnums[i] = this.nums[i+1];
                }
            }
            this.nums = newnums;
        }
        new midnode(this.nodes[0], this.nodes[1], this.global_n);
        double len = this.dmatrix[0][1]/2;
        this.nodes[0].setLen(len - this.nodes[0].getDistance());
        this.nodes[1].setLen(len - this.nodes[1].getDistance());
        int[] treelist = {this.nodes[0].getNum(), this.nodes[1].getNum(), this.global_n};
        this.TreeList.add(treelist.clone());
    }
}