package psa;

import java.util.Arrays;

/**
 * This is the abstract class of the kband dynamic programming.
 */
public abstract class kb {

    protected int ms = 7, mis = -3;
    protected int d = 11, e = 2;

    /**
     * judge the idx whether in the matrix
     * @param i
     * @param j
     * @param k
     * @param diff
     * @return boolean
     */
    protected boolean InsiderStrip(int i, int j, int k, int diff){
        return (-k <= (j - i)) && ((j - i )<= (k + diff));
    }

    /**
     * init the dp matrix
     * @param pm
     * @param k
     * @param diff
     */
    protected void Init(float[][][] pm, int k, int diff) {
        for (int i1 = 0; i1 < pm.length; ++i1) {
            for (int i2 = 0; i2 < pm[0].length; ++i2) {
                Arrays.fill(pm[i1][i2], Float.NEGATIVE_INFINITY);
            }
        }
        pm[0][0][k] = 0;
        for(int j = 1; j < k+1+diff; ++j) {
            pm[1][0][j+k] = -this.d - this.e*(j-1);
        }
        for(int i = 1; i < k+1; ++i) {
            pm[2][i][k-i] = -this.d - this.e * (i-1);
        }
    }

    /**
     * choose one to trace back
     * @param p0
     * @param p1
     * @param p2
     * @return 0/1/2
     */
    protected int ChooseMax(float p0, float p1, float p2) {
        return p0 >= p1 ? (p0 >= p2 ? 0 : 2) : (p1 >= p2 ? 1 : 2);
    }

    /**
     * compare the three float and return the max one
     * @param p0
     * @param p1
     * @param p2
     * @return
     */
    protected float Maxfloat3(float p0, float p1, float p2) {
        return Math.max(Math.max(p0, p1), p2);
    }

    protected abstract void TraceBack(int k);
    protected abstract int Align();
}
