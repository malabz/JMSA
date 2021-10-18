package psa;

import java.util.Arrays;

/**
 * This is the abstract class of the kband dynamic programming.
 */
public abstract class kb {

    protected final int ms = 7, mis = -3;
    protected final int d = 13, e = 2;

    /**
     * judge the idx whether in the matrix
     * @return boolean
     */
    protected boolean InsiderStrip(int i, int j, int k, int diff){
        return (-k <= (j - i)) && ((j - i )<= (k + diff));
    }

    /**
     * init the dp matrix
     */
    protected void Init(float[][][] pm, int k, int diff) {
        for (float[][] floats : pm) {
            for (int i2 = 0; i2 < pm[0].length; ++i2) {
                Arrays.fill(floats[i2], Float.NEGATIVE_INFINITY);
            }
        }
        pm[0][0][k] = 0;
        for(int j = 1; j < k+1+diff; ++j) {
            pm[1][0][j+k] = -this.d - this.e * (j-1);
        }
        for(int i = 1; i < k+1; ++i) {
            pm[2][i][k-i] = -this.d - this.e * (i-1);
        }
    }

    /**
     * choose one to trace back
     * @return 0/1/2
     */
    protected int ChooseMax(float p0, float p1, float p2) {
        return p0 >= p1 ? (p0 >= p2 ? 0 : 2) : (p1 >= p2 ? 1 : 2);
    }

    /**
     * compare the three float and return the max one
     */
    protected float Maxfloat3(float p0, float p1, float p2) {
        return Math.max(Math.max(p0, p1), p2);
    }

    protected abstract void TraceBack(float[][][] pm, int k);
}
