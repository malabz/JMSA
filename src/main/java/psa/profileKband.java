package psa;

import io.str;

import java.util.Arrays;


public class profileKband extends kb {
    private final String[] A, B;
    private String[] alignA, alignB;
    private final int rnumA, rnumB;
    private int numA = 0, numB = 0;
    private int[][] alphA, alphB;
    private boolean state = false;


    /**
     *
     * @param A
     * @param B
     * @param alphA
     * @param alphB
     * @param kk
     * @param rnumA
     * @param rnumB
     */
    public profileKband(String[] A, String[] B, int[][] alphA, int[][] alphB, int kk, int rnumA, int rnumB) {
        if (A[0].length() > B[0].length() && B[0].length() > 0) {
            this.A = B;
            this.B = A;
            this.alphA = alphB;
            this.alphB = alphA;
            this.rnumA = rnumB;
            this.rnumB = rnumA;
            this.state = true;
        }
        else {
            this.A = A;
            this.B = B;
            this.alphA = alphA;
            this.alphB = alphB;
            this.rnumA = rnumA;
            this.rnumB = rnumB;
        }
        this.InitalAB();
        this.Align(kk);
    }

    /**
     *
     * @param A
     * @param B
     * @param alphA
     * @param alphB
     * @param rnumA
     * @param rnumB
     */
    public profileKband(String[] A, String[] B, int[][] alphA, int[][] alphB, int rnumA, int rnumB) {
        if (A[0].length() > B[0].length() && B[0].length() > 0) {
            this.A = B;
            this.B = A;
            this.alphA = alphB;
            this.alphB = alphA;
            this.rnumA = rnumB;
            this.rnumB = rnumA;
            this.state = true;
        }
        else {
            this.A = A;
            this.B = B;
            this.alphA = alphA;
            this.alphB = alphB;
            this.rnumA = rnumA;
            this.rnumB = rnumB;
        }
        this.InitalAB();
        this.Align(1);
    }

    public int[][] getAlpha() {
        int[][] alphC = new int[alphA.length][];
        if (numA + numB > 128) {
            int len = Math.max(Integer.toBinaryString(numA+numB).length() - 7, 0);
            for (int i = 0; i < alphC.length; i++) {
                for (int j = 0; j < alphC[i].length; j++) { 
                    alphC[i][j] = ((alphA[i][j] + alphB[i][j]) >> len);
                }
            }
        }
        else {
            for (int i = 0; i < alphC.length; i++) {
                for (int j = 0; j < alphC[i].length; j++) { 
                    alphC[i][j] = alphA[i][j] + alphB[i][j];
                }
            }
        }
        
        return alphC;
    }

    /**
     * To get the aligned results.
     */
    public String[][] getStrAlign() {
        if (state) return new String[][]{alignB, alignA};
        return new String[][]{alignA, alignB};
    }

    public String[] getStrsAlign() {
        String[] res = new String[alignA.length + alignB.length];
        if (state) { 
            System.arraycopy(alignB, 0, res, 0, alignB.length);
            System.arraycopy(alignA, 0, res, alignB.length, alignA.length);
        }
        else {
            System.arraycopy(alignA, 0, res, 0, alignA.length);
            System.arraycopy(alignB, 0, res, alignA.length, alignB.length);
        }
        return res;
    }

    /**
     * init the result arrays
     */
    private void InitalAB() {
        // numA 目前存在的数量
        // rnumA 实际上的数量
        numA = this.A[0].length();
        numB = this.B[0].length();
        this.alignA = new String[numA];
        this.alignB = new String[numB];
        for (int i = 0; i < numA; i++) this.alignA[i] = "";
        for (int i = 0; i < numB; i++) this.alignB[i] = "";
        for (int i = 0; i < alphA[0].length; i++) this.numA += alphA[0][i];
        for (int i = 0; i < alphB[0].length; i++) this.numB += alphB[0][i];
    }

    private int Match(int idxm, int idxn) {
        int[] tempA = alphA[idxm], tempB = alphB[idxn];
        int results = 0, len = tempA.length;
        if (rnumA > 50 || rnumB > 50) {
            for (int i = 0; i < len - 2; i++) results += (tempA[i] * tempB[i]);
            results *= (this.ms - this.mis);
            results += (numA - tempA[len-1] - tempA[len-2]) * (numB - tempB[len-1] - tempB[len-2]) * this.mis;
            results += (numA - tempA[len-2]) * (tempB[len-1] * this.ms - tempB[len-2] * this.e);
            results -= tempA[len-2] * (numB - tempB[len-2]) * this.e;
            results += tempA[len-1] * (numB - tempB[len-1] - tempB[len-2]) * this.ms;
            return results/(numA*numB);
        }
        else {
            int sqA = 0, sqB = 0;
            for (int i = 0; i < len; i++) {
                sqA += Math.pow(tempA[i], 2);
                sqB += Math.pow(tempB[i], 2);
                results += (tempA[i] * tempB[i]); 
            }
            double res = results / Math.sqrt(sqA * sqB);
            return (int) (res * ms);
        }
    }

    @Override
    protected void TraceBack(float[][][] pm, int k) {
        int m = this.A[0].length(), n = this.B[0].length();
        int diff = n - m;

        int i = m, bj = n, j = diff + k;
        int channel = ChooseMax(pm[0][i][j], pm[1][i][j], pm[2][i][j]);

        StringBuilder[] alignA = new StringBuilder[numA];
        StringBuilder[] alignB = new StringBuilder[numB];
        for (int sba = 0; sba < numA; sba++) { alignA[sba] = new StringBuilder(); }
        for (int sbb = 0; sbb < numB; sbb++) { alignB[sbb] = new StringBuilder(); }

        while (i > 0 || j > k) {
            if (channel == 1 && j > 0) {
                channel = -1;
                if (pm[1][i][j] == pm[1][i][j-1] - e) { channel = 1; }
                else if (i >= 1 && pm[1][i][j] == pm[0][i][j-1] - d ) { channel = 0; }
                for (int idxA = 0; idxA < numA; idxA++) {
                    alignA[idxA].insert(0, "-");
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    alignB[idxB].insert(0, this.B[idxB].charAt(bj-1));
                }
                bj--;
                j--;
            }
            else if (channel == 0 && i > 0 && j >= 0) {
                channel = -1;
                int match = Match(i-1, bj-1);
                if (i > 1 && pm[0][i][j] == pm[0][i-1][j] + match) {
                    channel = 0;
                }
                else if (j > 0 && pm[0][i][j] == pm[1][i-1][j] + match) {
                    channel = 1;
                }
                else if (i > 1 && pm[0][i][j] == pm[2][i-1][j] + match) {
                    channel = 2;
                }
                for (int idxA = 0; idxA < numA; idxA++) {
                    alignA[idxA].insert(0, this.A[idxA].charAt(i-1));
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    alignB[idxB].insert(0, this.B[idxB].charAt(bj-1));
                }
                i--;
                bj--;
            }

            else if (channel == 2 && i > 0 && (j + 1) <= (2 * k + diff)) {
                channel = -1;
                if (pm[2][i][j] == pm[2][i-1][j+1] - e) { channel = 2; }
                else if (i > 1 && pm[2][i][j] == pm[0][i - 1][j + 1] - d) { channel = 0; }
                for (int idxA = 0; idxA < numA; idxA++) {
                    alignA[idxA].insert(0, this.A[idxA].charAt(i-1));
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    alignB[idxB].insert(0, "-");
                }
                i--;
                j++;
            }
            else {
                throw new IllegalStateException("channel = " + channel);
            }
        }
        for (i = 0; i < numA; i++) { this.alignA[i] = alignA[i].toString(); }
        for (j = 0; j < numB; j++) { this.alignB[j] = alignB[j].toString(); }
    }

    private void Align(int kk) {
        int m = this.A[0].length(), n = this.B[0].length();
        int diff = n - m, k = kk;

        if (m == 0 && n == 0) { return; }
        else if (m == 0) {
            Arrays.fill(this.alignA, str.repeat("-", n));
            this.alignB = this.B;
            return;
        }
        else if (n == 0) {
            Arrays.fill(this.alignB, str.repeat("-", m));
            this.alignA = this.A;
            return;
        }

        float valueOld = Float.NEGATIVE_INFINITY, valueNew;
        float[][][] pm = new float[3][m+1][diff+2*k+1];

        int maxk = Math.min(m, Math.max(m/5, kk));
        while (k <= maxk) {
            this.Init(pm, k, diff);

            for(int i = 1; i<m+1; ++i) {
                for(int ii = -k; ii<diff+k+1; ++ii) {
                    int j = ii;
                    if (1<=j+i && j+i<=n) {
                        j += k;
                        pm[0][i][j] = Maxfloat3(pm[0][i-1][j], pm[1][i-1][j], pm[2][i-1][j]) + Match(i-1, j+i-k-1);
                        
                        if (InsiderStrip(i, j+i-k-1, k, diff)) {
                            // p[1] : B[j] ~ -
                            pm[1][i][j] = Math.max(pm[0][i][j-1] - d, pm[1][i][j-1] - e);
                        }

                        if (InsiderStrip(i-1, j+i-k, k, diff)) {
                            // p[2] : A[j] ~ -
                            pm[2][i][j] = Math.max(pm[0][i-1][j+1] - d, pm[2][i-1][j+1] - e);
                        }
                    }
                }
            }
            valueNew = Maxfloat3(pm[0][m][diff+k], pm[1][m][diff+k], pm[2][m][diff+k]);
            if ((int) valueNew == (int) valueOld) break;
            else {
                valueOld = valueNew;
                k *= 2;
                if (k <= maxk) pm = new float[3][m+1][diff+2*k+1];
                else { k /= 2; break; }
            }
        }
        TraceBack(pm, k);
    }
}