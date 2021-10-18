package psa;

import io.str;

import java.util.Arrays;
import java.util.HashMap;


public class multiKband extends kb {
    private final String[] A, B;
    private String[] alignA, alignB;
    private final int numA, numB;
    // private int lenA, lenB, rnumB, rnumA;
    private int[][] alphA, alphB;
    private boolean state = false;
    HashMap<Character, Integer> mapAlph;
    

    /**
     * 
     * @param A
     * @param B
     * @param alphabet
     */
    public multiKband(String[] A, String[] B, char[] alphabet, int kk) {
        if (A[0].length() > B[0].length() && B[0].length() > 0) {
            this.A = B;
            this.B = A;
            this.numA = B.length;
            this.numB = A.length;
            this.state = true;
        }
        else {
            this.A = A;
            this.B = B;
            this.numA = A.length;
            this.numB = B.length;
        }
        this.InitalAB();
        this.CountAlphabet(alphabet);
        this.Align(kk);
    }

    /**
     * 
     * @param A
     * @param B
     * @param alphabet
     */
    public multiKband(String[] A, String[] B, char[] alphabet) {
        if (A[0].length() > B[0].length() && B[0].length() > 0) {
            this.A = B;
            this.B = A;
            this.numA = B.length;
            this.numB = A.length;
            this.state = true;
        }
        else {
            this.A = A;
            this.B = B;
            this.numA = A.length;
            this.numB = B.length;
        }
        this.InitalAB();
        this.CountAlphabet(alphabet);
        this.Align(1);
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

    private void CountAlphabet(char[] alphabet) {
        mapAlph = new HashMap<>();
        for (int i = 0; i < alphabet.length; i++) mapAlph.put(alphabet[i], i);
        mapAlph.put('-', alphabet.length);
        mapAlph.put('n', alphabet.length + 1);
        alphA = new int[A[0].length()][alphabet.length + 2];
        alphB = new int[B[0].length()][alphabet.length + 2];
        // lenA = Math.max(Integer.toBinaryString(numA).length() - 7, 0);
        // lenB = Math.max(Integer.toBinaryString(numB).length() - 7, 0);
        // rnumA = numA >> lenA;
        // rnumB = numB >> lenB;
        for (int i = 0; i < A[0].length(); i++) {
            for (String s : A) {
                char c = mapAlph.containsKey(s.charAt(i)) ? s.charAt(i) : 'n';
                alphA[i][mapAlph.get(c)] += 1;
            }
            // for (int j = 0; j < mapAlph.size(); j++) { alphA[i][j] >>= lenA; }
        }
        for (int i = 0; i < B[0].length(); i++) {
            for (String s : B) {
                char c = mapAlph.containsKey(s.charAt(i)) ? s.charAt(i) : 'n';
                alphB[i][mapAlph.get(c)] += 1;
            }
            // for (int j = 0; j < mapAlph.size(); j++) { alphB[i][j] >>= lenB; }
        }
    }

    /**
     * init the result arrays
     */
    private void InitalAB() {
        this.alignA = new String[numA];
        this.alignB = new String[numB];
        for (int i = 0; i < numA; i++) this.alignA[i] = "";
        for (int i = 0; i < numB; i++) this.alignB[i] = "";
    }

    private int Match(int idxm, int idxn) {
        int[] tempA = alphA[idxm], tempB = alphB[idxn];
        int results = 0, len = tempA.length;
        if (numA > 50 || numB > 50) {
            // int numA = rnumA, numB = rnumB;
            for (int i = 0; i < len - 2; i++) results += (tempA[i] * tempB[i]);
            results *= (this.ms - this.mis);
            results += (numA - tempA[len-1] - tempA[len-2]) * (numB - tempB[len-1] - tempB[len-2]) * this.mis;
            results += (numA - tempA[len-2]) * (tempB[len-1] * this.ms - tempB[len-2] * this.e);
            results -= tempA[len-2] * (numB - tempB[len-2]) * this.e;
            results += tempA[len-1] * (numB - tempB[len-1] - tempB[len-2]) * this.ms;
            // numA = Math.max(numA, 1);
            // numB = Math.max(numB, 1);
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