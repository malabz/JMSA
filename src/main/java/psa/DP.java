package psa;

import java.util.Arrays;


/** 
 * Pair sequence alignment (PSA) 
 * Affine gap penalty
*/
public class DP {
    
    private final int ms = 7, mis = -3;
    private final int d = 13, e = 2;
    private final String A, B;
    private String alignA = "", alignB = "";
    private boolean state = false;

    /**
     * input string A and string B
     * @param A
     * @param B
     */
    public DP(String A, String B) {
        if (A.length() > B.length() && B.length() > 0) {
            this.A = B;
            this.B = A;
            this.state = true;
        }
        else{
            this.A = A;
            this.B = B;
        }
        Align();
    }

    /**
     * To get the aligned results.
     */
    public String[] getStrAlign() {
        return new String[]{alignA, alignB};
    }

    private int Match(char a, char b) {
        return a == b ? this.ms : this.mis;
    }

    private void TraceBack(float[][][] p, int m, int n) {
        int channel = p[1][m][n] >= p[0][m][n] ? (p[1][m][n] >= p[2][m][n] ? 1 : 2) : (p[0][m][n] >= p[2][m][n] ? 0 : 2);
        StringBuilder alA = new StringBuilder();
        StringBuilder alB = new StringBuilder();
        int i = m, j = n;
        while (i > 0 || j > 0) {
            if (channel == 1 && j > 0) {
                channel = -1;
                if (p[1][i][j] == p[1][i][j-1] - e) channel = 1;
                else if (i > 0 && j > 1 && p[1][i][j] == p[0][i][j-1] - d) channel = 0;
                alA.insert(0, "-");
                alB.insert(0, B.charAt(j-1));
                j--;
            }
            else if (channel == 0 && i > 0 && j > 0) {
                channel = -1;
                if (i > 1 && j > 1 && p[0][i][j] == p[0][i-1][j-1] + Match(A.charAt(i-1), B.charAt(j-1))) {
                    channel = 0;
                }
                else if (j > 1 && p[0][i][j] == p[1][i-1][j-1] + Match(A.charAt(i-1), B.charAt(j-1))) {
                    channel = 1;
                }
                else if (i > 1 && p[0][i][j] == p[2][i-1][j-1] + Match(A.charAt(i-1), B.charAt(j-1))) {
                    channel = 2;
                }
                alA.insert(0, A.charAt(i-1));
                alB.insert(0, B.charAt(j-1));
                i--;
                j--;
            }
            else if (channel == 2 && i > 0) {
                channel = -1;
                if (p[2][i][j] == p[2][i-1][j] - e) channel = 2;
                else if (i > 1 && j > 0 && p[2][i][j] == p[0][i-1][j] - d) channel = 0;
                alA.insert(0, A.charAt(i-1));
                alB.insert(0, "-");
                i--;
            }
            else {
                System.out.println("(DP) Trace Back is wrong!");
                System.exit(0);
            }
        }
        alignA = alA.toString();
        alignB = alB.toString();
        if (state) {
            alignA = alB.toString();
            alignB = alA.toString();
        }
    }

    private float[][][] Init(int m, int n) {
        float[][][] p = new float[3][m+1][n+1];
        for (float[][] l : p) {
            for (int i = 0; i < m + 1; i++) {
                Arrays.fill(l[i], Float.NEGATIVE_INFINITY);
            }
        }
        p[0][0][0] = 0;
        for (int j = 1; j < n + 1; j++) {
            p[1][0][j] = - e * j;
        }
        for (int i = 1; i < m + 1; i++) {
            p[2][i][0] = - e * i;
        }
        return p;
    }

    private float Max3(float p0, float p1, float p2) {
        return Math.max(p0, Math.max(p1, p2));
    }

    private void Align() {
        // n >= m
        int m = A.length();
        int n = B.length();
        // special instances
        if (m == 0 && n == 0) { return; }
        else if (m == 0) {
            alignA = "-".repeat(n);
            alignB = B;
            return;
        }
        else if (n == 0) {
            alignB = "-".repeat(m);
            alignA = A;
            return;
        }
        float[][][] p = Init(m, n);
        for (int i = 1; i < m + 1; i++) {
            for (int j = 1; j < n + 1; j++) {
                // p[0] : A[i] ~ B[j]
                p[0][i][j] = Max3(p[0][i-1][j-1], p[1][i-1][j-1], p[2][i-1][j-1]) + Match(A.charAt(i-1), B.charAt(j-1));
                // p[1] : B[j] ~ -
                p[1][i][j] = Math.max(p[0][i][j-1] - d, p[1][i][j-1] - e);
                // p[2] : A[j] ~ -
                p[2][i][j] = Math.max(p[0][i-1][j] - d, p[2][i-1][j] - e);
            }
        }
        TraceBack(p, m, n);
    }
}
