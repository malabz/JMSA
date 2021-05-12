package psa;

/** 
 * Pair sequence alignment (PSA) 
 * Affine gap penalty + Kband
 * 
*/
public class Kband extends kb {
    
    private String A, B;
    private String alignA = "", alignB = "";
    private boolean state = false;

    private int score;
    private float[][][] pm = new float[3][][];

    /**
     * input string A and string B
     * @param A
     * @param B
     */
    public Kband(String A, String B) {
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
     * To get the align results.
     * @return
     */
    public String[] getStrAlign() {
        String[] res = new String[2];
        res[0] = alignA;
        res[1] = alignB;
        return res;
    }

    public int getScore() {
        return score;
    }

    private int Match(char a, char b) {
        return a==b ? this.ms : this.mis;
    }

    @Override
    protected void TraceBack(int k) {
        int n = this.B.length(), m = this.A.length();
        int diff = n - m;
        int i = m, bj = n, j = diff + k;
        int channel = ChooseMax(pm[0][i][j], pm[1][i][j], pm[2][i][j]);
        while (i > 0 || j > k) {
            if (channel == 0 && i > 0 && j >= 0) {
                if (pm[0][i][j] == pm[0][i-1][j] + Match(A.charAt(i-1), B.charAt(bj-1)) && (i > 1) && (j >= 0)) {
                    channel = 0;
                }
                else if (pm[0][i][j] == pm[2][i-1][j] + Match(A.charAt(i-1), B.charAt(bj-1)) && i > 1) {
                    channel = 2;
                }
                else if (pm[0][i][j] == pm[1][i-1][j] + Match(A.charAt(i-1), B.charAt(bj-1)) && j > 0) {
                    channel = 1;
                }
                this.alignA += A.charAt(i-1);
                this.alignB += B.charAt(bj-1);
                i -= 1;
                bj -= 1;
            }
            else if (channel == 1 && j > 0) {
                if (pm[1][i][j] == pm[1][i][j-1] -e) {
                    channel = 1;
                }
                else if (pm[1][i][j] == pm[0][i][j-1] - d && i >= 1) {
                    channel = 0;
                }
                this.alignA += "-";
                this.alignB += B.charAt(bj-1);
                bj -= 1;
                j -= 1;
            }
            else if (channel == 2 && i > 0 && (j + 1) <= (2 * k + diff)) {
                if (pm[2][i][j] == pm[2][i-1][j+1] - e) {
                    channel = 2;
                }
                else if (pm[2][i][j] == pm[0][i-1][j+1] - d && i > 1 && j >= 0) {
                    channel = 0;
                }
                this.alignA += A.charAt(i-1);
                this.alignB += "-";
                i -= 1;
                j += 1;
            }
            else {
                System.out.println("Trace Back is wrong!");
                System.exit(0);
            }
        }

        this.alignA = new StringBuilder(this.alignA).reverse().toString();
        this.alignB = new StringBuilder(this.alignB).reverse().toString();

        if (this.state) {
            String temp = this.alignA;
            this.alignA = this.alignB;
            this.alignB = temp;
        }
    }

    @Override
    protected int Align() {
        int m = this.A.length(), n = this.B.length();
        int diff = n - m, k = 1;

        // len(A)=0 or len(B)=0
        if (m == 0 && 0== n) { 
            this.score = 0; 
            return 0; 
        }
        else if (m == 0) {
            // this.alignA = this.StrsMultiply(n);
            this.alignA = "-".repeat(n);
            this.alignB = this.B;
            return 0;
        }
        else if (n == 0) {
            // this.alignB = this.StrsMultiply(m);
            this.alignB = "-".repeat(m);
            this.alignA = this.A;
            return 0;
        }

        float valueOld = Float.POSITIVE_INFINITY, valueNew = 0;

        float[][][] pm = new float[3][m+1][diff+2*k+1];

        int maxk = Math.min(m, Math.max(m/5, 10));
        while (k <= maxk) {
            // init
            this.Init(pm, k, diff);

            for(int i = 1; i<m+1; ++i) {
                for(int ii = -k; ii<diff+k+1; ++ii) {
                    int j = ii;
                    if (1<=j+i && j+i<=n) {
                        j += k;
                        // p[0] : A[i] ~ B[j]
                        pm[0][i][j] = Maxfloat3(pm[0][i-1][j], pm[1][i-1][j], pm[2][i-1][j]) + Match(A.charAt(i-1), B.charAt(j+i-k-1));
                        
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
            if ((int)valueNew == (int)valueOld) {break;}
            else {
                valueOld = valueNew;
                k *= 2;
                if (k <= maxk) {pm = new float[3][m+1][diff+2*k+1];}
                else {
                    k /= 2;
                    break;
                }
            }
        }
        // System.out.println(m/k);
        this.pm = pm;
        this.score = (int)valueNew;
        TraceBack(k);
        return 1;
    }
}