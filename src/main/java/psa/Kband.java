package psa;

/** 
 * Pair sequence alignment (PSA) 
 * Affine gap penalty + Kband
 * 
*/
public class Kband extends kb {
    
    private final String A, B;
    private String alignA = "", alignB = "";
    private boolean state = false;


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
     * To get the aligned results.
     */
    public String[] getStrAlign() {
        return new String[]{alignA, alignB};
    }

    private int Match(char a, char b) {
        return (a == b || a == 'N' || b == 'N') ? this.ms : this.mis;
    }

    @Override
    protected void TraceBack(float[][][] pm, int k) {
        int n = this.B.length(), m = this.A.length();
        int diff = n - m;
        int i = m, bj = n, j = diff + k;
        int channel = ChooseMax(pm[0][i][j], pm[1][i][j], pm[2][i][j]);
        StringBuilder alignA = new StringBuilder();
        StringBuilder alignB = new StringBuilder();

        while (i > 0 || j > k) {
            if (channel == 1 && j > 0) {
                channel = -1;
                if (pm[1][i][j] == pm[1][i][j-1] - e) { channel = 1; }
                else if (i >= 1 && pm[1][i][j] == pm[0][i][j-1] - d) { channel = 0; }
                alignA.insert(0, "-");
                alignB.insert(0, B.charAt(bj-1));
                bj--;
                j--;
            }
            else if (channel == 0 && i > 0 && j >= 0) {
                channel = -1;
                if (i > 1 && pm[0][i][j] == pm[0][i-1][j] + Match(A.charAt(i-1), B.charAt(bj-1))) {
                    channel = 0;
                }
                else if (j > 0 && pm[0][i][j] == pm[1][i-1][j] + Match(A.charAt(i-1), B.charAt(bj-1))) {
                    channel = 1;
                }
                else if (i > 1 && pm[0][i][j] == pm[2][i-1][j] + Match(A.charAt(i-1), B.charAt(bj-1))) {
                    channel = 2;
                }
                alignA.insert(0, A.charAt(i-1));
                alignB.insert(0, B.charAt(bj-1));
                i--;
                bj--;
            }
            else if (channel == 2 && i > 0 && (j + 1) <= (2 * k + diff)) {
                channel = -1;
                if (pm[2][i][j] == pm[2][i-1][j+1] - e) { channel = 2; }
                else if (i > 1 && pm[2][i][j] == pm[0][i-1][j+1] - d) { channel = 0; }
                alignA.insert(0, A.charAt(i-1));
                alignB.insert(0, "-");
                i--;
                j++;
            }
            else {
                System.out.println("(Kband) Trace Back is wrong!");
                System.exit(0);
            }
        }

        this.alignA = alignA.toString();
        this.alignB = alignB.toString();

        if (this.state) {
            this.alignA = alignB.toString();
            this.alignB = alignA.toString();
        }
    }

    private void Align() {
        int m = this.A.length(), n = this.B.length();
        int diff = n - m, k = 1;

        // len(A)=0 or len(B)=0
        if (m == 0 && 0 == n) { return; }
        else if (m == 0) {
            this.alignA = "-".repeat(n);
            this.alignB = this.B;
            return;
        }
        else if (n == 0) {
            this.alignB = "-".repeat(m);
            this.alignA = this.A;
            return;
        }

        float valueOld = Float.POSITIVE_INFINITY, valueNew;
        float[][][] pm = new float[3][m+1][diff+2*k+1];

        int maxk = Math.min(m, Math.max(m/5, 10));
        while (k <= maxk) {
            // init
            this.Init(pm, k, diff);

            for(int i = 1; i < m + 1; ++i) {
                for(int ii = -k; ii < diff + k + 1; ++ii) {
                    int j = ii;
                    if (1 <= j + i && j + i <= n) {
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
        TraceBack(pm, k);
    }
}