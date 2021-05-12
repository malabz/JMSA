package psa;

public abstract class dp {
    protected String A, B;
    protected String[] alignAB;
    protected float[][] pm;
    protected int ms = 7, mis = -3;
    protected int d = 5;
    protected boolean state;
    
    public dp(String A, String B) {
        if (A.length() < B.length()) {
            this.B = A;
            this.A = B;
            this.state = true;
        }
        else {
            this.A = A;
            this.B = B;
            this.state = false;
        }
    }

    protected void Init() {
        pm[0][0] = 0;
        for (int i = 1; i < pm.length; i++) {
            pm[i][0] = pm[i-1][0] - d;
        }
        for (int j = 1; j < pm[0].length; j++) {
            pm[0][j] = pm[0][j-1] - d;
        }
    }

    protected int Match(int i, int j) {
        return A.charAt(i) == B.charAt(j) ? ms : mis;
    }

    protected float Maxfloat3(float p0, float p1, float p2) {
        return Math.max(Math.max(p0, p1), Math.max(p2, 0));
    }

    public String[] getStrAlign() { return alignAB; }
    protected abstract void TraceBack();
    protected abstract void Align();
}
