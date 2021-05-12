package psa;

public class smithWaterman extends dp {

    public smithWaterman(String A, String B) {
        super(A, B);
        Align();
    }
    
    private int[] findMaxidx(int m, int n) {
        int[] idx = new int[]{1,1};
        for (int i = 1; i < m; i++) {
            for (int j = 1; j < n; j++) {
                if (pm[i][j] > pm[idx[0]][idx[1]]) {
                    idx[0] = i;
                    idx[1] = j;
                }
            }
        }
        return idx;
    }

    @Override
    protected void TraceBack() {
        int[] idx = findMaxidx(pm.length, pm[0].length);
        int pi = idx[0];
        int pj = idx[1];
        StringBuilder alA = new StringBuilder();
        StringBuilder alB = new StringBuilder();

        while (pi > 0 && pj > 0) {
            if (pm[pi][pj] == 0) { break; }
            if (pm[pi][pj] == pm[pi-1][pj-1] + Match(pj-1, pi-1)) {
                alA.append(A.charAt(--pj));
                alB.append(B.charAt(--pi));
            }
            else if (pm[pi][pj] == pm[pi][pj-1] - d) { 
                alA.append('-');
                alB.append(B.charAt(--pi));
            }
            else if (pm[pi][pj] == pm[pi-1][pj] - d) { 
                alA.append(A.charAt(--pj));
                alB.append('-');
            }
        }
        assert(pm[pi][pj] == 0);
    }

    @Override
    protected void Align() {
        // n >= m
        int n = A.length();
        int m = B.length();

        pm = new float[m+1][n+1];
        // Init();

        for (int j = 1; j <= m; j++) {
            for (int i = 1; i <= n; i++) {
                pm[j][i] = Maxfloat3(pm[j-1][i-1] + Match(i-1, j-1), pm[j-1][i] - d, pm[j][i-1] - d);
            }
        }
        TraceBack();
    }
}
