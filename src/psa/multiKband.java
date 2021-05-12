package psa;

import java.util.HashMap;


public class multiKband extends kb {
    private String[] A, B;
    private char[] alphabet;
    private int numA, numB;
    public String[] alignA, alignB;
    private boolean state = false;
    private float[][][] pm = new float[3][][];
    
    /**
     * 
     * @param A
     * @param B
     */
    public multiKband(String[] A, String[] B, char[] alphabet) {
        if (A[0].length() > B[0].length() && B[0].length() > 0) {
            this.A = B;
            this.B = A;
            this.numA = B.length;
            this.numB = A.length;
            InitalAB();
            this.alphabet = alphabet;
            this.state = true;
            Align();
        }
        else {
            this.A = A;
            this.B = B;
            this.numA = A.length;
            this.numB = B.length;
            InitalAB();
            this.alphabet = alphabet;
            Align();
        }
    }

    /**
     * To get the align results.
     * @return
     */
    public String[][] getStrAlign() {
        String[][] res = new String[2][alignA.length];
        res[0] = alignA;
        res[1] = alignB;
        return res;
    }

    /**
     * init the result arrays
     */
    private void InitalAB() {
        // System.out.print(numA + ":" + numB);
        this.alignA = new String[numA];
        this.alignB = new String[numB];
        for (int i = 0; i < numA; i++) {
            this.alignA[i] = "";
        }
        for (int i = 0; i < numB; i++) {
            this.alignB[i] = "";
        }
    }

    private int Match(int idxm, int idxn) {
        int results = 0;
        
        HashMap<Character, Integer> mapAlph = new HashMap<>();
        for (char c : alphabet) { mapAlph.put(c, 0); }
        HashMap<Character, Integer> score = new HashMap<>();
        for (char c : alphabet) { score.put(c, Integer.MIN_VALUE); }
        score.put('-', Integer.MIN_VALUE);
        score.put('n', Integer.MIN_VALUE);
        
        if (numA < numB) {
            for (int j = 0; j < numB; j++) {
                char bj = this.B[j].charAt(idxn);
                bj = (mapAlph.containsKey(bj) || bj == '-') ? bj : 'n'; 
                if (score.get(bj) > Integer.MIN_VALUE) {
                    results += score.get(bj);
                    continue;
                }
                int temp = 0;
                for (int i = 0; i < numA; i++) {
                    char ai = this.A[i].charAt(idxm);
                    if (mapAlph.containsKey(ai) && mapAlph.containsKey(bj)) {
                        temp += ai == bj ? this.ms : this.mis;
                    }
                    else if (ai == '-' || bj == '-') {
                        temp -= ai == bj ? 0 : this.e;
                    }
                    else {
                        temp += this.ms;
                    }
                }
                score.put(bj, temp);
                results += temp;
            }
            return results/(numA*numB);
        }
        for (int i = 0; i < numA; i++) {
            char ai = this.A[i].charAt(idxm);
            ai = (mapAlph.containsKey(ai) || ai == '-') ? ai : 'n'; 
            if (score.get(ai) > Integer.MIN_VALUE) {
                results += score.get(ai);
                continue;
            }
            int temp = 0;
            for (int j = 0; j < numB; j++) {
                char bj = this.B[j].charAt(idxn);
                if (mapAlph.containsKey(ai) && mapAlph.containsKey(bj)) {
                    temp += ai == bj ? this.ms : this.mis;
                }
                else if (ai == '-' || bj == '-') {
                    temp -= ai == bj ? 0 : this.e;
                }
                else {
                    temp += this.ms;
                }
            }
            score.put(ai, temp);
            results += temp;
        }
        return results/(numA*numB);
        /*
        // if (numA * numB < 0) {
        //     int results = 0;
        //     for (int i = 0; i < numA; i++) {
        //         char ai = this.A[i].charAt(idxm);
        //         for (int j = 0; j < numB; j++) {
        //             char bj = this.B[j].charAt(idxn);
        //             if (alNumA.containsKey(ai) && alNumA.containsKey(bj)) {
        //                 results += ai == bj ? this.ms : this.mis;
        //             }
        //             else if (ai == '-' || bj == '-') {
        //                 results -= ai == bj ? 0 : this.e;
        //             }
        //             else {
        //                 results += this.ms;
        //             }
        //         }
        //         results += results;
        //     }
        //     return results/(numA*numB);
        // }
        // else if (numA > 0) {
        //     HashMap<Character, Integer> score = new HashMap<>();
        //     for (char c : alphabet) { score.put(c, Integer.MIN_VALUE); }
        //     score.put('-', Integer.MIN_VALUE);
        //     score.put('n', Integer.MIN_VALUE);
        //     int results = 0;
            
        //     if (numA < numB) {
        //         for (int j = 0; j < numB; j++) {
        //             char bj = this.B[j].charAt(idxn);
        //             bj = (alNumA.containsKey(bj) || bj == '-') ? bj : 'n'; 
        //             if (score.get(bj) > Integer.MIN_VALUE) {
        //                 results += score.get(bj);
        //                 continue;
        //             }
        //             int temp = 0;
        //             for (int i = 0; i < numA; i++) {
        //                 char ai = this.A[i].charAt(idxm);
        //                 if (alNumA.containsKey(ai) && alNumA.containsKey(bj)) {
        //                     temp += ai == bj ? this.ms : this.mis;
        //                 }
        //                 else if (ai == '-' || bj == '-') {
        //                     temp -= ai == bj ? 0 : this.e;
        //                 }
        //                 else {
        //                     temp += this.ms;
        //                 }
        //             }
        //             score.put(bj, temp);
        //             results += temp;
        //         }
        //         return results/(numA*numB);
        //     }
        //     for (int i = 0; i < numA; i++) {
                
        //         char ai = this.A[i].charAt(idxm);
        //         ai = (alNumA.containsKey(ai) || ai == '-') ? ai : 'n'; 
        //         if (score.get(ai) > Integer.MIN_VALUE) {
        //             results += score.get(ai);
        //             continue;
        //         }
        //         int temp = 0;

        //         for (int j = 0; j < numB; j++) {
        //             char bj = this.B[j].charAt(idxn);
        //             if (alNumA.containsKey(ai) && alNumA.containsKey(bj)) {
        //                 temp += ai == bj ? this.ms : this.mis;
        //             }
        //             else if (ai == '-' || bj == '-') {
        //                 temp -= ai == bj ? 0 : this.e;
        //             }
        //             else {
        //                 temp += this.ms;
        //             }
        //         }
        //         score.put(ai, temp);
        //         results += temp;
        //     }
        //     return results/(numA*numB);
        // }
        
        // if (alNumA.containsKey('n') == false) { alNumA.put('n', 0); }
        // if (alNumA.containsKey('-') == false) { alNumA.put('-', 0); }
        // for (int i = 0; i < numA; i++) {
        //     char c = this.A[i].charAt(idxm);
        //     if (alNumA.containsKey(c)) { alNumA.put(c, alNumA.get(c) + 1); }
        //     else { alNumA.put('n', alNumA.get('n') + 1); }
        // }

        // HashMap<Character, Integer> alNumB = new HashMap<>();
        // for (char c : alphabet) { alNumB.put(c, 0); }
        // if (alNumB.containsKey('n') == false) { alNumB.put('n', 0); }
        // if (alNumB.containsKey('-') == false) { alNumB.put('-', 0); }
        // for (int i = 0; i < numB; i++) {
        //     char c = this.B[i].charAt(idxn);
        //     if (alNumB.containsKey(c)) { alNumB.put(c, alNumB.get(c) + 1); } 
        //     else { alNumB.put('n', alNumB.get('n') + 1); }
        // }

        // int results = 0, allB = 0;
        // for (char c : alphabet) { allB += alNumB.get(c); }
        // allB += alNumB.get('n');
        // allB += alNumB.get('-');

        // for (char c : alphabet) {
        //     int tempA = alNumA.get(c);
        //     int tempB = alNumB.get(c);
        //     int tempB_ = alNumB.get('-');
        //     results += tempA * ((ms - mis) * tempB + mis * allB - (mis + this.e) * tempB_);
        // }
        // results += alNumA.get('n') * allB * ms;
        // results -= alNumA.get('-') * (allB - alNumB.get('-')) * this.e;
        // results /= (numA * numB);

        // return results;
        */
    }

    @Override
    protected void TraceBack(int k) {
        int m = this.A[0].length(), n = this.B[0].length();
        int diff = n - m;

        int i = m, bj = n, j = diff + k;
        int channel = ChooseMax(pm[0][i][j], pm[1][i][j], pm[2][i][j]);

        while (i > 0 || j > k) {
            if (channel == 0 && i > 0 && j >= 0) {
                if (pm[0][i][j] == pm[0][i-1][j] + Match(i-1, bj-1) && (i > 1) && (j >= 0)) {
                    channel = 0;
                }
                else if (pm[0][i][j] == pm[2][i-1][j] + Match(i-1, bj-1) && i > 1) {
                    channel = 2;
                }
                else if (pm[0][i][j] == pm[1][i-1][j] + Match(i-1, bj-1) && j > 0) {
                    channel = 1;
                }
                for (int idxA = 0; idxA < numA; idxA++) {
                    this.alignA[idxA] += this.A[idxA].charAt(i-1);
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    this.alignB[idxB] += this.B[idxB].charAt(bj-1);
                }
                i -= 1;
                bj -= 1;
            }

            else if (channel == 1 && j > 0) {
                if (pm[1][i][j] == pm[1][i][j-1] - e) {
                    channel = 1;
                }
                else if (pm[1][i][j] == pm[0][i][j-1] - d && i >= 1) {
                    channel = 0;
                }
                for (int idxA = 0; idxA < numA; idxA++) {
                    this.alignA[idxA] += "-";
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    this.alignB[idxB] += this.B[idxB].charAt(bj-1);
                }
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
                for (int idxA = 0; idxA < numA; idxA++) {
                    this.alignA[idxA] += this.A[idxA].charAt(i-1);
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    this.alignB[idxB] += "-";
                }
                i -= 1;
                j += 1;
            }
            else {
                System.out.println("Trace Back is wrong!");
                System.exit(0);
            }
        }

        for (int idxA = 0; idxA < numA; idxA++) {
            this.alignA[idxA] = new StringBuilder(this.alignA[idxA]).reverse().toString();
        }
        for (int idxB = 0; idxB < numB; idxB++) {
            this.alignB[idxB] = new StringBuilder(this.alignB[idxB]).reverse().toString();
        }

        if (this.state) {
            String[] temp = this.alignA;
            this.alignA = this.alignB;
            this.alignB = temp;
        }
    }

    @Override
    protected int Align() {
        long startTime = System.currentTimeMillis();
        int m = this.A[0].length(), n = this.B[0].length();
        int diff = n - m, k = 1;

        if (m == 0 && n == 0) { return 0; }
        else if (m == 0) {
            for (int i = 0; i < this.alignA.length; i++) {
                this.alignA[i] = "-".repeat(n);
            }
            this.alignB = this.B;
            return 0;
        }
        else if (n == 0) {
            for (int i = 0; i < this.alignB.length; i++) {
                this.alignB[i] = "-".repeat(m);
            }
            this.alignA = this.A;
            return 0;
        }

        float valueOld = Float.NEGATIVE_INFINITY, valueNew;
        float[][][] pm = new float[3][m+1][diff+2*k+1];

        int maxk = Math.min(m, Math.max(m/5, 10));
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
            if ( (int) valueNew == (int) valueOld ) {break;}
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
        this.pm = pm;
        TraceBack(k);
        long endTime = System.currentTimeMillis();
        return (int) (startTime - endTime);
    }
}