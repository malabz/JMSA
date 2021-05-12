package measure;

/** 
 * Compute score sp
 * 
 * @author Juntao chen
*/
public class ScoreSP {
    
    protected int ms = 7, mis = -3;
    protected int d = 11, e = 2;

    /**
     * To compute the spscore of two aligned strings.
     * @param A
     * @param B
     * @return score
     */
    public double ComputeTwo(String A, String B) {
        int score = 0;
        int gap = 0;
        for (int i = 0; i < A.length(); ++i) {
            if (A.charAt(i) != '-' && B.charAt(i) != '-') {
                if (gap != 0) {gap=0;}
                if (A.charAt(i) == B.charAt(i)) {score += ms;}
                else {score += mis;}
            }
            else if (A.charAt(i) != '-' || B.charAt(i) != '-') {
                if (gap == 0) {
                    score -= d;
                    gap = 1;
                }
                else {score -= e;}
            }
        }
        return (double) score/(A.length() * this.ms);
    }

    /**
     * 
     * @param As
     * @param Bs
     * @return
     */
    public long ComputeN(String[] As, String[] Bs) {
        int score = 0;
        for (int i= 0; i < As.length; ++i) {
            for (int j = 0; j < i; ++j) {
                score += ComputeTwo(As[i], Bs[j]);
            }
        }
        return score;
    }

    /**
     * 
     * @param As
     * @return
     */
    public double ComputeN(String[] As) {
        double score = 0;
        for (int i = 0; i < As.length; ++i) {
            for (int j = 0; j < i; ++j) {
                score += ComputeTwo(As[i], As[j]);
            }
        }
        return score/(As.length * (As.length - 1) / 2);
    }
}
