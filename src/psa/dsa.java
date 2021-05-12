package psa;

public class dsa {
    private String A;
    private String B;
    private String mode;
    private String[] alignAB;

    /**
     * mode : "kand" or "suffix"
     * @param A
     * @param B
     * @param mode
     */
    public dsa (String A, String B, String mode) {
        this.A = A;
        this.B = B;
        this.mode = mode;
        Align();
    }

    /**
     * @param A
     * @param B
     */
    public dsa (String A, String B) {
        this.A = A;
        this.B = B;
        this.mode = "suffix";
        Align();
    }

    /**
     * To get the align results.
     * @return
     */
    public String[] getStrAlign() {
        return alignAB;
    }

    private int Align() {
        if (this.A.length() == 0 || this.B.length() == 0) {
            this.mode = "kband";
        }
        if (mode == "kband") {
            Kband kalign = new Kband(A, B);
            alignAB = kalign.getStrAlign();
        }
        else if (mode == "suffix") {
            if (this.B.length() > this.A.length()) {
                STAlign salign = new STAlign(B, A);
                alignAB = new String[2];
                alignAB[0] = salign.getStrAlign()[1];
                alignAB[1] = salign.getStrAlign()[0];
            }
            else {
                STAlign salign = new STAlign(A, B);
                alignAB = salign.getStrAlign();
            }
        }
        else {
            System.out.println("The mode is wrong!");
            return -1;
        }
        return 1;
    }
}
