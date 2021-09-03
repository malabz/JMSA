package psa;

public class dsa {
    private final String A;
    private final String B;
    private String mode;
    private String[] alignAB;

    /**
     * mode : "kand" or "suffix" or "fmindex"
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
     * To get the aligned results.
     */
    public String[] getStrAlign() {
        return alignAB;
    }

    private void Align() {
        if (this.A.length() == 0 || this.B.length() == 0) {
            this.mode = "kband";
        }
        switch (mode) {
            case "kband" -> {
                Kband kb = new Kband(A, B);
                alignAB = kb.getStrAlign();
            }
            case "suffix" -> {
                if (this.B.length() > this.A.length()) {
                    STAlign stAlign = new STAlign(B, A);
                    alignAB = new String[2];
                    alignAB[0] = stAlign.getStrAlign()[1];
                    alignAB[1] = stAlign.getStrAlign()[0];
                } else {
                    STAlign stAlign = new STAlign(A, B);
                    alignAB = stAlign.getStrAlign();
                }
            }
            case "fmindex" -> {
                if (this.B.length() > this.A.length()) {
                    FMAlign fmAlign = new FMAlign(B, A);
                    alignAB = new String[2];
                    alignAB[0] = fmAlign.getStrAlign()[1];
                    alignAB[1] = fmAlign.getStrAlign()[0];
                } else {
                    FMAlign fmAlign = new FMAlign(A, B);
                    alignAB = fmAlign.getStrAlign();
                }
            }
            default -> {
                System.out.println("(dsa) mode is wrong!");
                System.exit(0);
            }
        }
    }
}
