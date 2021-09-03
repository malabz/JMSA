package hierCluster;

public class midnode extends node {
    private final node a, b;
    private int number;

    public midnode (node a, node b, int num) {
        this.a = a;
        this.b = b;
        this.number = num;
    }

    public int getNum() {
        return number;
    }

    public String toString() {
        return "(" + a + ":" + a.getLen() + ", " + b + ":" + b.getLen() + ")";
        // return "(" + a + ", " + b  + ")";

    }

    @Override
    public double getDistance() {
        return a.getLen() + a.getDistance();
    } 
}
