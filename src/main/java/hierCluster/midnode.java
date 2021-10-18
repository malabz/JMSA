package hierCluster;

public class midnode extends node {
    private final node a, b;
    private final int number;

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
    }

    @Override
    public double getDistance() {
        return a.getLen() + a.getDistance();
    } 
}
