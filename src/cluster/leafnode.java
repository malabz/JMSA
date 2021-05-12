package cluster;

public class leafnode extends node{
    private String name;
    public int number;

    public leafnode (String name, int num) {
        this.name = name;
        this.number = num;
    }

    public int getNum() {
        return number;
    }

    public String toString() {
        return name;
    }

    @Override
    public double getDistance() {
        return 0.0;
    }
}
