package cluster;

public abstract class node {
    private double branchlength;

    public abstract int getNum();
    public abstract double getDistance();
    public double getLen() {return branchlength;}
    public void setLen(double length) {branchlength = length;}
}
