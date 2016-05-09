package cz.sange;

/**
 * Created by sange on 25/04/16.
 */
public class Solution {

    private final int numWarehouses;
    private final int numCustomers;
    private int status;
    private double objVal;
    private int [][] y;
    private double [][] x;
    private double [][] w;
    private double [] u;
    private double [] v;

    public Solution(int numWarehouses, int numCustomers) {
        this.numWarehouses = numWarehouses;
        this.numCustomers = numCustomers;

        y = new int[numWarehouses][numCustomers];
        x = new double[numWarehouses][numCustomers];
        w = new double[numWarehouses][numCustomers];
        u = new double[numWarehouses];
        v = new double[numCustomers];
    }

    public int getStatus() {
        return status;
    }

    public void setStatus(int status) {
        this.status = status;
    }

    public double getObjVal() {
        return objVal;
    }

    public void setObjVal(double objVal) {
        this.objVal = objVal;
    }

    public int[][] getY() {
        return y;
    }

    public void setY(int[][] y) {
        this.y = y;
    }

    public double[][] getX() {
        return x;
    }

    public void setX(double[][] x) {
        this.x = x;
    }

    public double[][] getW() {
        return w;
    }

    public void setW(double[][] w) {
        this.w = w;
    }

    public double[] getU() {
        return u;
    }

    public void setU(double[] u) {
        this.u = u;
    }

    public double[] getV() {
        return v;
    }

    public void setV(double[] v) {
        this.v = v;
    }

    public void reset() {
        this.status = 0;
        this.objVal = 0;

        for (int i = 0; i < numWarehouses; i++) {
            u[i] = 0;
            for (int j = 0; j < numCustomers; j++) {
                v[j] = 0;
                x[i][j] = 0;
                y[i][j] = 0;
                w[i][j] = 0;
            }
        }
    }
}
