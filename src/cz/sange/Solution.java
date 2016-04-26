package cz.sange;

/**
 * Created by sange on 25/04/16.
 */
public class Solution {

    private int status;
    private double objVal;
    private double unbdRay;
    private int [][] y;
    private double [][] x;
    private double [][] w;
    private double [] u;
    private double [] v;

    public Solution() {
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

    public double getUnbdRay() {
        return unbdRay;
    }

    public void setUnbdRay(double unbdRay) {
        this.unbdRay = unbdRay;
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
}
