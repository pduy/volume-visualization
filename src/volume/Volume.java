/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import util.VectorMath;

import java.io.File;
import java.io.IOException;

/**
 * @author michel
 */
public class Volume {

    public Volume(int xd, int yd, int zd) {
        data = new short[xd * yd * zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }

    public Volume(File file) {

        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }

    }


    public short getVoxel(double x, double y, double z) {
        int xFloor = (int) Math.floor(x);
        int xCeiling = (int) Math.ceil(x);
        int yFloor = (int) Math.floor(y);
        int yCeiling = (int) Math.ceil(y);
        int zFloor = (int) Math.floor(z);
        int zCeiling = (int) Math.ceil(z);

        //Find 8 corner points
        double[] x0 = new double[3];
        VectorMath.setVector(x0, xFloor, yFloor, zFloor);
        double[] x1 = new double[3];
        VectorMath.setVector(x1, xCeiling, yFloor, zFloor);
        double[] x2 = new double[3];
        VectorMath.setVector(x2, xFloor, yFloor, zCeiling);
        double[] x3 = new double[3];
        VectorMath.setVector(x3, xCeiling, xFloor, zCeiling);
        double[] x4 = new double[3];
        VectorMath.setVector(x4, xFloor, yCeiling, zFloor);
        double[] x5 = new double[3];
        VectorMath.setVector(x5, xCeiling, yCeiling, zFloor);
        double[] x6 = new double[3];
        VectorMath.setVector(x6, xFloor, yCeiling, zCeiling);
        double[] x7 = new double[3];
        VectorMath.setVector(x7, xCeiling, yCeiling, zCeiling);

        //Alpha, beta. gamma
        double alpha = Math.abs((x - xFloor) / (x - xCeiling));
        double beta = Math.abs((z - zFloor) / (z - zCeiling));
        double gamma = Math.abs((y - yFloor) / (y - yCeiling));

        //Find 8 intensities
        double sX0 = 1.0 * data[(int)x0[0] + dimX * ((int)x0[1] + dimY * (int)x0[2])];
        double sX1 = 1.0 * data[(int)x1[0] + dimX * ((int)x1[1] + dimY * (int)x1[2])];
        double sX2 = 1.0 * data[(int)x2[0] + dimX * ((int)x2[1] + dimY * (int)x2[2])];
        double sX3 = 1.0 * data[(int)x3[0] + dimX * ((int)x3[1] + dimY * (int)x3[2])];
        double sX4 = 1.0 * data[(int)x4[0] + dimX * ((int)x4[1] + dimY * (int)x4[2])];
        double sX5 = 1.0 * data[(int)x5[0] + dimX * ((int)x5[1] + dimY * (int)x5[2])];
        double sX6 = 1.0 * data[(int)x6[0] + dimX * ((int)x6[1] + dimY * (int)x6[2])];
        double sX7 = 1.0 * data[(int)x7[0] + dimX * ((int)x7[1] + dimY * (int)x7[2])];


        return (short)((1 - alpha) * (1 - beta) * (1 - gamma) * sX0 +
                alpha * (1 - beta) * (1 - gamma) * sX1 +
                (1 - alpha) * beta * (1 - gamma) * sX2 +
                alpha * beta * (1 - gamma) * sX3 +
                (1 - alpha) * (1 - beta) * gamma * sX4 +
                alpha * (1 - beta) * gamma * sX5 +
                (1 - alpha) * beta * gamma * sX6 +
                alpha * beta * gamma * sX7);
    }

    public void setVoxel(int x, int y, int z, short value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, short value) {
        data[i] = value;
    }

    public short getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    public short getMinimum() {
        short minimum = data[0];
        for (int i = 0; i < data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }

    public short getMaximum() {
        short maximum = data[0];
        for (int i = 0; i < data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }

    public int[] getHistogram() {
        return histogram;
    }

    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i = 0; i < data.length; i++) {
            histogram[data[i]]++;
        }
    }

    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;
}
