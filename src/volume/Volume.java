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
    public boolean interpolate = false;

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
        if (interpolate) {
            return(getVoxelWithInterpolation(x, y, z));
        } else {
            try {
                // or should I floor x,y,z here?
                return data[(int)x + (int)dimX * ((int)y + (int)dimY * (int)z)];
            } catch (ArrayIndexOutOfBoundsException e) {
                return 0;
            }
        }
    }

    public short getVoxelWithInterpolation(double x, double y, double z) {
        int xFloor = (int) Math.floor(x);
        int xCeiling = (int) Math.ceil(x);
        int yFloor = (int) Math.floor(y);
        int yCeiling = (int) Math.ceil(y);
        int zFloor = (int) Math.floor(z);
        int zCeiling = (int) Math.ceil(z);

        //Alpha, beta. gamma
        double alpha = (xCeiling - xFloor) == 0 ? 1 : (x - xFloor) / (xCeiling - xFloor);
        double beta = (zCeiling - zFloor) == 0 ? 1 : (z - zFloor) / (zCeiling - zFloor);
        double gamma = (yCeiling - yFloor) == 0 ? 1 : (y - yFloor) / (yCeiling - yFloor);

        int yzFloor = yFloor + dimY * zFloor;
        int yFloorZCeiling = yFloor + dimY * zCeiling;
        int yCeilingZFloor = yCeiling + dimY * zFloor;
        int yzCeiling = yCeiling + dimY * zCeiling;

        int sX0Index = xFloor + dimX * yzFloor;
        int sX1Index = xCeiling + dimX * yzFloor;
        int sX2Index = xFloor + dimX * yFloorZCeiling;
        int sX3Index = xCeiling + dimX * yFloorZCeiling;
        int sX4Index = xFloor + dimX * yCeilingZFloor;
        int sX5Index = xCeiling + dimX * yCeilingZFloor;
        int sX6Index = xFloor + dimX * yzCeiling;
        int sX7Index = xCeiling + dimX * yzCeiling;

        if (sX0Index >= data.length ||
                sX1Index >= data.length ||
                sX2Index >= data.length ||
                sX3Index >= data.length ||
                sX4Index >= data.length ||
                sX5Index >= data.length ||
                sX6Index >= data.length ||
                sX7Index >= data.length) {
            return 0;
        }
        //Find 8 intensities
        double sX0 = 1.0 * data[sX0Index];
        double sX1 = 1.0 * data[sX1Index];
        double sX2 = 1.0 * data[sX2Index];
        double sX3 = 1.0 * data[sX3Index];
        double sX4 = 1.0 * data[sX4Index];
        double sX5 = 1.0 * data[sX5Index];
        double sX6 = 1.0 * data[sX6Index];
        double sX7 = 1.0 * data[sX7Index];

        return (short) ((1 - alpha) * (1 - beta) * (1 - gamma) * sX0 +
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
        if (i < 0) return 0;
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
