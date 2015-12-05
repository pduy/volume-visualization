/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;

import java.awt.image.BufferedImage;
import java.lang.reflect.Array;
import java.util.ArrayList;
import javax.media.opengl.GL;
import javax.media.opengl.GL2;

import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {
    private static int BOUNDING_GRADIENT_THRESHOLD = 15;

    private Volume volume = null;
    private GradientVolume gradients = null;
    public String rendererType = "slicer";
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());

        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();


        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());

        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);
//        alphas = computeAlpha();

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public Volume getVolume() {
        return volume;
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }

    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }


    short getVoxel(double[] coord) {

        if (coord[0] <= 0 || coord[0] >= volume.getDimX() || coord[1] <= 0 || coord[1] >= volume.getDimY()
                || coord[2] <= 0 || coord[2] >= volume.getDimZ()) {
            return 0;
        }

        return volume.getVoxel(coord[0], coord[1], coord[2]);
    }

    //MIP approach by brute force
    int getValueByBruteForce(int i, int j, double[] viewVec, double[] uVec, double[] vVec, double[] volumeCenter, int imageCenter, int nLoops) {
        int maxVal = 0;
        double[] pixelCoord = new double[3];
        for (int t = 0; t < nLoops; ++t) {

            double[] currentPosition = new double[3];
            VectorMath.setVector(currentPosition, t * viewVec[0] + volumeCenter[0], t * viewVec[1] + volumeCenter[1], t * viewVec[2] + volumeCenter[2]);

            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                    + currentPosition[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                    + currentPosition[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                    + currentPosition[2];

            int val = getVoxel(pixelCoord);
//            System.out.println(maxVal);
            maxVal = val > maxVal ? val : maxVal;
        }

        return maxVal;
    }

    int getValueByFindingIntersections(int i, int j, double[] viewVec, double[] uVec, double[] vVec, double[] volumeCenter, int imageCenter) {
        int maxVal = 0;
        double[] pixelCoord = new double[3];

        int[] traversalRange = getTraversalRange(i, j, viewVec, uVec, vVec, imageCenter);

        //swap the 2 points if they are not in the correct order
        if (traversalRange[0] > traversalRange[1]) {
            int temp = traversalRange[0];
            traversalRange[0] = traversalRange[1];
            traversalRange[1] = temp;

        }

        for (int t = traversalRange[0]; t <= traversalRange[1]; ++t) {

            double[] currentPosition = new double[3];
            VectorMath.setVector(currentPosition, t * viewVec[0] + volumeCenter[0], t * viewVec[1] + volumeCenter[1], t * viewVec[2] + volumeCenter[2]);

            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                    + currentPosition[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                    + currentPosition[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                    + currentPosition[2];

            int val = getVoxel(pixelCoord);
//            System.out.println(maxVal);
            maxVal = val > maxVal ? val : maxVal;
        }

        return maxVal;
    }

    TFColor getColorByCompositing(int i, int j, double[] viewVec, double[] uVec, double[] vVec, double[] volumeCenter, int imageCenter) {
        double[] pixelCoord = new double[3];

        int[] traversalRange = getTraversalRange(i, j, viewVec, uVec, vVec, imageCenter);

        //swap the 2 points if they are not in the correct order
        if (traversalRange[0] > traversalRange[1]) {
            int temp = traversalRange[0];
            traversalRange[0] = traversalRange[1];
            traversalRange[1] = temp;

        }

        TFColor tempColor = new TFColor(0.0, 0.0, 0.0, 1.0);
        ArrayList<VoxelGradient> boundaryGradCandidates = new ArrayList<>();
        for (int t = traversalRange[1]; t >= traversalRange[0]; --t) {

            //Get the current pixel
            double[] currentPosition = new double[3];
            VectorMath.setVector(currentPosition, t * viewVec[0] + volumeCenter[0], t * viewVec[1] + volumeCenter[1], t * viewVec[2] + volumeCenter[2]);

            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                    + currentPosition[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                    + currentPosition[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                    + currentPosition[2];

            //Extract color
            int val = getVoxel(pixelCoord);
            TFColor pixelColor = tFunc.getColor(val);
            VoxelGradient currGradient = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);

            tempColor.a = 1 - pixelColor.a * tempColor.a;
            tempColor.r = pixelColor.a * pixelColor.r + (1 - pixelColor.a) * tempColor.r;
            tempColor.g = pixelColor.a * pixelColor.g + (1 - pixelColor.a) * tempColor.g;
            tempColor.b = pixelColor.a * pixelColor.b + (1 - pixelColor.a) * tempColor.b;

            if (currGradient.mag >= BOUNDING_GRADIENT_THRESHOLD) {
                boundaryGradCandidates.add(new VoxelGradient(currGradient.x, currGradient.y, currGradient.z));
            }
        }

        if (getVolume().shading && boundaryGradCandidates.size() > 0)  {
            double oldAlpha = tempColor.a;
            tempColor = shade(boundaryGradCandidates.get(boundaryGradCandidates.size() - 1), tempColor, viewVec, 0.1, 0.7, 0.2, 10);
            tempColor.a = oldAlpha;
        }

        return tempColor;

    }

    TFColor getColorBy2DTransferFunction(int i, int j, double[] viewVec, double[] uVec, double[] vVec, double[] volumeCenter, int imageCenter) {
        double[] pixelCoord = new double[3];
        VoxelGradient maximumGradient = new VoxelGradient();

        int[] traversalRange = getTraversalRange(i, j, viewVec, uVec, vVec, imageCenter);

        //swap the 2 points if they are not in the correct order
        if (traversalRange[0] > traversalRange[1]) {
            int temp = traversalRange[0];
            traversalRange[0] = traversalRange[1];
            traversalRange[1] = temp;

        }

        TFColor tempColor = new TFColor(0.0, 0.0, 0.0, 1.0);
        ArrayList<VoxelGradient> boundaryGradCandidates = new ArrayList<>();
        for (int t = traversalRange[1]; t >= traversalRange[0]; --t) {
//        for (int t = 100; t >= -200; --t) {

            //Get the current pixel
            double[] currentPosition = new double[3];
            VectorMath.setVector(currentPosition, t * viewVec[0] + volumeCenter[0], t * viewVec[1] + volumeCenter[1], t * viewVec[2] + volumeCenter[2]);

            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                    + currentPosition[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                    + currentPosition[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                    + currentPosition[2];

            //Extract color
            TFColor pixelColor = tfEditor2D.triangleWidget.color;

            double alpha = getWeightedAlpha(pixelCoord);
            VoxelGradient currGradient = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);
            if (currGradient.mag > maximumGradient.mag && t < traversalRange[1] / 2) {
                maximumGradient = currGradient;
            }

            tempColor.a = 1 - alpha * tempColor.a;
            tempColor.r = alpha * pixelColor.r + (1 - alpha) * tempColor.r;
            tempColor.g = alpha * pixelColor.g + (1 - alpha) * tempColor.g;
            tempColor.b = alpha * pixelColor.b + (1 - alpha) * tempColor.b;

            if (currGradient.mag >= BOUNDING_GRADIENT_THRESHOLD) {
                boundaryGradCandidates.add(new VoxelGradient(currGradient.x, currGradient.y, currGradient.z));
            }
        }

        if (getVolume().shading && boundaryGradCandidates.size() > 0) {
            double oldAlpha = tempColor.a;
            tempColor = shade(boundaryGradCandidates.get(boundaryGradCandidates.size() - 1), tempColor, viewVec, 0.1, 0.7, 0.2, 10);
            tempColor.a = oldAlpha;
        }

        return tempColor;

    }

    double getWeightedAlpha(double[] pixelCoord) {
        double alpha = 0;
        double r = tfEditor2D.triangleWidget.radius;
        int fv = tfEditor2D.triangleWidget.baseIntensity;
        double av = tfEditor2D.triangleWidget.color.a;
        double fi = getVoxel(pixelCoord);
        float deltaFi = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]).mag;

        if (deltaFi == 0 && fi == fv) {
            alpha = av;
        } else if (deltaFi > 0 && fi - r * deltaFi <= fv && fi + r * deltaFi >= fv) {
            alpha = av * (1 - Math.abs(fv - fi) / (r * deltaFi));
        }

        return alpha;
    }

    TFColor shade(VoxelGradient gradient, TFColor pixelColor, double[] viewVector, double kAmbient, double kDiff, double kSpec, double alpha) {

        viewVector[0] /= VectorMath.length(viewVector);
        viewVector[1] /= VectorMath.length(viewVector);
        viewVector[2] /= VectorMath.length(viewVector);

        double[] gradientVector = new double[]{gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag};
        double dotProduct = VectorMath.dotproduct(gradientVector, viewVector);

        if (dotProduct < 0) dotProduct = -dotProduct;
        double r = kAmbient
                + pixelColor.r * kDiff * dotProduct
                + kSpec * Math.pow(dotProduct, alpha);
        double g = kAmbient
                + pixelColor.g * kDiff * dotProduct
                + kSpec * Math.pow(dotProduct, alpha);
        double b = kAmbient
                + pixelColor.b * kDiff * dotProduct
                + kSpec * Math.pow(dotProduct, alpha);
        return new TFColor(r, g, b, 1);
    }

    int[] getTraversalRange(int i, int j, double[] viewVec, double[] uVec, double[] vVec, int imageCenter) {
        double[] pixelCoord = new double[3];

        //Getting the original point v_0
        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter);
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter);
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter);

        //Find 6 intersections
        double[] k = new double[6];

        k[0] = viewVec[0] != 0 ? (-volume.getDimX() / 2 - pixelCoord[0]) / viewVec[0] : Double.MAX_VALUE;
        k[1] = viewVec[1] != 0 ? (-volume.getDimY() / 2 - pixelCoord[1]) / viewVec[1] : Double.MAX_VALUE;
        k[2] = viewVec[2] != 0 ? (-volume.getDimZ() / 2 - pixelCoord[2]) / viewVec[2] : Double.MAX_VALUE;
        k[3] = viewVec[0] != 0 ? ((volume.getDimX() / 2 - pixelCoord[0]) / viewVec[0]) : Double.MAX_VALUE;
        k[4] = viewVec[1] != 0 ? ((volume.getDimY() / 2 - pixelCoord[1]) / viewVec[1]) : Double.MAX_VALUE;
        k[5] = viewVec[2] != 0 ? ((volume.getDimZ() / 2 - pixelCoord[2]) / viewVec[2]) : Double.MAX_VALUE;

        //Check for the valid intersections (which are inside the volume)
        int[] intersections = new int[2];
        int count = 0;
        for (int t = 0; t < 6; ++t) {
            if (k[t] == Double.MAX_VALUE) continue;

            double[] p = new double[3];
            p[0] = pixelCoord[0] + k[t] * viewVec[0];
            p[1] = pixelCoord[1] + k[t] * viewVec[1];
            p[2] = pixelCoord[2] + k[t] * viewVec[2];

            if (p[0] >= -volume.getDimX() / 2 && p[0] <= volume.getDimX() / 2 &&
                    p[1] >= -volume.getDimY() / 2 && p[1] <= volume.getDimY() / 2 &&
                    p[2] >= -volume.getDimZ() / 2 && p[2] <= volume.getDimZ() / 2) {
                intersections[count] = ((int) Math.ceil(k[t]));
                count++;
            }
        }

        return intersections;
    }

    //Original method (Center Slice)
    int getValueByCenter(int i, int j, double[] uVec, double[] vVec, double[] volumeCenter, int imageCenter) {
        double[] pixelCoord = new double[3];

        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                + volumeCenter[0];
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                + volumeCenter[1];
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                + volumeCenter[2];

        return getVoxel(pixelCoord);
    }

    void render(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data

        TFColor voxelColor = new TFColor();

        //Calculating max
        double max = volume.getMaximum();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                int val;

                switch (this.rendererType) {
                    case "MIP":
                        //                int val = getValueByBruteForce(i, j, viewVec, uVec, vVec, volumeCenter, imageCenter, diagonal);
                        val = getValueByFindingIntersections(i, j, viewVec, uVec, vVec, volumeCenter, imageCenter);

                        // Map the intensity to a grey value by linear scaling
                        voxelColor.r = val / max;
                        voxelColor.g = voxelColor.r;
                        voxelColor.b = voxelColor.r;
                        voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                        break;
                    case "compositing":
                        voxelColor = getColorByCompositing(i, j, viewVec, uVec, vVec, volumeCenter, imageCenter);
                        break;
                    case "tf2d":
                        voxelColor = getColorBy2DTransferFunction(i, j, viewVec, uVec, vVec, volumeCenter, imageCenter);
                        break;
                    case "slicer":
                    default:
                        val = getValueByCenter(i, j, uVec, vVec, volumeCenter, imageCenter);
                        // Map the intensity to a grey value by linear scaling
                        voxelColor.r = val / max;
                        voxelColor.g = voxelColor.r;
                        voxelColor.b = voxelColor.r;
                        voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                        break;
                }


                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        render(viewMatrix);

        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }

    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
//        System.out.println("Recomputing alphas");
//        computeAlpha();
//        System.out.println("Done");

        for (int i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
