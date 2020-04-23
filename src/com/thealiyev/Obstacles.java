package com.thealiyev;

import java.util.ArrayList;

public class Obstacles {
    private ArrayList<ArrayList<ArrayList<Double>>> obstacles = new ArrayList<>();
    private ArrayList<ArrayList<Double>> obstacle1 = null;

    public Obstacles() {
    }

    public double getR(ArrayList<ArrayList<Double>> obstacle) {
        double R;
        R = findEuclideanDistance(obstacle.get(0), obstacle.get(2));
        R = R / 2;
        return R;
    }

    public void setObstacle1() {
        this.obstacle1 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(10.0);
        startCoordinate.add(10.0);
        startCoordinate.add(0.0);
        obstacle1.add(startCoordinate);

        endCoordinate.add(20.0);
        endCoordinate.add(20.0);
        endCoordinate.add(1.5);
        obstacle1.add(endCoordinate);

        obstacles.add(obstacle1);
    }

    public ArrayList<ArrayList<ArrayList<Double>>> getObstacles() {
        return obstacles;
    }

    public ArrayList<ArrayList<Double>> create3DObstacle(ArrayList<ArrayList<Double>> obstacle) {
        ArrayList<ArrayList<Double>> threeDimensionalObstacle = new ArrayList<>();
        ArrayList<Double> stCoordinate = new ArrayList<>();
        ArrayList<Double> ndCoordinate = new ArrayList<>();
        ArrayList<Double> rdCoordinate = new ArrayList<>();
        ArrayList<Double> fourthCoordinate = new ArrayList<>();
        ArrayList<Double> zCoorditanes = new ArrayList<>();
        double x, y, z;
        //st coordinate
        x = obstacle.get(0).get(0);
        y = obstacle.get(0).get(1);
        stCoordinate.add(x);
        stCoordinate.add(y);
        threeDimensionalObstacle.add(stCoordinate);
        //nd coordinate
        x = obstacle.get(0).get(0);
        y = obstacle.get(1).get(1);
        ndCoordinate.add(x);
        ndCoordinate.add(y);
        threeDimensionalObstacle.add(ndCoordinate);
        //rd coordinate
        x = obstacle.get(1).get(0);
        y = obstacle.get(1).get(1);
        rdCoordinate.add(x);
        rdCoordinate.add(y);
        threeDimensionalObstacle.add(rdCoordinate);
        //fourth coordinate
        x = obstacle.get(1).get(0);
        y = obstacle.get(0).get(1);
        fourthCoordinate.add(x);
        fourthCoordinate.add(y);
        threeDimensionalObstacle.add(fourthCoordinate);
        //z coordinates
        z = obstacle.get(0).get(2);
        zCoorditanes.add(z);
        z = obstacle.get(1).get(2);
        zCoorditanes.add(z);
        threeDimensionalObstacle.add(zCoorditanes);

        return threeDimensionalObstacle;
    }

    public ArrayList<Double> findCenterOfObstacle(ArrayList<ArrayList<Double>> obstacle) {
        ArrayList<Double> centerCoordinate = new ArrayList<>();
        double d, x = obstacle.get(0).get(0), y = obstacle.get(0).get(1), z;
        double cosAlpha, alphaDegree;

        d = findEuclideanDistance(obstacle.get(0), obstacle.get(2));
        cosAlpha = (obstacle.get(2).get(0) - obstacle.get(0).get(0)) / d;
        alphaDegree = Math.toDegrees(Math.acos(cosAlpha));

        d = d / 2;
        x = x + (d * Math.cos(Math.toRadians(alphaDegree)));
        y = y + (d * Math.sin(Math.toRadians(alphaDegree)));
        z = obstacle.get(4).get(1) - obstacle.get(4).get(0);
        z = z / 2;

        centerCoordinate.add(x);
        centerCoordinate.add(y);
        centerCoordinate.add(z);

        return centerCoordinate;
    }

    private static double findEuclideanDistance(ArrayList<Double> firstPoint, ArrayList<Double> secondPoint) {
        double euclideanDistance = 0.0;

        for (int stCounter = 0; stCounter < secondPoint.size(); stCounter = stCounter + 1) {
            euclideanDistance = euclideanDistance + Math.pow((secondPoint.get(stCounter) - firstPoint.get(stCounter)), 2);
        }

        euclideanDistance = Math.sqrt(euclideanDistance);

        return euclideanDistance;
    }
}