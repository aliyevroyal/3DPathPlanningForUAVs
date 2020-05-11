package com.thealiyev;

import java.util.ArrayList;

public class Obstacles {
    private ArrayList<ArrayList<ArrayList<Double>>> obstacles = new ArrayList<>();
    private ArrayList<ArrayList<Double>> obstacle1 = null;
    private ArrayList<ArrayList<Double>> obstacle2 = null;
    private ArrayList<ArrayList<Double>> obstacle3 = null;
    private ArrayList<ArrayList<Double>> obstacle4 = null;
    private ArrayList<ArrayList<Double>> obstacle5 = null;
    private ArrayList<ArrayList<Double>> obstacle6 = null;
    private ArrayList<ArrayList<Double>> obstacle7 = null;
    private ArrayList<ArrayList<Double>> obstacle8 = null;
    private ArrayList<ArrayList<Double>> obstacle9 = null;
    private ArrayList<ArrayList<Double>> obstacle10 = null;
    private ArrayList<ArrayList<Double>> obstacle11 = null;
    private ArrayList<ArrayList<Double>> obstacle12 = null;
    private ArrayList<ArrayList<Double>> obstacle13 = null;
    private ArrayList<ArrayList<Double>> obstacle14 = null;
    private ArrayList<ArrayList<Double>> obstacle15 = null;
    private ArrayList<ArrayList<Double>> obstacle16 = null;

    private ArrayList<ArrayList<Double>> obstacle17 = null;
    private ArrayList<ArrayList<Double>> obstacle18 = null;
    private ArrayList<ArrayList<Double>> obstacle19 = null;
    private ArrayList<ArrayList<Double>> obstacle20 = null;
    private ArrayList<ArrayList<Double>> obstacle21 = null;
    private ArrayList<ArrayList<Double>> obstacle22 = null;
    private ArrayList<ArrayList<Double>> obstacle23 = null;
    private ArrayList<ArrayList<Double>> obstacle24 = null;

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
        startCoordinate.add(5.0);
        startCoordinate.add(7.5);
        startCoordinate.add(0.0);
        obstacle1.add(startCoordinate);

        endCoordinate.add(10.0);
        endCoordinate.add(20.0);
        endCoordinate.add(100.0);
        obstacle1.add(endCoordinate);

        obstacles.add(obstacle1);
    }

    public void setObstacle2() {
        this.obstacle2 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(20.0);
        startCoordinate.add(5.0);
        startCoordinate.add(0.0);
        obstacle2.add(startCoordinate);

        endCoordinate.add(44.0);
        endCoordinate.add(44.0);
        endCoordinate.add(100.0);
        obstacle2.add(endCoordinate);

        obstacles.add(obstacle2);
    }

    public void setObstacle3() {
        this.obstacle3 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(8.0);
        startCoordinate.add(30.0);
        startCoordinate.add(0.0);
        obstacle3.add(startCoordinate);

        endCoordinate.add(10.0);
        endCoordinate.add(34.0);
        endCoordinate.add(100.0);
        obstacle3.add(endCoordinate);

        obstacles.add(obstacle3);
    }

    public void setObstacle4() {
        this.obstacle4 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(17.0);
        startCoordinate.add(14.0);
        startCoordinate.add(0.0);
        obstacle4.add(startCoordinate);

        endCoordinate.add(19.0);
        endCoordinate.add(16.0);
        endCoordinate.add(100.0);
        obstacle4.add(endCoordinate);

        obstacles.add(obstacle4);
    }

    public void setObstacle5() {
        this.obstacle5 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(85.0);
        startCoordinate.add(85.0);
        startCoordinate.add(0.0);
        obstacle5.add(startCoordinate);

        endCoordinate.add(90.0);
        endCoordinate.add(90.0);
        endCoordinate.add(100.0);
        obstacle5.add(endCoordinate);

        obstacles.add(obstacle5);
    }

    public void setObstacle6() {
        this.obstacle6 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(91.0);
        startCoordinate.add(91.0);
        startCoordinate.add(0.0);
        obstacle6.add(startCoordinate);

        endCoordinate.add(95.0);
        endCoordinate.add(95.0);
        endCoordinate.add(100.0);
        obstacle6.add(endCoordinate);

        obstacles.add(obstacle6);
    }

    public void setObstacle7() {
        this.obstacle7 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(96.0);
        startCoordinate.add(96.0);
        startCoordinate.add(0.0);
        obstacle7.add(startCoordinate);

        endCoordinate.add(98.0);
        endCoordinate.add(98.0);
        endCoordinate.add(100.0);
        obstacle7.add(endCoordinate);

        obstacles.add(obstacle7);
    }

    public void setObstacle8() {
        this.obstacle8 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(60.0);
        startCoordinate.add(15.0);
        startCoordinate.add(0.0);
        obstacle8.add(startCoordinate);

        endCoordinate.add(75.0);
        endCoordinate.add(40.0);
        endCoordinate.add(100.0);
        obstacle8.add(endCoordinate);

        obstacles.add(obstacle8);
    }

    public void setObstacle17() {
        this.obstacle17 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(0.0);
        startCoordinate.add(2.0);
        startCoordinate.add(0.0);
        obstacle17.add(startCoordinate);

        endCoordinate.add(10.0);
        endCoordinate.add(2.5);
        endCoordinate.add(1.5);
        obstacle17.add(endCoordinate);

        obstacles.add(obstacle17);
    }

    public void setObstacle18() {
        this.obstacle18 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(0.0);
        startCoordinate.add(2.0);
        startCoordinate.add(4.5);
        obstacle18.add(startCoordinate);

        endCoordinate.add(10.0);
        endCoordinate.add(2.5);
        endCoordinate.add(6.0);
        obstacle18.add(endCoordinate);

        obstacles.add(obstacle18);
    }

    public void setObstacle19() {
        this.obstacle19 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(0.0);
        startCoordinate.add(2.0);
        startCoordinate.add(1.5);
        obstacle19.add(startCoordinate);

        endCoordinate.add(3.0);
        endCoordinate.add(2.5);
        endCoordinate.add(4.5);
        obstacle19.add(endCoordinate);

        obstacles.add(obstacle19);
    }

    public void setObstacle20() {
        this.obstacle20 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(7.0);
        startCoordinate.add(2.0);
        startCoordinate.add(1.5);
        obstacle20.add(startCoordinate);

        endCoordinate.add(10.0);
        endCoordinate.add(2.5);
        endCoordinate.add(4.5);
        obstacle20.add(endCoordinate);

        obstacles.add(obstacle20);
    }

    public void setObstacle21() {
        this.obstacle21 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(3.0);
        startCoordinate.add(0.0);
        startCoordinate.add(2.4);
        obstacle21.add(startCoordinate);

        endCoordinate.add(7.0);
        endCoordinate.add(0.5);
        endCoordinate.add(4.5);
        obstacle21.add(endCoordinate);

        obstacles.add(obstacle21);
    }

    public void setObstacle22() {
        this.obstacle22 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(0.0);
        startCoordinate.add(15.0);
        startCoordinate.add(0.0);
        obstacle22.add(startCoordinate);

        endCoordinate.add(10.0);
        endCoordinate.add(20.0);
        endCoordinate.add(1.0);
        obstacle22.add(endCoordinate);

        obstacles.add(obstacle22);
    }

    public void setObstacle23() {
        this.obstacle23 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(0.0);
        startCoordinate.add(15.0);
        startCoordinate.add(1.0);
        obstacle23.add(startCoordinate);

        endCoordinate.add(10.0);
        endCoordinate.add(16.0);
        endCoordinate.add(3.5);
        obstacle23.add(endCoordinate);

        obstacles.add(obstacle23);
    }

    public void setObstacle24() {
        this.obstacle24 = new ArrayList<>();
        ArrayList<Double> startCoordinate = new ArrayList<>();
        ArrayList<Double> endCoordinate = new ArrayList<>();
        startCoordinate.add(0.0);
        startCoordinate.add(18.0);
        startCoordinate.add(4.5);
        obstacle24.add(startCoordinate);

        endCoordinate.add(10.0);
        endCoordinate.add(19.0);
        endCoordinate.add(6.0);
        obstacle24.add(endCoordinate);

        obstacles.add(obstacle24);
    }

    public ArrayList<ArrayList<ArrayList<Double>>> getObstacles() {
        return obstacles;
    }

    public ArrayList<ArrayList<ArrayList<Double>>> getSortedObstacles(ArrayList<Double> station) {
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> duplicatedObstacles = new ArrayList<>(obstacles);
        ArrayList<ArrayList<Double>> obstacle;
        double distance;

        for (int stCounter = 0; stCounter < obstacles.size(); stCounter = stCounter + 1) {
            distance = obstacleAvoider.findEuclidean2DDistance(station, findCenterOfObstacle(create3DObstacle(duplicatedObstacles.get(0))));
            obstacle = duplicatedObstacles.get(0);
            for (int ndCounter = 0; ndCounter < duplicatedObstacles.size(); ndCounter = ndCounter + 1) {
                if (obstacleAvoider.findEuclidean2DDistance(station, findCenterOfObstacle(create3DObstacle(duplicatedObstacles.get(ndCounter)))) < distance) {
                    distance = obstacleAvoider.findEuclidean2DDistance(station, findCenterOfObstacle(create3DObstacle(duplicatedObstacles.get(ndCounter))));
                    obstacle = duplicatedObstacles.get(ndCounter);
                }
            }
            sortedObstacles.add(obstacle);
            duplicatedObstacles.remove(new ArrayList<>(obstacle));
        }

        return sortedObstacles;
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