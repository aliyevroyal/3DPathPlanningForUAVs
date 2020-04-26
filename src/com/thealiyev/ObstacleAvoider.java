package com.thealiyev;

import java.util.ArrayList;

public class ObstacleAvoider {

    public static void main(String[] args) {
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<ArrayList<Double>> path = new ArrayList<>();
        ArrayList<Double> position = new ArrayList<>();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        boolean isPointInsideOfObstacle, didPointCollideWithObstacle;
        ArrayList<ArrayList<Double>> positionsMatrixWithCollisions = new ArrayList<>();

        Obstacles obstacles = new Obstacles();
        obstacles.setObstacle1();
        obstacles.setObstacle2();
        obstacles.setObstacle3();
        obstacles.setObstacle4();
        obstacles.setObstacle5();
        obstacles.setObstacle6();
        obstacles.setObstacle7();
        obstacles.setObstacle8();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;

        position.add(45.968439064027976);
        position.add(46.82734025683803);
        path.add(position);
        position = new ArrayList<>();

        position.add(50.0);
        position.add(50.0);
        path.add(position);
        position = new ArrayList<>();

        position.add(1.0);
        position.add(1.0);
        path.add(position);

        for (int stCounter = 0; stCounter < path.size() - 1; stCounter = stCounter + 1) {
            ObstacleAvoidanceCurrentStation = path.get(stCounter);
            ObstacleAvoidanceNextStation = path.get(stCounter + 1);
            positionsMatrixWithCollisions.add(ObstacleAvoidanceCurrentStation);
            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                if (isPointInsideOfObstacle) {
                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                    didPointCollideWithObstacle = obstacleAvoider.didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                    if (didPointCollideWithObstacle) {
                        nearestCornerToCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        positionsMatrixWithCollisions.add(nearestCornerToCurrentStation);
                        newCornerForNextStation = obstacleAvoider.findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                        positionsMatrixWithCollisions.add(newCornerForNextStation);
                        ObstacleAvoidanceCurrentStation = newCornerForNextStation;
                    }
                } else {
                    didPointCollideWithObstacle = obstacleAvoider.didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                    if (didPointCollideWithObstacle) {
                        nearestCornerToCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        positionsMatrixWithCollisions.add(nearestCornerToCurrentStation);
                        newCornerForNextStation = obstacleAvoider.findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                        positionsMatrixWithCollisions.add(newCornerForNextStation);
                        ObstacleAvoidanceCurrentStation = newCornerForNextStation;
                    }
                }
            }
            positionsMatrixWithCollisions.add(ObstacleAvoidanceNextStation);
        }

        System.out.println(positionsMatrixWithCollisions);
    }

    public boolean isPointInsideOfObstacle(ArrayList<Double> p, ArrayList<ArrayList<Double>> obstacle) {
        boolean isPointInsideOfObstacle = false;
        double obstacleXMin, obstacleXMax, obstacleYMin, obstacleYMax;

        obstacleXMin = obstacle.get(0).get(0);
        obstacleXMax = obstacle.get(1).get(0);
        obstacleYMin = obstacle.get(0).get(1);
        obstacleYMax = obstacle.get(1).get(1);

        if (p.get(0) < obstacleXMax && p.get(0) > obstacleXMin) {
            if (p.get(1) < obstacleYMax && p.get(1) > obstacleYMin) {
                isPointInsideOfObstacle = true;
            }
        }

        return isPointInsideOfObstacle;
    }

    public boolean didPointCollideWithObstacle(ArrayList<Double> p1, ArrayList<Double> p2, ArrayList<ArrayList<Double>> obstacle) {
        Obstacles obstacles = new Obstacles();
        boolean didCollide = false;
        ArrayList<Double> O = new ArrayList<>();
        double A, h, R;

        O.add(obstacles.findCenterOfObstacle(obstacles.create3DObstacle(obstacle)).get(0));
        O.add(obstacles.findCenterOfObstacle(obstacles.create3DObstacle(obstacle)).get(1));

        A = calculateTriangleArea(p1, p2, O);
        h = (2 * A) / findEuclidean2DDistance(p2, p1);
        R = obstacles.getR(obstacles.create3DObstacle(obstacle));

        if (h <= R && h >= 0) {
            if (findEuclidean2DDistance(p1, p2) >= findEuclidean2DDistance(O, p2) && findEuclidean2DDistance(p1, p2) >= findEuclidean2DDistance(O, p1)) {
                didCollide = true;
            }
        }

        return didCollide;
    }

    public double calculateTriangleArea(ArrayList<Double> p1, ArrayList<Double> p2, ArrayList<Double> obstacle) {
        double a, b, c;
        double A, P;

        a = findEuclidean2DDistance(obstacle, p1);
        b = findEuclidean2DDistance(obstacle, p2);
        c = findEuclidean2DDistance(p1, p2);

        P = a + b + c;
        P = P / 2;
        A = P * (P - a) * (P - b) * (P - c);
        if (A < 0) {
            A = -1 * A;
        }
        A = Math.sqrt(A);

        return A;
    }

    public double findEuclidean2DDistance(ArrayList<Double> firstPoint, ArrayList<Double> secondPoint) {
        double euclideanDistance = 0.0;

        for (int stCounter = 0; stCounter < 2; stCounter = stCounter + 1) {
            euclideanDistance = euclideanDistance + Math.pow((secondPoint.get(stCounter) - firstPoint.get(stCounter)), 2);
        }

        euclideanDistance = Math.sqrt(euclideanDistance);

        return euclideanDistance;
    }

    public ArrayList<Double> findNearestCorner(ArrayList<Double> p, ArrayList<ArrayList<Double>> obstacle) {
        Obstacles obstacles = new Obstacles();
        ArrayList<Double> nearestPoint;
        double d = findEuclidean2DDistance(p, obstacles.create3DObstacle(obstacle).get(0));
        int index = 0;

        for (int counter = 0; counter < obstacles.create3DObstacle(obstacle).size() - 1; counter = counter + 1) {
            if (findEuclidean2DDistance(p, obstacles.create3DObstacle(obstacle).get(counter)) < d) {
                d = findEuclidean2DDistance(p, obstacles.create3DObstacle(obstacle).get(counter));
                index = counter;
            }
        }
        nearestPoint = obstacles.create3DObstacle(obstacle).get(index);

        return nearestPoint;
    }

    public ArrayList<Double> findPathToOppositeSide(ArrayList<Double> firstPoint, ArrayList<Double> secondPoint, ArrayList<ArrayList<Double>> obstacle) {
        Obstacles obstacles = new Obstacles();
        ArrayList<Double> nearestCornerToP1 = findNearestCorner(firstPoint, obstacle);
        int nearestCornersIndexToP1 = obstacles.create3DObstacle(obstacle).indexOf(nearestCornerToP1);
        double d1, d2;
        ArrayList<Double> newCornerForP2 = new ArrayList<>();

        if (nearestCornersIndexToP1 == 0 || nearestCornersIndexToP1 == 2) {
            d1 = findEuclidean2DDistance(secondPoint, obstacles.create3DObstacle(obstacle).get(1));
            d2 = findEuclidean2DDistance(secondPoint, obstacles.create3DObstacle(obstacle).get(3));
            if (d1 <= d2) {
                newCornerForP2 = obstacles.create3DObstacle(obstacle).get(1);
            } else {
                newCornerForP2 = obstacles.create3DObstacle(obstacle).get(3);
            }
        } else if (nearestCornersIndexToP1 == 1 || nearestCornersIndexToP1 == 3) {
            d1 = findEuclidean2DDistance(secondPoint, obstacles.create3DObstacle(obstacle).get(0));
            d2 = findEuclidean2DDistance(secondPoint, obstacles.create3DObstacle(obstacle).get(2));
            if (d1 <= d2) {
                newCornerForP2 = obstacles.create3DObstacle(obstacle).get(0);
            } else {
                newCornerForP2 = obstacles.create3DObstacle(obstacle).get(2);
            }
        }

        return newCornerForP2;
    }
}