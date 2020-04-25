package com.thealiyev;

import java.util.ArrayList;

public class ObstacleAvoider {

    public static void main(String[] args) {
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation = new ArrayList<>();
        ArrayList<Double> ObstacleAvoidanceNextStation = new ArrayList<>();
        Obstacles obstacles = new Obstacles();
        obstacles.setObstacle1();
        obstacles.setObstacle2();
        obstacles.setObstacle3();
        boolean isPointInsideOfObstacle = false;
        boolean didPointCollideWithObstacle = false;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;

        ObstacleAvoidanceCurrentStation.add(5.0);
        ObstacleAvoidanceCurrentStation.add(25.0);

        ObstacleAvoidanceNextStation.add(12.0);
        ObstacleAvoidanceNextStation.add(18.0);


        ArrayList<Double> destinationStation = new ArrayList<>();
        destinationStation.add(99.9);
        destinationStation.add(99.9);
        destinationStation.add(99.9);


        for (int fourthCounter = 0; fourthCounter < obstacles.getObstacles().size(); fourthCounter = fourthCounter + 1) {
            System.out.println("Obstacle " + fourthCounter);
            isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceNextStation, obstacles, fourthCounter);
            if (isPointInsideOfObstacle) {
                System.out.println("Point inside do something!");
                ObstacleAvoidanceNextStation = obstacleAvoider.findNearestCorner(destinationStation, obstacles, fourthCounter);
                System.out.println(ObstacleAvoidanceNextStation);
            } else {
                didPointCollideWithObstacle = obstacleAvoider.didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, obstacles, fourthCounter);
                if (didPointCollideWithObstacle) {
                    System.out.println("Point collided do something!");
                    nearestCornerToCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, obstacles, fourthCounter);
                    newCornerForNextStation = obstacleAvoider.findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, obstacles, fourthCounter);
                    System.out.println(ObstacleAvoidanceCurrentStation + " " + nearestCornerToCurrentStation + " to " + newCornerForNextStation + " " + ObstacleAvoidanceNextStation);
                }
            }
        }
    }

    public boolean isPointInsideOfObstacle(ArrayList<Double> p, Obstacles obstacles, int turn) {
        boolean isPointInsideOfObstacle = false;
        double obstacleXMin, obstacleXMax, obstacleYMin, obstacleYMax;
        obstacleXMin = obstacles.create3DObstacle((obstacles.getObstacles().get(turn))).get(0).get(0);
        obstacleXMax = obstacles.create3DObstacle((obstacles.getObstacles().get(turn))).get(2).get(0);
        obstacleYMin = obstacles.create3DObstacle((obstacles.getObstacles().get(turn))).get(0).get(1);
        obstacleYMax = obstacles.create3DObstacle((obstacles.getObstacles().get(turn))).get(2).get(1);

        if (p.get(0) < obstacleXMax && p.get(0) > obstacleXMin) {
            if (p.get(1) < obstacleYMax && p.get(1) > obstacleYMin) {
                isPointInsideOfObstacle = true;
            }
        }

        return isPointInsideOfObstacle;
    }

    public boolean didPointCollideWithObstacle(ArrayList<Double> p1, ArrayList<Double> p2, Obstacles obstacles, int turn) {
        boolean didCollide = false;
        ArrayList<Double> O = new ArrayList<>();
        double A, h, R;

        O.add(obstacles.findCenterOfObstacle(obstacles.create3DObstacle(obstacles.getObstacles().get(turn))).get(0));
        O.add(obstacles.findCenterOfObstacle(obstacles.create3DObstacle(obstacles.getObstacles().get(turn))).get(1));

        A = calculateTriangleArea(p1, p2, O);
        h = (2 * A) / findEuclidean2DDistance(p2, p1);
        R = obstacles.getR(obstacles.create3DObstacle(obstacles.getObstacles().get(turn)));

        if (h < R && h >= 0) {
            if (findEuclidean2DDistance(p1, p2) > findEuclidean2DDistance(O, p2) && findEuclidean2DDistance(p1, p2) > findEuclidean2DDistance(O, p1)) {
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
        A = Math.sqrt(P * (P - a) * (P - b) * (P - c));
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

    public ArrayList<Double> findNearestCorner(ArrayList<Double> p, Obstacles obstacles, int turn) {
        ArrayList<Double> nearestPoint;
        double d = findEuclidean2DDistance(p, obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(0));
        int index = 0;

        for (int counter = 0; counter < obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).size() - 1; counter = counter + 1) {
            if (findEuclidean2DDistance(p, obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(counter)) < d) {
                d = findEuclidean2DDistance(p, obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(counter));
                index = counter;
            }
        }
        nearestPoint = obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(index);

        return nearestPoint;
    }

    public ArrayList<Double> findPathToOppositeSide(ArrayList<Double> firstPoint, ArrayList<Double> secondPoint, Obstacles obstacles, int turn) {
        ArrayList<Double> nearestCornerToP1 = findNearestCorner(firstPoint, obstacles, turn);
        int nearestCornersIndexToP1 = obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).indexOf(nearestCornerToP1);
        double d1, d2;
        ArrayList<Double> newCornerForP2 = new ArrayList<>();

        if (nearestCornersIndexToP1 == 0 || nearestCornersIndexToP1 == 2) {
            d1 = findEuclidean2DDistance(secondPoint, obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(1));
            d2 = findEuclidean2DDistance(secondPoint, obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(3));
            if (d1 <= d2) {
                newCornerForP2 = obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(1);
            } else {
                newCornerForP2 = obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(3);
            }
        } else if (nearestCornersIndexToP1 == 1 || nearestCornersIndexToP1 == 3) {
            d1 = findEuclidean2DDistance(secondPoint, obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(0));
            d2 = findEuclidean2DDistance(secondPoint, obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(2));
            if (d1 <= d2) {
                newCornerForP2 = obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(0);
            } else {
                newCornerForP2 = obstacles.create3DObstacle(obstacles.getObstacles().get(turn)).get(2);
            }
        }

        return newCornerForP2;
    }
}