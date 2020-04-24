package com.thealiyev;

import java.util.ArrayList;

public class ObstacleAvoider {

    public static void main(String[] args) {
        ArrayList<Double> ObstacleAvoidanceCurrentStation = new ArrayList<>();
        ArrayList<Double> ObstacleAvoidanceNextStation = new ArrayList<>();
        Obstacles obstacles = new Obstacles();
        obstacles.setObstacle1();
        obstacles.setObstacle2();
        boolean isPointInsideOfObstacle = false;
        boolean didPointCollideWithObstacle = false;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;

        ObstacleAvoidanceCurrentStation.add(40.0);
        ObstacleAvoidanceCurrentStation.add(38.0);

        ObstacleAvoidanceNextStation.add(60.0);
        ObstacleAvoidanceNextStation.add(55.0);


        isPointInsideOfObstacle = isPointInsideOfObstacle(ObstacleAvoidanceNextStation, obstacles, 1);
        if (isPointInsideOfObstacle) {
            System.out.println("Point inside do something!");
            ObstacleAvoidanceNextStation = findNearestPoint(ObstacleAvoidanceNextStation, obstacles, 1);
            System.out.println(ObstacleAvoidanceNextStation);
        } else {
            didPointCollideWithObstacle = didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, obstacles, 1);
            if (didPointCollideWithObstacle) {
                System.out.println("Point collided do something!");
                nearestCornerToCurrentStation = findNearestPoint(ObstacleAvoidanceCurrentStation, obstacles, 1);
                newCornerForNextStation = findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, obstacles, 1);
                System.out.println(ObstacleAvoidanceCurrentStation + " " + nearestCornerToCurrentStation + " to " + newCornerForNextStation + " " + ObstacleAvoidanceNextStation);
            }
        }
    }

    public static boolean isPointInsideOfObstacle(ArrayList<Double> p, Obstacles obstacles, int turn) {
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

    public static boolean didPointCollideWithObstacle(ArrayList<Double> p1, ArrayList<Double> p2, Obstacles obstacles, int turn) {
        boolean didCollide = false;
        ArrayList<Double> O = new ArrayList<>();
        double A, h, R;
        double cosAlpha, alphaDegree;

        O.add(obstacles.findCenterOfObstacle(obstacles.create3DObstacle(obstacles.getObstacles().get(turn))).get(0));
        O.add(obstacles.findCenterOfObstacle(obstacles.create3DObstacle(obstacles.getObstacles().get(turn))).get(1));

        A = calculateTriangleArea(p1, p2, O);
        h = (2 * A) / findEuclidean2DDistance(p2, p1);
        R = obstacles.getR(obstacles.create3DObstacle(obstacles.getObstacles().get(turn)));
        cosAlpha = (Math.pow(findEuclidean2DDistance(O, p1), 2) + Math.pow(findEuclidean2DDistance(p1, p2), 2) - Math.pow(findEuclidean2DDistance(O, p2), 2)) / (2 * findEuclidean2DDistance(O, p1) * findEuclidean2DDistance(p1, p2));
        alphaDegree = Math.toDegrees(Math.acos(cosAlpha));

        double Xmin = obstacles.create3DObstacle(obstacles.getObstacles().get(1)).get(0).get(0), Ymin = obstacles.create3DObstacle(obstacles.getObstacles().get(1)).get(0).get(1);
        double x, y;

        x = p1.get(0);
        y = p1.get(1);

        if (x >= Xmin && y >= Ymin) {
            x = p2.get(0);
            y = p2.get(1);
            if (x >= Xmin && y >= Ymin) {
                didCollide = false;
            }
        } else {
            if (h <= R) {
                didCollide = true;
            }
        }


        /*if (alphaDegree > 90) {
            didCollide = false;
        } else {
            if (h <= R) {
                didCollide = true;
            }
        }*/

        return didCollide;
    }

    public static double calculateTriangleArea(ArrayList<Double> p1, ArrayList<Double> p2, ArrayList<Double> obstacle) {
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

    public static double findEuclidean2DDistance(ArrayList<Double> firstPoint, ArrayList<Double> secondPoint) {
        double euclideanDistance = 0.0;

        for (int stCounter = 0; stCounter < 2; stCounter = stCounter + 1) {
            euclideanDistance = euclideanDistance + Math.pow((secondPoint.get(stCounter) - firstPoint.get(stCounter)), 2);
        }

        euclideanDistance = Math.sqrt(euclideanDistance);

        return euclideanDistance;
    }

    public static ArrayList<Double> findNearestPoint(ArrayList<Double> p, Obstacles obstacles, int turn) {
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

    public static ArrayList<Double> findPathToOppositeSide(ArrayList<Double> firstPoint, ArrayList<Double> secondPoint, Obstacles obstacles, int turn) {
        ArrayList<Double> nearestCornerToP1 = findNearestPoint(firstPoint, obstacles, turn);
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