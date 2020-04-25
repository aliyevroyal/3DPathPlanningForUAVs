package com.thealiyev;

import sun.jvm.hotspot.runtime.StubRoutines;

import java.util.ArrayList;
import java.util.Random;

public class GWOPathPlanningDynamicStationaryMethod {
    private static Random random = null;

    public static void main(String[] args) {
        GWOPathPlanningDynamicStationaryMethod gwoPathPlanningDynamicStationaryMethod = new GWOPathPlanningDynamicStationaryMethod();
        gwoPathPlanningDynamicStationaryMethod.GWO();
    }

    private void GWO() {
        random = new Random();
        //Gray Wolf Optimization and Path Planning start here...
        //Boundaries of map
        ArrayList<Double> Xboundaries = new ArrayList<>(), Yboundaries = new ArrayList<>(), Zboundaries = new ArrayList<>();
        //X boundaries
        Xboundaries.add(0.0);
        Xboundaries.add(100.0);
        //Y boundaries
        Yboundaries.add(0.0);
        Yboundaries.add(100.0);
        //Z boundaries
        Zboundaries.add(0.0);
        Zboundaries.add(100.0);
        //Source station coordinates
        ArrayList<Double> sourceStation = new ArrayList<>();
        sourceStation.add(0.1);
        sourceStation.add(0.1);
        sourceStation.add(0.1);
        //Destination station coordinates
        ArrayList<Double> destinationStation = new ArrayList<>();
        destinationStation.add(99.9);
        destinationStation.add(99.9);
        destinationStation.add(99.9);

        //Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        double A1, A2, A3;
        double C1, C2, C3;
        double x, X1, X2, X3;
        ArrayList<Double> P;
        double Dalpha, Dbeta, Ddelta;
        double Xalpha, Xbeta, Xdelta;
        int population = 100, dimension = 5;
        int iteration = 100;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions = createRandomPositionsMatrix(population, dimension, Xboundaries, Yboundaries, Zboundaries);
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = calculateFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        //Path Planning components initialization starts here...
        double r = (findEuclideanDistance(sourceStation, destinationStation)) / (dimension + 1);
        double xCurrent, yCurrent;
        double xDestination = destinationStation.get(0);
        double yDestination = destinationStation.get(1);
        double xNext, yNext;
        double d, cosAlpha, sinAlpha, alphaDegree;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        Obstacles obstacles = new Obstacles();
        obstacles.setObstacle1();
        obstacles.setObstacle2();
        obstacles.setObstacle3();
        obstacles.setObstacle4();
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        //Gray Wolf Optimization iterations start here...
        System.out.println("Initialization, alpha's fitness value: " + sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * stCounter / iteration;
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    x = optimizationMatrix.get(ndCounter).get(rdCounter);

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A1 = 2 * a * r1 - a;
                    C1 = 2 * r2;
                    Xalpha = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(rdCounter);
                    Dalpha = C1 * Xalpha - x;
                    if (Dalpha < 0) {
                        Dalpha = Dalpha * -1;
                    }
                    X1 = Xalpha - A1 * Dalpha;

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A2 = 2 * a * r1 - a;
                    C2 = 2 * r2;
                    Xbeta = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(1))).get(rdCounter);
                    Dbeta = C2 * Xbeta - x;
                    if (Dbeta < 0) {
                        Dbeta = Dbeta * -1;
                    }
                    X2 = Xbeta - A2 * Dbeta;

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A3 = 2 * a * r1 - a;
                    C3 = 2 * r2;
                    Xdelta = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(2))).get(rdCounter);
                    Ddelta = C3 * Xdelta - x;
                    if (Ddelta < 0) {
                        Ddelta = Ddelta * -1;
                    }
                    X3 = Xdelta - A3 * Ddelta;

                    x = (X1 + X2 + X3) / 3;
                    if (x > r) {
                        x = r + (x - r) * random.nextDouble();
                    } else {
                        x = x + (r - x) * random.nextDouble();
                    }

                    //Update Position
                    if (rdCounter > 0) {
                        xCurrent = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(0);
                        yCurrent = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(1);
                        d = Math.sqrt(Math.pow((xDestination - xCurrent), 2) + Math.pow((yDestination - yCurrent), 2));
                        cosAlpha = (xDestination - xCurrent) / d;
                        sinAlpha = (yDestination - yCurrent) / d;
                        alphaDegree = Math.toDegrees(Math.acos(cosAlpha));
                        alphaDegree = (alphaDegree - 15 * a) + ((alphaDegree + 15 * a) - (alphaDegree - 15 * a)) * random.nextDouble();
                        xNext = xCurrent + (x * Math.cos(Math.toRadians(alphaDegree)));
                        alphaDegree = Math.toDegrees(Math.asin(sinAlpha));
                        alphaDegree = (alphaDegree - 15 * a) + ((alphaDegree + 15 * a) - (alphaDegree - 15 * a)) * random.nextDouble();
                        yNext = yCurrent + (x * Math.sin(Math.toRadians(alphaDegree)));
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, xNext);
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, yNext);
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        for (int fourthCounter = 0; fourthCounter < obstacles.getObstacles().size(); fourthCounter = fourthCounter + 1) {
                            if (ObstacleAvoidanceNextStation.get(2) >= obstacles.getObstacles().get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= obstacles.getObstacles().get(fourthCounter).get(1).get(2)) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceNextStation, obstacles, fourthCounter);
                                if (isPointInsideOfObstacle) {
                                    if (fourthCounter == 0) {
                                        ObstacleAvoidanceNextStation = obstacleAvoider.findNearestCorner(destinationStation, obstacles, fourthCounter);
                                    } else {
                                        ObstacleAvoidanceNextStation = obstacleAvoider.findNearestCorner(obstacles.getObstacles().get(fourthCounter - 1).get(1), obstacles, fourthCounter);
                                    }
                                    ObstacleAvoidanceNextStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).get(2));
                                    positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceNextStation.get(0));
                                    positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceNextStation.get(1));
                                } else {
                                    didPointCollideWithObstacle = obstacleAvoider.didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, obstacles, fourthCounter);
                                    if (didPointCollideWithObstacle) {
                                        nearestCornerToCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, obstacles, fourthCounter);
                                        nearestCornerToCurrentStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).get(2));
                                        newCornerForNextStation = obstacleAvoider.findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, obstacles, fourthCounter);
                                        newCornerForNextStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).get(2));
                                        pathWithCollisiions.add(nearestCornerToCurrentStation);
                                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, newCornerForNextStation.get(0));
                                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, newCornerForNextStation.get(1));
                                    }
                                }
                            }
                        }
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(2, positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(2));
                    } else if (rdCounter == 0) {
                        xCurrent = sourceStation.get(0);
                        yCurrent = sourceStation.get(1);
                        d = Math.sqrt(Math.pow((xDestination - xCurrent), 2) + Math.pow((yDestination - yCurrent), 2));
                        cosAlpha = (xDestination - xCurrent) / d;
                        sinAlpha = (yDestination - yCurrent) / d;
                        alphaDegree = Math.toDegrees(Math.acos(cosAlpha));
                        xNext = xCurrent + (x * Math.cos(Math.toRadians(alphaDegree)));
                        alphaDegree = Math.toDegrees(Math.asin(sinAlpha));
                        yNext = yCurrent + (x * Math.sin(Math.toRadians(alphaDegree)));
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, xNext);
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, yNext);
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        for (int fourthCounter = 0; fourthCounter < obstacles.getObstacles().size(); fourthCounter = fourthCounter + 1) {
                            if (ObstacleAvoidanceNextStation.get(2) >= obstacles.getObstacles().get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= obstacles.getObstacles().get(fourthCounter).get(1).get(2)) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceNextStation, obstacles, fourthCounter);
                                if (isPointInsideOfObstacle) {
                                    ObstacleAvoidanceNextStation = obstacleAvoider.findNearestCorner(destinationStation, obstacles, fourthCounter);
                                    ObstacleAvoidanceNextStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).get(2));
                                    positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceNextStation.get(0));
                                    positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceNextStation.get(1));
                                } else {
                                    didPointCollideWithObstacle = obstacleAvoider.didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, obstacles, fourthCounter);
                                    if (didPointCollideWithObstacle) {
                                        nearestCornerToCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, obstacles, fourthCounter);
                                        nearestCornerToCurrentStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).get(2));
                                        newCornerForNextStation = obstacleAvoider.findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, obstacles, fourthCounter);
                                        newCornerForNextStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).get(2));
                                        pathWithCollisiions.add(nearestCornerToCurrentStation);
                                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, newCornerForNextStation.get(0));
                                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, newCornerForNextStation.get(1));
                                    }
                                }
                            }
                        }
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }

                    P = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    //Controls Boundaries
                    if (P.get(0) < Xboundaries.get(0)) {
                        P.set(0, Xboundaries.get(0));
                    } else if (P.get(0) > Xboundaries.get(1)) {
                        P.set(0, Xboundaries.get(1));
                    }
                    if (P.get(1) < Yboundaries.get(0)) {
                        P.set(1, Yboundaries.get(0));
                    } else if (P.get(1) > Yboundaries.get(1)) {
                        P.set(1, Yboundaries.get(1));
                    }
                    if (P.get(2) < Zboundaries.get(0)) {
                        P.set(2, Zboundaries.get(0));
                    } else if (P.get(2) > Zboundaries.get(1)) {
                        P.set(2, Zboundaries.get(1));
                    }
                    //Updates Boundaries
                    positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, P);
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int counter = 0; counter < pathWithCollisiions.size(); counter = counter + 1) {
                    if (counter > 0) {
                        if (pathWithCollisiions.get(counter).get(0).equals(pathWithCollisiions.get(counter - 1).get(0))) {
                            if (pathWithCollisiions.get(counter).get(1).equals(pathWithCollisiions.get(counter - 1).get(1))) {
                                if (pathWithCollisiions.get(counter).get(2).equals(pathWithCollisiions.get(counter - 1).get(2))) {
                                    pathWithCollisiions.remove(new ArrayList<>(pathWithCollisiions.get(counter)));
                                }
                            }
                        }
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
                optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
            }
            fitnessValues = calculateFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(stCounter + " iteration, alpha's fitness value: " + sortedFitnessValues.get(0));
        }

        System.out.println("Alpha path without collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
    }

    private ArrayList<ArrayList<ArrayList<Double>>> createRandomPositionsMatrix(int population, int dimension,
                                                                                ArrayList<Double> Xboundaries,
                                                                                ArrayList<Double> Yboundaries,
                                                                                ArrayList<Double> Zboundaries) {
        random = new Random();

        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrix = new ArrayList<>();
        ArrayList<ArrayList<Double>> positionsVector = new ArrayList<>();
        ArrayList<Double> position = new ArrayList<>();
        double X, Y, Z;
        //double Xmin, Ymin, Zmin;
        //double Xmax, Ymax, Zmax;

        for (int stCounter = 0; stCounter < population; stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < dimension; ndCounter = ndCounter + 1) {
                //More deterministic way
                /*Xmin = ((Xboundaries.get(1) - Xboundaries.get(0)) / dimension) * ndCounter;
                Xmax = ((Xboundaries.get(1) - Xboundaries.get(0)) / dimension) * (ndCounter + 1);
                Ymin = ((Yboundaries.get(1) - Yboundaries.get(0)) / dimension) * ndCounter;
                Ymax = ((Yboundaries.get(1) - Yboundaries.get(0)) / dimension) * (ndCounter + 1);
                Zmin = ((Zboundaries.get(1) - Zboundaries.get(0)) / dimension) * ndCounter;
                Zmax = ((Zboundaries.get(1) - Zboundaries.get(0)) / dimension) * (ndCounter + 1);
                X = Xmin + (Xmax - Xmin) * random.nextDouble();
                Y = Ymin + (Ymax - Ymin) * random.nextDouble();
                Z = Zmin + (Zmax - Zmin) * random.nextDouble();*/
                X = Xboundaries.get(0) + (Xboundaries.get(1) - Xboundaries.get(0)) * random.nextDouble();
                Y = Yboundaries.get(0) + (Yboundaries.get(1) - Yboundaries.get(0)) * random.nextDouble();
                Z = Zboundaries.get(0) + (Zboundaries.get(1) - Zboundaries.get(0)) * random.nextDouble();

                position.add(X);
                position.add(Y);
                position.add(Z);

                positionsVector.add(position);
                position = new ArrayList<>();
            }
            positionsMatrix.add(positionsVector);
            positionsVector = new ArrayList<>();
        }

        return positionsMatrix;
    }

    private ArrayList<ArrayList<Double>> createOptimizationMatrix(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions, ArrayList<Double> sourceStation) {
        random = new Random();

        ArrayList<ArrayList<Double>> optimizationMatrix = new ArrayList<>();
        ArrayList<Double> optimizationVector = new ArrayList<>();
        double euclideanDistance;

        for (int stCounter = 0; stCounter < positionsMatrixWithoutCollisions.size(); stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(stCounter).size(); ndCounter = ndCounter + 1) {
                if (ndCounter == 0) {
                    euclideanDistance = findEuclideanDistance(sourceStation, positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter));
                    optimizationVector.add(euclideanDistance);
                } else {
                    euclideanDistance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter - 1), positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter));
                    optimizationVector.add(euclideanDistance);
                }
            }
            optimizationMatrix.add(optimizationVector);
            optimizationVector = new ArrayList<>();
        }

        return optimizationMatrix;
    }

    private ArrayList<Double> calculateFitnessValues(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrix, ArrayList<Double> sourceStation, ArrayList<Double> destinatioStation) {
        ArrayList<Double> fitnessValues = new ArrayList<>();
        double sum = 0.0;
        double euclideanDistance;

        for (int stCounter = 0; stCounter < positionsMatrix.size(); stCounter = stCounter + 1) {
            euclideanDistance = findEuclideanDistance(sourceStation, positionsMatrix.get(stCounter).get(0));
            sum = sum + euclideanDistance;

            euclideanDistance = findEuclideanDistance(positionsMatrix.get(stCounter).get(positionsMatrix.get(stCounter).size() - 1), destinatioStation);
            sum = sum + euclideanDistance;

            for (int ndCounter = 0; ndCounter < positionsMatrix.get(stCounter).size() - 1; ndCounter = ndCounter + 1) {
                euclideanDistance = findEuclideanDistance(positionsMatrix.get(stCounter).get(ndCounter), positionsMatrix.get(stCounter).get(ndCounter + 1));
                sum = sum + euclideanDistance;
            }
            fitnessValues.add(sum);
        }

        return fitnessValues;
    }

    private double findEuclideanDistance(ArrayList<Double> firstPoint, ArrayList<Double> secondPoint) {
        double euclideanDistance = 0.0;

        for (int stCounter = 0; stCounter < 3; stCounter = stCounter + 1) {
            euclideanDistance = euclideanDistance + Math.pow((secondPoint.get(stCounter) - firstPoint.get(stCounter)), 2);
        }

        euclideanDistance = Math.sqrt(euclideanDistance);

        return euclideanDistance;
    }

    private ArrayList<Double> sortFitnessValues(ArrayList<Double> fitnessValues) {
        ArrayList<Double> duplicatedFitnessValues = new ArrayList<>(fitnessValues);
        ArrayList<Double> sortedFitnessValues = new ArrayList<>();

        double min;
        for (int stCounter = 0; stCounter < fitnessValues.size(); stCounter = stCounter + 1) {
            min = duplicatedFitnessValues.get(0);
            for (int ndCounter = 0; ndCounter < duplicatedFitnessValues.size(); ndCounter = ndCounter + 1) {
                if (duplicatedFitnessValues.get(ndCounter) < min) {
                    min = duplicatedFitnessValues.get(ndCounter);
                }
            }
            sortedFitnessValues.add(min);
            duplicatedFitnessValues.remove(new Double(min));
        }

        return sortedFitnessValues;
    }
}