package com.thealiyev;

import java.util.ArrayList;
import java.util.Random;

public class IGWOPathPlanningDynamicStationaryMethod {
    private static Random random = null;

    public static void main(String[] args) {
        IGWOPathPlanningDynamicStationaryMethod igwoPathPlanningDynamicStationaryMethod = new IGWOPathPlanningDynamicStationaryMethod();
        igwoPathPlanningDynamicStationaryMethod.IGWO();
    }

    private void IGWO() {
        random = new Random();
        //Gray Wolf Optimization and Path Planning start here...
        //Boundaries of map
        ArrayList<Double> Xboundaries = new ArrayList<>(), Yboundaries = new ArrayList<>(), Zboundaries = new ArrayList<>();
        //X boundaries
        Xboundaries.add(0.0);
        Xboundaries.add(50.0);
        //Y boundaries
        Yboundaries.add(0.0);
        Yboundaries.add(50.0);
        //Z boundaries
        Zboundaries.add(0.0);
        Zboundaries.add(50.0);
        //Source station coordinates
        ArrayList<Double> sourceStation = new ArrayList<>();
        sourceStation.add(1.0);
        sourceStation.add(1.0);
        sourceStation.add(1.0);
        //Destination station coordinates
        ArrayList<Double> destinationStation = new ArrayList<>();
        destinationStation.add(49.0);
        destinationStation.add(49.0);
        destinationStation.add(49.0);
        //Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        ArrayList<Double> A;
        ArrayList<Double> C;
        ArrayList<Double> D;
        ArrayList<Double> X;
        double x;
        ArrayList<Double> P;
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
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        double hMin, hMax, z;
        //Gray Wolf Optimization iterations start here...
        System.out.println("Initialization, alpha's fitness value: " + sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * stCounter / iteration;
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    A = new ArrayList<>();
                    C = new ArrayList<>();
                    D = new ArrayList<>();
                    X = new ArrayList<>();
                    x = optimizationMatrix.get(ndCounter).get(rdCounter);

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A.add(2 * a * r1 - a);
                    C.add(2 * r2);
                    D.add(C.get(0) * optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(rdCounter) - x);
                    if (D.get(0) < 0) {
                        D.set(0, -1 * D.get(0));
                    }
                    X.add(optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(rdCounter) - A.get(0) * D.get(0));

                    if (ndCounter > 0) {
                        for (int fourthCounter = 0; fourthCounter < ndCounter; fourthCounter = fourthCounter + 1) {
                            r1 = random.nextDouble();
                            r2 = random.nextDouble();
                            A.add(2 * a * r1 - a);
                            C.add(2 * r2);
                            D.add(C.get(fourthCounter) * optimizationMatrix.get(ndCounter - 1).get(rdCounter) - x);
                            if (D.get(fourthCounter) < 0) {
                                D.set(fourthCounter, -1 * D.get(fourthCounter));
                            }
                            X.add(optimizationMatrix.get(ndCounter - 1).get(rdCounter) - A.get(fourthCounter) * D.get(fourthCounter));
                        }
                    }

                    x = 0;
                    for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                        x = x + X.get(fifthCounter);
                    }
                    x = x / X.size();
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
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(2, positionsMatrixWithoutCollisions.get(0).get(rdCounter - 1).get(2));
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        z = ObstacleAvoidanceCurrentStation.get(2);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            hMin = sortedObstacles.get(fourthCounter).get(0).get(2);
                            hMax = sortedObstacles.get(fourthCounter).get(1).get(2);
                            if (z >= hMin && z <= hMax) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                if (isPointInsideOfObstacle) {
                                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                    positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                                    positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                                }
                            }
                        }
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
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(2, positionsMatrixWithoutCollisions.get(0).get(0).get(2));
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        z = ObstacleAvoidanceCurrentStation.get(2);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            hMin = sortedObstacles.get(fourthCounter).get(0).get(2);
                            hMax = sortedObstacles.get(fourthCounter).get(1).get(2);
                            if (z >= hMin && z <= hMax) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                if (isPointInsideOfObstacle) {
                                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                    positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                                    positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                                }
                            }
                        }
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
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        z = ObstacleAvoidanceCurrentStation.get(2);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            hMin = sortedObstacles.get(fourthCounter).get(0).get(2);
                            hMax = sortedObstacles.get(fourthCounter).get(1).get(2);
                            if (z >= hMin && z <= hMax) {
                                didPointCollideWithObstacle = obstacleAvoider.didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                                if (didPointCollideWithObstacle) {
                                    nearestCornerToCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                    nearestCornerToCurrentStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(2));
                                    pathWithCollisiions.add(nearestCornerToCurrentStation);
                                    newCornerForNextStation = obstacleAvoider.findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                                    newCornerForNextStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(2));
                                    pathWithCollisiions.add(newCornerForNextStation);
                                    ObstacleAvoidanceCurrentStation = newCornerForNextStation;
                                }
                            }
                        }
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            z = ObstacleAvoidanceCurrentStation.get(2);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                                hMin = sortedObstacles.get(fourthCounter).get(0).get(2);
                                hMax = sortedObstacles.get(fourthCounter).get(1).get(2);
                                if (z >= hMin && z <= hMax) {
                                    didPointCollideWithObstacle = obstacleAvoider.didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                                    if (didPointCollideWithObstacle) {
                                        nearestCornerToCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                        nearestCornerToCurrentStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(2));
                                        pathWithCollisiions.add(nearestCornerToCurrentStation);
                                        newCornerForNextStation = obstacleAvoider.findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                                        newCornerForNextStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(2));
                                        pathWithCollisiions.add(newCornerForNextStation);
                                        ObstacleAvoidanceCurrentStation = newCornerForNextStation;
                                    }
                                }
                            }
                        }
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        z = ObstacleAvoidanceCurrentStation.get(2);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            hMin = sortedObstacles.get(fourthCounter).get(0).get(2);
                            hMax = sortedObstacles.get(fourthCounter).get(1).get(2);
                            if (z >= hMin && z <= hMax) {
                                didPointCollideWithObstacle = obstacleAvoider.didPointCollideWithObstacle(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                                if (didPointCollideWithObstacle) {
                                    nearestCornerToCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                    nearestCornerToCurrentStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(0).get(2));
                                    pathWithCollisiions.add(nearestCornerToCurrentStation);
                                    newCornerForNextStation = obstacleAvoider.findPathToOppositeSide(ObstacleAvoidanceCurrentStation, ObstacleAvoidanceNextStation, sortedObstacles.get(fourthCounter));
                                    newCornerForNextStation.add(positionsMatrixWithoutCollisions.get(ndCounter).get(0).get(2));
                                    pathWithCollisiions.add(newCornerForNextStation);
                                    ObstacleAvoidanceCurrentStation = newCornerForNextStation;
                                }
                            }
                        }
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
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
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
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

        for (int stCounter = 0; stCounter < population; stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < dimension; ndCounter = ndCounter + 1) {
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

    private ArrayList<ArrayList<Double>> createOptimizationMatrix(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrix, ArrayList<Double> sourceStation) {
        random = new Random();

        ArrayList<ArrayList<Double>> optimizationMatrix = new ArrayList<>();
        ArrayList<Double> optimizationVector = new ArrayList<>();
        double euclideanDistance;

        for (int stCounter = 0; stCounter < positionsMatrix.size(); stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < positionsMatrix.get(stCounter).size(); ndCounter = ndCounter + 1) {
                if (ndCounter == 0) {
                    euclideanDistance = findEuclideanDistance(sourceStation, positionsMatrix.get(stCounter).get(ndCounter));
                    optimizationVector.add(euclideanDistance);
                } else {
                    euclideanDistance = findEuclideanDistance(positionsMatrix.get(stCounter).get(ndCounter - 1), positionsMatrix.get(stCounter).get(ndCounter));
                    optimizationVector.add(euclideanDistance);
                }
            }
            optimizationMatrix.add(optimizationVector);
            optimizationVector = new ArrayList<>();
        }

        return optimizationMatrix;
    }

    private ArrayList<Double> calculateFitnessValues(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions, ArrayList<Double> sourceStation, ArrayList<Double> destinatioStation) {
        ArrayList<Double> fitnessValues = new ArrayList<>();
        double sum = 0.0;
        double euclideanDistance;

        for (int stCounter = 0; stCounter < positionsMatrixWithCollisions.size(); stCounter = stCounter + 1) {
            euclideanDistance = findEuclideanDistance(sourceStation, positionsMatrixWithCollisions.get(stCounter).get(0));
            sum = sum + euclideanDistance;

            euclideanDistance = findEuclideanDistance(positionsMatrixWithCollisions.get(stCounter).get(positionsMatrixWithCollisions.get(stCounter).size() - 1), destinatioStation);
            sum = sum + euclideanDistance;

            for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(stCounter).size() - 1; ndCounter = ndCounter + 1) {
                euclideanDistance = findEuclideanDistance(positionsMatrixWithCollisions.get(stCounter).get(ndCounter), positionsMatrixWithCollisions.get(stCounter).get(ndCounter + 1));
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