package com.thealiyev;

import java.util.ArrayList;
import java.util.Random;

public class AdaptedRRTExperimentResults {
    private Random random = null;

    public static void main(String[] args) {
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
        //Initialization starts here...
        int population = 100, dimension = 5;
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
        AdaptedRRTExperimentResults adaptedRRTExperimentResults = new AdaptedRRTExperimentResults();
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions = adaptedRRTExperimentResults.createRandomPositionsMatrix(population, dimension, Xboundaries, Yboundaries, Zboundaries);
        ArrayList<ArrayList<ArrayList<Double>>> GWOPositionsMatrixWithoutCollisions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> IGWOPositionsMatrixWithoutCollisions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> ExGWOPositionsMatrixWithoutCollisions = new ArrayList<>();
        ArrayList<Double> vector = new ArrayList();
        ArrayList<ArrayList<Double>> matrix = new ArrayList<>();

        for (int stCounter = 0; stCounter < positionsMatrixWithoutCollisions.size(); stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(stCounter).size(); ndCounter = ndCounter + 1) {
                for (int rdCounter = 0; rdCounter < positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    vector.add(positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter).get(rdCounter));
                }
                matrix.add(vector);
                vector = new ArrayList<>();
            }
            GWOPositionsMatrixWithoutCollisions.add(matrix);
            matrix = new ArrayList<>();
        }

        for (int stCounter = 0; stCounter < positionsMatrixWithoutCollisions.size(); stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(stCounter).size(); ndCounter = ndCounter + 1) {
                for (int rdCounter = 0; rdCounter < positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    vector.add(positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter).get(rdCounter));
                }
                matrix.add(vector);
                vector = new ArrayList<>();
            }
            IGWOPositionsMatrixWithoutCollisions.add(matrix);
            matrix = new ArrayList<>();
        }

        for (int stCounter = 0; stCounter < positionsMatrixWithoutCollisions.size(); stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(stCounter).size(); ndCounter = ndCounter + 1) {
                for (int rdCounter = 0; rdCounter < positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    vector.add(positionsMatrixWithoutCollisions.get(stCounter).get(ndCounter).get(rdCounter));
                }
                matrix.add(vector);
                vector = new ArrayList<>();
            }
            ExGWOPositionsMatrixWithoutCollisions.add(matrix);
            matrix = new ArrayList<>();
        }

        AdaptedRRTGWOExperimentResults gwoExperimentResults = new AdaptedRRTGWOExperimentResults();
        gwoExperimentResults.GWO(GWOPositionsMatrixWithoutCollisions);

        AdaptedRRTIGWOExperimentResults igwoExperimentResults = new AdaptedRRTIGWOExperimentResults();
        igwoExperimentResults.IGWO(IGWOPositionsMatrixWithoutCollisions);

        AdaptedRRTExGWOExperimentResults exGWOExperimentResults = new AdaptedRRTExGWOExperimentResults();
        exGWOExperimentResults.ExGWO(ExGWOPositionsMatrixWithoutCollisions);
    }

    public ArrayList<ArrayList<ArrayList<Double>>> createRandomPositionsMatrix(int population, int dimension,
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

    public ArrayList<ArrayList<Double>> createOptimizationMatrix
            (ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions, ArrayList<Double> sourceStation) {
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

    public ArrayList<Double> calculateFitnessValues
            (ArrayList<ArrayList<ArrayList<Double>>> positionsMatrix, ArrayList<Double> sourceStation, ArrayList<Double> destinatioStation) {
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

    public double findEuclideanDistance(ArrayList<Double> firstPoint, ArrayList<Double> secondPoint) {
        double euclideanDistance = 0.0;

        for (int stCounter = 0; stCounter < 3; stCounter = stCounter + 1) {
            euclideanDistance = euclideanDistance + Math.pow((secondPoint.get(stCounter) - firstPoint.get(stCounter)), 2);
        }

        euclideanDistance = Math.sqrt(euclideanDistance);

        return euclideanDistance;
    }

    public ArrayList<Double> sortFitnessValues(ArrayList<Double> fitnessValues) {
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


class AdaptedRRTGWOExperimentResults {
    private Random random = null;
    private AdaptedRRTExperimentResults adaptedRRTExperimentResults = null;

    public void GWO(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        adaptedRRTExperimentResults = new AdaptedRRTExperimentResults();
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
        double A1, A2, A3;
        double C1, C2, C3;
        double x, X1, X2, X3;
        ArrayList<Double> P;
        double Dalpha, Dbeta, Ddelta;
        double Xalpha, Xbeta, Xdelta;
        int population = 100, dimension = 5;
        int iteration = 100;
        ArrayList<ArrayList<Double>> optimizationMatrix = adaptedRRTExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = adaptedRRTExperimentResults.calculateFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = adaptedRRTExperimentResults.sortFitnessValues(fitnessValues);
        //Path Planning components initialization starts here...
        double r = (adaptedRRTExperimentResults.findEuclideanDistance(sourceStation, destinationStation)) / (dimension + 1);
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
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(2, positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(2));
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                if (isPointInsideOfObstacle) {
                                    pathWithCollisiions.remove(pathWithCollisiions.size() - 1);
                                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                } else {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                                if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                    isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                    if (isPointInsideOfObstacle) {
                                        ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                    } else {
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
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                if (isPointInsideOfObstacle) {
                                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                } else {
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
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
                optimizationMatrix = adaptedRRTExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
            }
            fitnessValues = adaptedRRTExperimentResults.calculateFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = adaptedRRTExperimentResults.sortFitnessValues(fitnessValues);
        }

        System.out.println("GWO Experiment Results: ");
        System.out.println("Alpha's fitness value: " + sortedFitnessValues.get(0));
        System.out.println("Alpha path without collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
    }
}

class AdaptedRRTIGWOExperimentResults {
    private Random random = null;
    private AdaptedRRTExperimentResults adaptedRRTExperimentResults = null;

    public void IGWO(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        adaptedRRTExperimentResults = new AdaptedRRTExperimentResults();
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
        ArrayList<ArrayList<Double>> optimizationMatrix = adaptedRRTExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = adaptedRRTExperimentResults.calculateFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = adaptedRRTExperimentResults.sortFitnessValues(fitnessValues);
        //Path Planning components initialization starts here...
        double r = (adaptedRRTExperimentResults.findEuclideanDistance(sourceStation, destinationStation)) / (dimension + 1);
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
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(2, positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(2));
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                if (isPointInsideOfObstacle) {
                                    pathWithCollisiions.remove(pathWithCollisiions.size() - 1);
                                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                } else {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                                if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                    isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                    if (isPointInsideOfObstacle) {
                                        ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                    } else {
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
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                if (isPointInsideOfObstacle) {
                                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                } else {
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
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
                optimizationMatrix = adaptedRRTExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
            }
            fitnessValues = adaptedRRTExperimentResults.calculateFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = adaptedRRTExperimentResults.sortFitnessValues(fitnessValues);
        }

        System.out.println("IGWO Experiment Results: ");
        System.out.println("Alpha's fitness value: " + sortedFitnessValues.get(0));
        System.out.println("Alpha path without collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
    }
}

class AdaptedRRTExGWOExperimentResults {
    private Random random = null;
    private AdaptedRRTExperimentResults adaptedRRTExperimentResults = null;

    public void ExGWO(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        adaptedRRTExperimentResults = new AdaptedRRTExperimentResults();

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
        ArrayList<ArrayList<Double>> optimizationMatrix = adaptedRRTExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = adaptedRRTExperimentResults.calculateFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = adaptedRRTExperimentResults.sortFitnessValues(fitnessValues);
        //Path Planning components initialization starts here...
        double r = (adaptedRRTExperimentResults.findEuclideanDistance(sourceStation, destinationStation)) / (dimension + 1);
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
        //Gray Wolf Optimization iterations start here...
        System.out.println("Initialization, alpha's fitness value: " + sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * Math.pow(stCounter, 5) / Math.pow(iteration, 5);
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

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A.add(2 * a * r1 - a);
                    C.add(2 * r2);
                    D.add(C.get(1) * optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(1))).get(rdCounter) - x);
                    if (D.get(1) < 0) {
                        D.set(1, -1 * D.get(1));
                    }
                    X.add(optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(1))).get(rdCounter) - A.get(1) * D.get(1));

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A.add(2 * a * r1 - a);
                    C.add(2 * r2);
                    D.add(C.get(2) * optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(2))).get(rdCounter) - x);
                    if (D.get(2) < 0) {
                        D.set(2, -1 * D.get(2));
                    }
                    X.add(optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(2))).get(rdCounter) - A.get(2) * D.get(2));

                    if (ndCounter > 2) {
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
                        positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(2, positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1).get(2));
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                if (isPointInsideOfObstacle) {
                                    pathWithCollisiions.remove(pathWithCollisiions.size() - 1);
                                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                } else {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                                if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                    isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                    if (isPointInsideOfObstacle) {
                                        ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                    } else {
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
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                            if (ObstacleAvoidanceNextStation.get(2) >= sortedObstacles.get(fourthCounter).get(0).get(2) && ObstacleAvoidanceNextStation.get(2) <= sortedObstacles.get(fourthCounter).get(1).get(2)) {
                                isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                                if (isPointInsideOfObstacle) {
                                    ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
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
                                } else {
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
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
                optimizationMatrix = adaptedRRTExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
            }
            fitnessValues = adaptedRRTExperimentResults.calculateFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = adaptedRRTExperimentResults.sortFitnessValues(fitnessValues);
        }

        System.out.println("ExGWO Experiment Results: ");
        System.out.println("Alpha's fitness value: " + sortedFitnessValues.get(0));
        System.out.println("Alpha path without collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
    }
}