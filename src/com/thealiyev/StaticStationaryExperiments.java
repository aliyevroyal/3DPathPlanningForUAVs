package com.thealiyev;

import java.util.ArrayList;
import java.util.Random;

public class StaticStationaryExperiments {
    private static StaticStationaryExperiments staticStationaryExperiments = null;
    private static Random random = null;
    //Boundaries of map
    private static ArrayList<Double> Xboundaries = null, Yboundaries = null, Zboundaries = null;
    //Source and Destination stations coordinates
    private static ArrayList<Double> sourceStation = null, destinationStation = null;
    //Obstacles
    private static Obstacles obstacles = null;
    //Optimization algorithms initialization
    private static int population, dimension;
    private static int iteration;
    private static int theNumberOfStations;

    public static void main(String[] args) {
        staticStationaryExperiments = new StaticStationaryExperiments();
        //X boundaries
        Xboundaries = new ArrayList<>();
        Xboundaries.add(0.0);
        Xboundaries.add(50.0);
        //Y boundaries
        Yboundaries = new ArrayList<>();
        Yboundaries.add(0.0);
        Yboundaries.add(50.0);
        //Z boundaries
        Zboundaries = new ArrayList<>();
        Zboundaries.add(0.0);
        Zboundaries.add(50.0);
        //Source station coordinates
        sourceStation = new ArrayList<>();
        sourceStation.add(1.0);
        sourceStation.add(1.0);
        sourceStation.add(1.0);
        //Destination station coordinates
        destinationStation = new ArrayList<>();
        destinationStation.add(49.0);
        destinationStation.add(49.0);
        destinationStation.add(49.0);
        //Obstacles
        obstacles = new Obstacles();
        obstacles.setObstacle1();
        obstacles.setObstacle2();
        obstacles.setObstacle3();
        obstacles.setObstacle4();
        //Optimization algorithms initialization
        population = 100;
        dimension = 5;
        iteration = 100;
        theNumberOfStations = 1000;
        //Variables for results of experiments
        int repeat = 10;
        ArrayList<ArrayList<Double>> stations;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions;
        ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> positionsMatricesWithoutCollisions;
        ArrayList<ArrayList<ArrayList<Double>>> matrices;
        ArrayList<ArrayList<Double>> matrix;
        ArrayList<Double> vector;
        long startTime, endTime, executionTime;
        ArrayList<Long> runtimesVector;
        ArrayList<ArrayList<Long>> runtimesMatrix = new ArrayList<>();
        ArrayList<Double> alphaValuesVector;
        ArrayList<ArrayList<Double>> alphalphaValuesMatrix = new ArrayList<>();

        for (int counter = 0; counter < repeat; counter = counter + 1) {
            System.out.println("Experimental Repeat No: " + counter);
            stations = staticStationaryExperiments.createRandomStations(theNumberOfStations, Xboundaries, Yboundaries, Zboundaries);
            positionsMatrixWithoutCollisions = staticStationaryExperiments.createRandomVisitedStations(population, dimension, stations);
            positionsMatricesWithoutCollisions = new ArrayList<>();
            matrices = new ArrayList<>();
            matrix = new ArrayList<>();
            vector = new ArrayList();
            runtimesVector = new ArrayList<>();
            alphaValuesVector = new ArrayList<>();
            for (int stCounter = 0; stCounter < 8; stCounter = stCounter + 1) {
                for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.size(); ndCounter = ndCounter + 1) {
                    for (int rdCounter = 0; rdCounter < positionsMatrixWithoutCollisions.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                        for (int fourthCounter = 0; fourthCounter < positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).size(); fourthCounter = fourthCounter + 1) {
                            vector.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).get(fourthCounter));
                        }
                        matrix.add(vector);
                        vector = new ArrayList<>();
                    }
                    matrices.add(matrix);
                    matrix = new ArrayList<>();
                }
                positionsMatricesWithoutCollisions.add(matrices);
                matrices = new ArrayList<>();
            }

            //Meta Heuristic Optimization Algorithms
            System.out.println("GWO");
            startTime = System.currentTimeMillis();
            alphaValuesVector.add(staticStationaryExperiments.GWO(stations, positionsMatricesWithoutCollisions.get(0)));
            endTime = System.currentTimeMillis();
            executionTime = (endTime - startTime);
            runtimesVector.add(executionTime);
            System.out.println("IGWO");
            startTime = System.currentTimeMillis();
            alphaValuesVector.add(staticStationaryExperiments.IGWO(stations, positionsMatricesWithoutCollisions.get(1)));
            endTime = System.currentTimeMillis();
            executionTime = (endTime - startTime);
            runtimesVector.add(executionTime);
            System.out.println("ExGWO");
            startTime = System.currentTimeMillis();
            alphaValuesVector.add(staticStationaryExperiments.ExGWO(stations, positionsMatricesWithoutCollisions.get(2)));
            endTime = System.currentTimeMillis();
            executionTime = (endTime - startTime);
            runtimesVector.add(executionTime);
            System.out.println("WOA");
            startTime = System.currentTimeMillis();
            alphaValuesVector.add(staticStationaryExperiments.WOA(stations, positionsMatricesWithoutCollisions.get(3)));
            endTime = System.currentTimeMillis();
            executionTime = (endTime - startTime);
            runtimesVector.add(executionTime);
            //Reinforcement Learning based Meta Heuristic Optimization Algorithms
            System.out.println("RLGWO");
            startTime = System.currentTimeMillis();
            alphaValuesVector.add(staticStationaryExperiments.RLGWO(stations, positionsMatricesWithoutCollisions.get(4)));
            endTime = System.currentTimeMillis();
            executionTime = (endTime - startTime);
            runtimesVector.add(executionTime);
            System.out.println("RLIGWO");
            startTime = System.currentTimeMillis();
            alphaValuesVector.add(staticStationaryExperiments.RLIGWO(stations, positionsMatricesWithoutCollisions.get(5)));
            endTime = System.currentTimeMillis();
            executionTime = (endTime - startTime);
            runtimesVector.add(executionTime);
            System.out.println("RLExGWO");
            startTime = System.currentTimeMillis();
            alphaValuesVector.add(staticStationaryExperiments.RLExGWO(stations, positionsMatricesWithoutCollisions.get(6)));
            endTime = System.currentTimeMillis();
            executionTime = (endTime - startTime);
            runtimesVector.add(executionTime);
            System.out.println("RLWOA");
            startTime = System.currentTimeMillis();
            alphaValuesVector.add(staticStationaryExperiments.RLWOA(stations, positionsMatricesWithoutCollisions.get(7)));
            endTime = System.currentTimeMillis();
            executionTime = (endTime - startTime);
            runtimesVector.add(executionTime);

            runtimesMatrix.add(runtimesVector);
            alphalphaValuesMatrix.add(alphaValuesVector);
        }
        double average;
        ArrayList<Double> runtimeAverages = new ArrayList<>();
        for (int stCounter = 0; stCounter < runtimesMatrix.get(0).size(); stCounter = stCounter + 1) {
            average = 0;
            for (int ndCounter = 0; ndCounter < runtimesMatrix.size(); ndCounter = ndCounter + 1) {
                average = average + runtimesMatrix.get(ndCounter).get(stCounter);
            }
            average = average / runtimesMatrix.size();
            runtimeAverages.add(average);
        }
        System.out.println("Averages of runtimes: " + runtimeAverages);

        double min, max;
        ArrayList<Double> mins = new ArrayList<>(), maxs = new ArrayList<>();
        ArrayList<Double> alphalphaValuesAverages = new ArrayList<>();
        for (int stCounter = 0; stCounter < alphalphaValuesMatrix.get(0).size(); stCounter = stCounter + 1) {
            average = 0;
            min = alphalphaValuesMatrix.get(0).get(stCounter);
            max = alphalphaValuesMatrix.get(0).get(stCounter);
            for (int ndCounter = 0; ndCounter < alphalphaValuesMatrix.size(); ndCounter = ndCounter + 1) {
                if (alphalphaValuesMatrix.get(ndCounter).get(stCounter) < min)
                    min = alphalphaValuesMatrix.get(ndCounter).get(stCounter);
                if (alphalphaValuesMatrix.get(ndCounter).get(stCounter) > max)
                    max = alphalphaValuesMatrix.get(ndCounter).get(stCounter);
                average = average + alphalphaValuesMatrix.get(ndCounter).get(stCounter);
            }
            mins.add(min);
            maxs.add(max);
            average = average / alphalphaValuesMatrix.size();
            alphalphaValuesAverages.add(average);
        }
        System.out.println(alphalphaValuesMatrix);
        System.out.println("Minimums: " + mins);
        System.out.println("Maximums: " + maxs);
        System.out.println("Averages of alpha values: " + alphalphaValuesAverages);
    }


    private double GWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        //Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        double A1, A2, A3;
        double C1, C2, C3;
        double x, X1, X2, X3;
        double Dalpha, Dbeta, Ddelta;
        double Xalpha, Xbeta, Xdelta;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Gray Wolf Optimization iterations start here...
        System.out.println(sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * stCounter / iteration;
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                visitingStations = new ArrayList<>();
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
                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    //Update Position
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Obstacle avoidance
                    ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                    for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                        isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        if (isPointInsideOfObstacle) {
                            ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                        }
                    }
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(sortedFitnessValues.get(0));
        }
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(0).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(0))
                    && positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(1).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }

        return sortedFitnessValues.get(0);
    }


    private double IGWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        //Iterative Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        ArrayList<Double> A;
        ArrayList<Double> C;
        ArrayList<Double> D;
        ArrayList<Double> X;
        double x;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Iterative Gray Wolf Optimization iterations start here...
        System.out.println(sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * stCounter / iteration;
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                visitingStations = new ArrayList<>();
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

                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    //Update Position
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Obstacle avoidance
                    ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                    for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                        isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        if (isPointInsideOfObstacle) {
                            ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                        }
                    }
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(sortedFitnessValues.get(0));
        }
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(0).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(0))
                    && positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(1).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }

        return sortedFitnessValues.get(0);
    }


    private double ExGWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        //Expanded Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        ArrayList<Double> A;
        ArrayList<Double> C;
        ArrayList<Double> D;
        ArrayList<Double> X;
        double x;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Expanded Gray Wolf Optimization iterations start here...
        System.out.println(sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * Math.pow(stCounter, 2) / Math.pow(iteration, 2);
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                visitingStations = new ArrayList<>();
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

                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    //Update Position
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Obstacle avoidance
                    ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                    for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                        isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        if (isPointInsideOfObstacle) {
                            ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                        }
                    }
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(sortedFitnessValues.get(0));
        }
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(0).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(0))
                    && positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(1).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }

        return sortedFitnessValues.get(0);
    }


    private double WOA(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        //Whale Optimization Algorithm initialization starts here...
        double a1, a2, r1, r2, A, C;
        boolean p;
        int randIndividual;
        double Xalpha, Dalpha, Xrand, Drand;
        double D;
        double b = 1, l;
        double x;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Whale Optimization Algorithm iterations start here...
        System.out.println(sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a1 = 2.0 - 2.0 * stCounter / iteration;
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                visitingStations = new ArrayList<>();
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    x = optimizationMatrix.get(ndCounter).get(rdCounter);

                    p = random.nextBoolean();
                    if (p) {
                        r1 = random.nextDouble();
                        A = 2 * a1 * r1 - a1;
                        r2 = random.nextDouble();
                        C = 2 * r2;
                        if (Math.abs(A) < 1) {
                            Xalpha = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(rdCounter);
                            Dalpha = C * Xalpha - x;
                            if (Dalpha < 0) {
                                Dalpha = Dalpha * -1;
                            }
                            x = Xalpha - A * Dalpha;
                        } else if (Math.abs(A) >= 1) {
                            randIndividual = random.nextInt(sortedFitnessValues.size());
                            Xrand = optimizationMatrix.get(randIndividual).get(rdCounter);
                            Drand = C * Xrand - x;
                            if (Drand < 0) {
                                Drand = Drand * -1;
                            }
                            x = Xrand - A * Drand;
                        }
                    } else {
                        Xalpha = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(rdCounter);
                        D = Xalpha - x;
                        if (D < 0) {
                            D = D * -1;
                        }
                        a2 = -1.0 + stCounter * ((-1.0) / iteration);
                        l = (a2 - 1.0) * random.nextDouble() + 1.0;
                        x = D * Math.pow(Math.E, (b * l)) * Math.cos(2 * Math.PI * l) + Xalpha;
                    }

                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    //Update Position
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Obstacle avoidance
                    ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                    for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                        isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        if (isPointInsideOfObstacle) {
                            ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                        }
                    }
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(sortedFitnessValues.get(0));
        }
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(0).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(0))
                    && positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(1).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }

        return sortedFitnessValues.get(0);
    }


    private double RLGWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        //Reinforcement Learning based Gray Wolf Optimization and Path Planning start here...
        //Reinforcement Learning Initialization Starts...
        ArrayList<ArrayList<Double>> QTable = new ArrayList<>();
        ArrayList<Double> vector = new ArrayList<>();

        vector.add(0.0);
        vector.add(0.0);
        QTable.add(vector);
        vector = new ArrayList<>();

        vector.add(0.0);
        vector.add(0.0);
        QTable.add(vector);

        double QValue, MaxQValue, reward;
        double alpha = 0.7, gamma = 0.8;
        //Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        double A, A1, A2, A3;
        double C1, C2, C3;
        double x, X1, X2, X3;
        double Dalpha, Dbeta, Ddelta;
        double Xalpha, Xbeta, Xdelta;
        double sigma1 = 0.1, sigma2 = 0.5, sigma3 = 0.9;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Gray Wolf Optimization iterations start here...
        System.out.println(sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * stCounter / iteration;
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                visitingStations = new ArrayList<>();
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

                    A = (A1 + A2 + A3) / 3;
                    if (A < 0) {
                        A = -1.0 * A;
                    }

                    if (A > 1) {
                        //State = exploration, work on 1st row of Q Table
                        if (QTable.get(0).get(0) > QTable.get(0).get(1)) {
                            //Action = exploration
                            QValue = QTable.get(0).get(0);
                            x = (X1 + X2 + X3) / 3;
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(0).get(0);
                            if (QTable.get(0).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(0).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(0).set(0, QValue);
                        } else {
                            //Action = exploitation
                            QValue = QTable.get(0).get(1);
                            x = X1 * sigma1 + X2 * sigma2 + X3 * sigma3;
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(1).get(0);
                            if (QTable.get(1).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(1).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(0).set(1, QValue);
                        }
                    } else {
                        //State = exploitation, work on 2nd row of Q Table
                        if (QTable.get(1).get(0) > QTable.get(1).get(1)) {
                            //Action = exploration
                            QValue = QTable.get(1).get(0);
                            x = (X1 + X2 + X3) / 3;
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(0).get(0);
                            if (QTable.get(0).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(0).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(1).set(0, QValue);
                        } else {
                            //Action = exploitation
                            QValue = QTable.get(1).get(1);
                            x = X1 * sigma1 + X2 * sigma2 + X3 * sigma3;
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(1).get(0);
                            if (QTable.get(1).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(1).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(1).set(1, QValue);
                        }
                    }

                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    //Update Position
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Obstacle avoidance
                    ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                    for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                        isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        if (isPointInsideOfObstacle) {
                            ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                        }
                    }
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(sortedFitnessValues.get(0));
        }
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(0).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(0))
                    && positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(1).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }

        return sortedFitnessValues.get(0);
    }


    private double RLIGWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        //Reinforcement Learning based Iterative Gray Wolf Optimization and Path Planning start here...
        //Reinforcement Learning Initialization Starts...
        ArrayList<ArrayList<Double>> QTable = new ArrayList<>();
        ArrayList<Double> vector = new ArrayList<>();

        vector.add(0.0);
        vector.add(0.0);
        QTable.add(vector);
        vector = new ArrayList<>();

        vector.add(0.0);
        vector.add(0.0);
        QTable.add(vector);

        double QValue, MaxQValue, reward;
        double alpha = 0.9, gamma = 0.8;
        ArrayList<Double> sigmas;
        double sum, averageA;
        //Iterative Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        ArrayList<Double> A;
        ArrayList<Double> C;
        ArrayList<Double> D;
        ArrayList<Double> X;
        double x;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Iterative Gray Wolf Optimization iterations start here...
        System.out.println(sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * stCounter / iteration;
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                visitingStations = new ArrayList<>();
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    A = new ArrayList<>();
                    C = new ArrayList<>();
                    D = new ArrayList<>();
                    X = new ArrayList<>();
                    sigmas = new ArrayList<>();
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

                    sum = 0;
                    for (double counter = 0; counter < ndCounter + 1; counter = counter + 1) {
                        sigmas.add(counter + 1);
                        sum = sum + (counter + 1);
                    }
                    for (int counter = 0; counter < sigmas.size(); counter = counter + 1) {
                        sigmas.set(counter, sigmas.get(counter) / sum);
                    }

                    averageA = 0;
                    for (int counter = 0; counter < A.size(); counter = counter + 1) {
                        averageA = averageA + A.get(counter);
                    }
                    averageA = averageA / A.size();

                    if (averageA > 1) {
                        //State = exploration, work on 1st row of Q Table
                        if (QTable.get(0).get(0) > QTable.get(0).get(1)) {
                            //Action = exploration
                            QValue = QTable.get(0).get(0);
                            x = 0;
                            for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                                x = x + X.get(fifthCounter);
                            }
                            x = x / X.size();
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(0).get(0);
                            if (QTable.get(0).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(0).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(0).set(0, QValue);
                        } else {
                            //Action = exploitation
                            QValue = QTable.get(0).get(1);
                            x = 0;
                            for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                                x = x + X.get(fifthCounter) * A.get(fifthCounter);
                            }
                            x = x / X.size();
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(1).get(0);
                            if (QTable.get(1).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(1).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(0).set(1, QValue);
                        }
                    } else {
                        //State = exploitation, work on 2nd row of Q Table
                        if (QTable.get(1).get(0) > QTable.get(1).get(1)) {
                            //Action = exploration
                            QValue = QTable.get(1).get(0);
                            x = 0;
                            for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                                x = x + X.get(fifthCounter);
                            }
                            x = x / X.size();
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(0).get(0);
                            if (QTable.get(0).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(0).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(1).set(0, QValue);
                        } else {
                            //Action = exploitation
                            QValue = QTable.get(1).get(1);
                            x = 0;
                            for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                                x = x + X.get(fifthCounter) * A.get(fifthCounter);
                            }
                            x = x / X.size();
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(1).get(0);
                            if (QTable.get(1).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(1).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(1).set(1, QValue);
                        }
                    }

                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    //Update Position
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Obstacle avoidance
                    ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                    for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                        isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        if (isPointInsideOfObstacle) {
                            ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                        }
                    }
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(sortedFitnessValues.get(0));
        }
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(0).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(0))
                    && positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(1).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }

        return sortedFitnessValues.get(0);
    }

    private double RLExGWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        //Reinforcement Learning based Expanded Gray Wolf Optimization and Path Planning start here...
        //Reinforcement Learning Initialization Starts...
        ArrayList<ArrayList<Double>> QTable = new ArrayList<>();
        ArrayList<Double> vector = new ArrayList<>();

        vector.add(0.0);
        vector.add(0.0);
        QTable.add(vector);
        vector = new ArrayList<>();

        vector.add(0.0);
        vector.add(0.0);
        QTable.add(vector);

        double QValue, MaxQValue, reward;
        double alpha = 0.9, gamma = 0.8;
        ArrayList<Double> sigmas;
        double sum, averageA;
        //Expanded Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        ArrayList<Double> A;
        ArrayList<Double> C;
        ArrayList<Double> D;
        ArrayList<Double> X;
        double x;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Expanded Gray Wolf Optimization iterations start here...
        System.out.println(sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a = 2.0 - 2.0 * Math.pow(stCounter, 2) / Math.pow(iteration, 2);
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                visitingStations = new ArrayList<>();
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    A = new ArrayList<>();
                    C = new ArrayList<>();
                    D = new ArrayList<>();
                    X = new ArrayList<>();
                    sigmas = new ArrayList<>();
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

                    sum = 0;
                    for (double counter = 0; counter < ndCounter + 1; counter = counter + 1) {
                        sigmas.add(counter + 1);
                        sum = sum + (counter + 1);
                    }
                    for (int counter = 0; counter < sigmas.size(); counter = counter + 1) {
                        sigmas.set(counter, sigmas.get(counter) / sum);
                    }

                    averageA = 0;
                    for (int counter = 0; counter < A.size(); counter = counter + 1) {
                        averageA = averageA + A.get(counter);
                    }
                    averageA = averageA / A.size();

                    if (averageA > 1) {
                        //State = exploration, work on 1st row of Q Table
                        if (QTable.get(0).get(0) > QTable.get(0).get(1)) {
                            //Action = exploration
                            QValue = QTable.get(0).get(0);
                            x = 0;
                            for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                                x = x + X.get(fifthCounter);
                            }
                            x = x / X.size();
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(0).get(0);
                            if (QTable.get(0).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(0).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(0).set(0, QValue);
                        } else {
                            //Action = exploitation
                            QValue = QTable.get(0).get(1);
                            x = 0;
                            for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                                x = x + X.get(fifthCounter) * A.get(fifthCounter);
                            }
                            x = x / X.size();
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(1).get(0);
                            if (QTable.get(1).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(1).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(0).set(1, QValue);
                        }
                    } else {
                        //State = exploitation, work on 2nd row of Q Table
                        if (QTable.get(1).get(0) > QTable.get(1).get(1)) {
                            //Action = exploration
                            QValue = QTable.get(1).get(0);
                            x = 0;
                            for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                                x = x + X.get(fifthCounter);
                            }
                            x = x / X.size();
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(0).get(0);
                            if (QTable.get(0).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(0).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(1).set(0, QValue);
                        } else {
                            //Action = exploitation
                            QValue = QTable.get(1).get(1);
                            x = 0;
                            for (int fifthCounter = 0; fifthCounter < X.size(); fifthCounter = fifthCounter + 1) {
                                x = x + X.get(fifthCounter) * A.get(fifthCounter);
                            }
                            x = x / X.size();
                            if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                reward = 1.0;
                            } else {
                                reward = -1.0;
                            }
                            //Finds Q max
                            MaxQValue = QTable.get(1).get(0);
                            if (QTable.get(1).get(1) > MaxQValue) {
                                MaxQValue = QTable.get(1).get(1);
                            }
                            //Calculate Q and update Q Table
                            QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                            QTable.get(1).set(1, QValue);
                        }
                    }

                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    //Update Position
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Obstacle avoidance
                    ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                    for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                        isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        if (isPointInsideOfObstacle) {
                            ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                        }
                    }
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(sortedFitnessValues.get(0));
        }
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(0).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(0))
                    && positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(1).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }

        return sortedFitnessValues.get(0);
    }


    private double RLWOA(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        //Reinforcement Learning based Whale Optimization Algorithm and Path Planning start here...
        //Reinforcement Learning Initialization Starts...
        ArrayList<ArrayList<Double>> QTable = new ArrayList<>();
        ArrayList<Double> vector = new ArrayList<>();

        vector.add(0.0);
        vector.add(0.0);
        QTable.add(vector);
        vector = new ArrayList<>();

        vector.add(0.0);
        vector.add(0.0);
        QTable.add(vector);

        double QValue, MaxQValue, reward;
        double alpha = 0.9, gamma = 0.8;
        //Whale Optimization Algorithm initialization starts here...
        double a1, a2, r1, r2, A, C;
        boolean p;
        int randIndividual;
        double Xalpha, Dalpha, Xrand, Drand;
        ArrayList<Double> P;
        double D;
        double b = 1, l;
        double x;
        double sigma1 = 0.6;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Whale Optimization Algorithm iterations start here...
        System.out.println(sortedFitnessValues.get(0));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
            positionsMatrixWithCollisions = new ArrayList<>();
            a1 = 2.0 - 2.0 * stCounter / iteration;
            for (int ndCounter = 0; ndCounter < optimizationMatrix.size(); ndCounter = ndCounter + 1) {
                visitingStations = new ArrayList<>();
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    x = optimizationMatrix.get(ndCounter).get(rdCounter);

                    r1 = random.nextDouble();
                    A = 2 * a1 * r1 - a1;
                    r2 = random.nextDouble();
                    C = 2 * r2;

                    p = random.nextBoolean();
                    if (p) {
                        if (Math.abs(A) < 1) {
                            Xalpha = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(rdCounter);
                            Dalpha = C * Xalpha - x;
                            if (Dalpha < 0) {
                                Dalpha = Dalpha * -1;
                            }

                            //State = exploitation, work on 2nd row of Q Table
                            if (QTable.get(1).get(0) > QTable.get(1).get(1)) {
                                //Action = exploration
                                QValue = QTable.get(1).get(0);
                                x = Xalpha - A * Dalpha;
                                if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                    reward = 1.0;
                                } else {
                                    reward = -1.0;
                                }
                                //Finds Q max
                                MaxQValue = QTable.get(0).get(0);
                                if (QTable.get(0).get(1) > MaxQValue) {
                                    MaxQValue = QTable.get(0).get(1);
                                }
                                //Calculate Q and update Q Table
                                QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                                QTable.get(1).set(0, QValue);
                            } else {
                                //Action = exploitation
                                QValue = QTable.get(1).get(1);
                                x = sigma1 * Xalpha - A * Dalpha;
                                if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                    reward = 1.0;
                                } else {
                                    reward = -1.0;
                                }
                                //Finds Q max
                                MaxQValue = QTable.get(1).get(0);
                                if (QTable.get(1).get(1) > MaxQValue) {
                                    MaxQValue = QTable.get(1).get(1);
                                }
                                //Calculate Q and update Q Table
                                QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                                QTable.get(1).set(1, QValue);
                            }
                        } else if (Math.abs(A) >= 1) {
                            randIndividual = random.nextInt(sortedFitnessValues.size());
                            Xrand = optimizationMatrix.get(randIndividual).get(rdCounter);
                            Drand = C * Xrand - x;
                            if (Drand < 0) {
                                Drand = Drand * -1;
                            }

                            //State = exploration, work on 1st row of Q Table
                            if (QTable.get(0).get(0) > QTable.get(0).get(1)) {
                                //Action = exploration
                                QValue = QTable.get(0).get(0);
                                x = Xrand - A * Drand;

                                if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                    reward = 1.0;
                                } else {
                                    reward = -1.0;
                                }
                                //Finds Q max
                                MaxQValue = QTable.get(0).get(0);
                                if (QTable.get(0).get(1) > MaxQValue) {
                                    MaxQValue = QTable.get(0).get(1);
                                }
                                //Calculate Q and update Q Table
                                QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                                QTable.get(0).set(0, QValue);
                            } else {
                                //Action = exploitation
                                QValue = QTable.get(0).get(0);
                                x = sigma1 * Xrand - A * Drand;

                                if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                    reward = 1.0;
                                } else {
                                    reward = -1.0;
                                }
                                //Finds Q max
                                MaxQValue = QTable.get(0).get(0);
                                if (QTable.get(0).get(1) > MaxQValue) {
                                    MaxQValue = QTable.get(0).get(1);
                                }
                                //Calculate Q and update Q Table
                                QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                                QTable.get(0).set(0, QValue);
                            }
                        }
                    } else {
                        Xalpha = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(rdCounter);
                        D = Xalpha - x;
                        if (D < 0) {
                            D = D * -1;
                        }
                        a2 = -1.0 + stCounter * ((-1.0) / iteration);
                        l = (a2 - 1.0) * random.nextDouble() + 1.0;

                        if (A > 1) {
                            //State = exploration, work on 1st row of Q Table
                            if (QTable.get(0).get(0) > QTable.get(0).get(1)) {
                                //Action = exploration
                                QValue = QTable.get(0).get(0);
                                x = D * Math.pow(Math.E, (b * l)) * Math.cos(2 * Math.PI * l) + Xalpha;
                                if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                    reward = 1.0;
                                } else {
                                    reward = -1.0;
                                }
                                //Finds Q max
                                MaxQValue = QTable.get(0).get(0);
                                if (QTable.get(0).get(1) > MaxQValue) {
                                    MaxQValue = QTable.get(0).get(1);
                                }
                                //Calculate Q and update Q Table
                                QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                                QTable.get(0).set(0, QValue);
                            } else {
                                //Action = exploitation
                                QValue = QTable.get(0).get(1);
                                x = sigma1 * D * Math.pow(Math.E, (b * l)) * Math.cos(2 * Math.PI * l) + Xalpha;
                                if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                    reward = 1.0;
                                } else {
                                    reward = -1.0;
                                }
                                //Finds Q max
                                MaxQValue = QTable.get(1).get(0);
                                if (QTable.get(1).get(1) > MaxQValue) {
                                    MaxQValue = QTable.get(1).get(1);
                                }
                                //Calculate Q and update Q Table
                                QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                                QTable.get(0).set(1, QValue);
                            }
                        } else {
                            //State = exploitation, work on 2nd row of Q Table
                            if (QTable.get(1).get(0) > QTable.get(1).get(1)) {
                                //Action = exploration
                                QValue = QTable.get(1).get(0);
                                x = D * Math.pow(Math.E, (b * l)) * Math.cos(2 * Math.PI * l) + Xalpha;
                                if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                    reward = 1.0;
                                } else {
                                    reward = -1.0;
                                }
                                //Finds Q max
                                MaxQValue = QTable.get(0).get(0);
                                if (QTable.get(0).get(1) > MaxQValue) {
                                    MaxQValue = QTable.get(0).get(1);
                                }
                                //Calculate Q and update Q Table
                                QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                                QTable.get(1).set(0, QValue);
                            } else {
                                //Action = exploitation
                                QValue = QTable.get(1).get(1);
                                x = sigma1 * D * Math.pow(Math.E, (b * l)) * Math.cos(2 * Math.PI * l) + Xalpha;
                                if (x < optimizationMatrix.get(ndCounter).get(rdCounter)) {
                                    reward = 1.0;
                                } else {
                                    reward = -1.0;
                                }
                                //Finds Q max
                                MaxQValue = QTable.get(1).get(0);
                                if (QTable.get(1).get(1) > MaxQValue) {
                                    MaxQValue = QTable.get(1).get(1);
                                }
                                //Calculate Q and update Q Table
                                QValue = QValue + alpha * (reward + gamma * MaxQValue - QValue);
                                QTable.get(1).set(1, QValue);
                            }
                        }
                    }

                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    //Update Position
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Obstacle avoidance
                    ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                    sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                    for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
                        isPointInsideOfObstacle = obstacleAvoider.isPointInsideOfObstacle(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                        if (isPointInsideOfObstacle) {
                            ObstacleAvoidanceCurrentStation = obstacleAvoider.findNearestCorner(ObstacleAvoidanceCurrentStation, sortedObstacles.get(fourthCounter));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(0, ObstacleAvoidanceCurrentStation.get(0));
                            positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter).set(1, ObstacleAvoidanceCurrentStation.get(1));
                        }
                    }
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    if (rdCounter > 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1);
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                        if (rdCounter == optimizationMatrix.get(ndCounter).size() - 1) {
                            ObstacleAvoidanceCurrentStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                            ObstacleAvoidanceNextStation = destinationStation;
                            sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                            for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                    } else if (rdCounter == 0) {
                        //Obstacle avoidance
                        ObstacleAvoidanceCurrentStation = sourceStation;
                        ObstacleAvoidanceNextStation = positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter);
                        sortedObstacles = new ArrayList<>(obstacles.getSortedObstacles(ObstacleAvoidanceCurrentStation));
                        for (int fourthCounter = 0; fourthCounter < sortedObstacles.size(); fourthCounter = fourthCounter + 1) {
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
                        pathWithCollisiions.add(ObstacleAvoidanceNextStation);
                    }
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(sortedFitnessValues.get(0));
        }
        for (int counter = 0; counter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size() - 1; counter = counter + 1) {
            if (positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(0).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(0))
                    && positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter).get(1).equals(positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(counter + 1).get(1))) {
                positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).remove(counter + 1);
                counter = counter - 1;
            }
        }
        System.out.println("Alpha path with collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }

        return sortedFitnessValues.get(0);
    }


    private ArrayList<ArrayList<Double>> createRandomStations(int theNumberOfStations, ArrayList<Double> Xboundaries,
                                                              ArrayList<Double> Yboundaries,
                                                              ArrayList<Double> Zboundaries) {
        random = new Random();
        ArrayList<ArrayList<Double>> stations = new ArrayList<>();
        ArrayList<Double> vector = new ArrayList<>();
        double x, y, z;

        for (int stCounter = 0; stCounter < theNumberOfStations; stCounter = stCounter + 1) {
            x = Xboundaries.get(0) + (Xboundaries.get(1) - Xboundaries.get(0)) * random.nextDouble();
            y = Yboundaries.get(0) + (Yboundaries.get(1) - Yboundaries.get(0)) * random.nextDouble();
            z = Zboundaries.get(0) + (Zboundaries.get(1) - Zboundaries.get(0)) * random.nextDouble();
            vector.add(x);
            vector.add(y);
            vector.add(z);

            stations.add(vector);
            vector = new ArrayList<>();
        }

        return stations;
    }

    private ArrayList<ArrayList<ArrayList<Double>>> createRandomVisitedStations(int population, int dimension, ArrayList<ArrayList<Double>> stations) {
        random = new Random();
        ArrayList<ArrayList<ArrayList<Double>>> visitedStationsMatrix = new ArrayList<>();
        ArrayList<ArrayList<Double>> visitedStationsVector = new ArrayList<>();
        int index;

        for (int stCounter = 0; stCounter < population; stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < dimension; ndCounter = ndCounter + 1) {
                index = 0 + random.nextInt(((stations.size() - 1) - 0) + 1);
                visitedStationsVector.add(stations.get(index));
            }
            visitedStationsMatrix.add(visitedStationsVector);
            visitedStationsVector = new ArrayList<>();
        }

        return visitedStationsMatrix;
    }

    private ArrayList<ArrayList<Double>> createOptimizationMatrix(ArrayList<ArrayList<ArrayList<Double>>> visitedStationsMatrix, ArrayList<Double> sourceStation) {
        ArrayList<ArrayList<Double>> optimizationMatrix = new ArrayList<>();
        ArrayList<Double> optimizationVector = new ArrayList<>();
        double euclideanDistance;

        for (int stCounter = 0; stCounter < visitedStationsMatrix.size(); stCounter = stCounter + 1) {
            for (int ndCounter = 0; ndCounter < visitedStationsMatrix.get(stCounter).size(); ndCounter = ndCounter + 1) {
                if (ndCounter == 0) {
                    euclideanDistance = findEuclideanDistance(sourceStation, visitedStationsMatrix.get(stCounter).get(ndCounter));
                    optimizationVector.add(euclideanDistance);
                } else {
                    euclideanDistance = findEuclideanDistance(visitedStationsMatrix.get(stCounter).get(ndCounter - 1), visitedStationsMatrix.get(stCounter).get(ndCounter));
                    optimizationVector.add(euclideanDistance);
                }
            }
            optimizationMatrix.add(optimizationVector);
            optimizationVector = new ArrayList<>();
        }

        return optimizationMatrix;
    }

    private ArrayList<Double> findFitnessValues(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions, ArrayList<Double> sourceStation, ArrayList<Double> destinationStation) {
        ArrayList<Double> fitnessValues = new ArrayList<>();
        double sum = 0.0;
        double euclideanDistance;

        for (int stCounter = 0; stCounter < positionsMatrixWithCollisions.size(); stCounter = stCounter + 1) {
            euclideanDistance = findEuclideanDistance(sourceStation, positionsMatrixWithCollisions.get(stCounter).get(0));
            sum = sum + euclideanDistance;

            euclideanDistance = findEuclideanDistance(positionsMatrixWithCollisions.get(stCounter).get(positionsMatrixWithCollisions.get(stCounter).size() - 1), destinationStation);
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
