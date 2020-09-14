package com.thealiyev;

import java.util.ArrayList;
import java.util.Random;

public class WOAPathPlanningStaticStationaryMethod {
    private static Random random = null;

    public static void main(String[] args) {
        WOAPathPlanningStaticStationaryMethod woaPathPlanningStaticStationaryMethod = new WOAPathPlanningStaticStationaryMethod();
        woaPathPlanningStaticStationaryMethod.WOA();
    }

    private void WOA() {
        random = new Random();
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
        //Whale Optimization Algorithm initialization starts here...
        double a1, a2, r1, r2, A, C;
        boolean p;
        int randIndividual;
        double Xalpha, Dalpha, Xrand, Drand;
        ArrayList<Double> P;
        double D;
        double b = 1, l;
        double x;
        int theNumberOfStations = 1000;
        int population = 100, dimension = 5;
        int iteration = 100;
        double sigma1 = 0.6;
        ArrayList<ArrayList<Double>> stations = createRandomStations(theNumberOfStations, Xboundaries, Yboundaries, Zboundaries);
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions = createRandomVisitedStations(population, dimension, stations);
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
        //Whale Optimization Algorithm iterations start here...
        System.out.println("Initialization, alpha's fitness value: " + sortedFitnessValues.get(0));
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
            System.out.println(stCounter + " iteration, alpha's fitness value: " + sortedFitnessValues.get(0));
        }
        System.out.println("After iterations, alpha's fitness value:" + sortedFitnessValues.get(0));
        System.out.println("Alpha path without collisions");
        for (int ndCounter = 0; ndCounter < positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + positionsMatrixWithoutCollisions.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
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