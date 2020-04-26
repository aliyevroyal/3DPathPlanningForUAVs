package com.thealiyev;

import java.util.ArrayList;
import java.util.Random;

public class IGWOPathPlanningStaticStationaryMethod {
    private static Random random = null;

    public static void main(String[] args) {
        IGWOPathPlanningStaticStationaryMethod igwoPathPlanningStaticStationaryMethod = new IGWOPathPlanningStaticStationaryMethod();
        igwoPathPlanningStaticStationaryMethod.GWO();
    }

    private void GWO() {
        random = new Random();
        //Boundaries of map
        ArrayList<Double> xBoundaries = new ArrayList<>(), yBoundaries = new ArrayList<>(), zBoundaries = new ArrayList<>();
        //X boundaries
        xBoundaries.add(0.0);
        xBoundaries.add(100.0);
        //Y boundaries
        yBoundaries.add(0.0);
        yBoundaries.add(100.0);
        //Z boundaries
        zBoundaries.add(0.0);
        zBoundaries.add(100.0);
        //Source station coordinates
        ArrayList<Double> sourceStation = new ArrayList<>();
        sourceStation.add(1.0);
        sourceStation.add(1.0);
        sourceStation.add(1.0);
        //Destination station coordinates
        ArrayList<Double> destinationStation = new ArrayList<>();
        destinationStation.add(99.0);
        destinationStation.add(99.0);
        destinationStation.add(99.0);
        //Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        ArrayList<Double> A;
        ArrayList<Double> C;
        ArrayList<Double> D;
        ArrayList<Double> X;
        double x;
        int theNumberOfStations = 1000;
        int population = 100, dimension = 30;
        int iteration = 100;
        ArrayList<ArrayList<Double>> stations = createRandomStations(theNumberOfStations, xBoundaries, yBoundaries, zBoundaries);
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<ArrayList<Double>>> visitedStationsMatrix = createRandomVisitedStations(population, dimension, stations);
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(visitedStationsMatrix, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(visitedStationsMatrix, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance = 0;
        //Obstacle detection, checking and avoidance components initialization starts here...
        ObstacleAvoider obstacleAvoider = new ObstacleAvoider();
        ArrayList<Double> ObstacleAvoidanceCurrentStation;
        ArrayList<Double> ObstacleAvoidanceNextStation;
        Obstacles obstacles = new Obstacles();
        obstacles.setObstacle1();
        obstacles.setObstacle2();
        obstacles.setObstacle3();
        obstacles.setObstacle4();
        obstacles.setObstacle5();
        obstacles.setObstacle6();
        obstacles.setObstacle7();
        obstacles.setObstacle8();
        boolean isPointInsideOfObstacle;
        boolean didPointCollideWithObstacle;
        ArrayList<Double> nearestCornerToCurrentStation;
        ArrayList<Double> newCornerForNextStation;
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pathWithCollisiions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> sortedObstacles;
        //Gray Wolf Optimization iterations start here...
        System.out.println("Initialization, alpha's fitness value: " + sortedFitnessValues.get(0));
        System.out.println("Initialization, alpha values: " + optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
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
                                distance = findEuclideanDistance(visitedStationsMatrix.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    if (distances.size() > 0) {
                        visitedStationsMatrix.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }

                    visitingStations.add(visitedStationsMatrix.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = createOptimizationMatrix(visitedStationsMatrix, sourceStation);
                }
            }
            fitnessValues = findFitnessValues(visitedStationsMatrix, sourceStation, destinationStation);
            sortedFitnessValues = sortFitnessValues(fitnessValues);
            System.out.println(stCounter + " iteration, alpha's fitness value: " + sortedFitnessValues.get(0));
            //System.out.println(stCounter + " iteration, alpha's values: " + optimizationMatrix.get(fitnessValues.indexOf(trio.get(0))));
        }
        System.out.println("After iterations, alpha's fitness value:" + sortedFitnessValues.get(0));
        System.out.println("After iterations, alpha's values: " + optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))));
        for (int ndCounter = 0; ndCounter < visitedStationsMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).size(); ndCounter = ndCounter + 1) {
            System.out.println(ndCounter + " " + visitedStationsMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(ndCounter));
        }
    }

    private ArrayList<ArrayList<Double>> createRandomStations(int theNumberOfStations, ArrayList<Double> xBoundaries,
                                                              ArrayList<Double> yBoundaries,
                                                              ArrayList<Double> zBoundaries) {
        random = new Random();
        ArrayList<ArrayList<Double>> stations = new ArrayList<>();
        ArrayList<Double> vector = new ArrayList<>();
        double x, y, z;

        for (int stCounter = 0; stCounter < theNumberOfStations; stCounter = stCounter + 1) {
            x = xBoundaries.get(0) + (xBoundaries.get(1) - xBoundaries.get(0)) * random.nextDouble();
            y = yBoundaries.get(0) + (yBoundaries.get(1) - yBoundaries.get(0)) * random.nextDouble();
            z = zBoundaries.get(0) + (zBoundaries.get(1) - zBoundaries.get(0)) * random.nextDouble();
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

    private ArrayList<Double> findFitnessValues(ArrayList<ArrayList<ArrayList<Double>>> optimizationMatrix, ArrayList<Double> sourceStation, ArrayList<Double> destinatioStation) {
        ArrayList<Double> fitnessValues = new ArrayList<>();
        double sum = 0.0;
        double euclideanDistance;

        for (int stCounter = 0; stCounter < optimizationMatrix.size(); stCounter = stCounter + 1) {
            euclideanDistance = findEuclideanDistance(sourceStation, optimizationMatrix.get(stCounter).get(0));
            sum = sum + euclideanDistance;

            euclideanDistance = findEuclideanDistance(optimizationMatrix.get(stCounter).get(optimizationMatrix.get(stCounter).size() - 1), destinatioStation);
            sum = sum + euclideanDistance;

            for (int ndCounter = 0; ndCounter < optimizationMatrix.get(stCounter).size() - 1; ndCounter = ndCounter + 1) {
                euclideanDistance = findEuclideanDistance(optimizationMatrix.get(stCounter).get(ndCounter), optimizationMatrix.get(stCounter).get(ndCounter + 1));
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