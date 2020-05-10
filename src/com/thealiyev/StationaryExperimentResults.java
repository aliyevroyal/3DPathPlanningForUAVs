package com.thealiyev;

import java.util.ArrayList;
import java.util.Random;

public class StationaryExperimentResults {
    private static Random random = null;

    public static void main(String[] args) {
        //Boundaries of map
        ArrayList<Double> xBoundaries = new ArrayList<>(), yBoundaries = new ArrayList<>(), zBoundaries = new ArrayList<>();
        //X boundaries
        xBoundaries.add(0.0);
        xBoundaries.add(50.0);
        //Y boundaries
        yBoundaries.add(0.0);
        yBoundaries.add(50.0);
        //Z boundaries
        zBoundaries.add(0.0);
        zBoundaries.add(50.0);
        StationaryExperimentResults stationaryExperimentResults = new StationaryExperimentResults();
        int theNumberOfStations = 1000;
        int population = 100, dimension = 5;
        ArrayList<ArrayList<Double>> stations = stationaryExperimentResults.createRandomStations(theNumberOfStations, xBoundaries, yBoundaries, zBoundaries);
        ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions = stationaryExperimentResults.createRandomVisitedStations(population, dimension, stations);
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

        StationaryGWOExperimentResults stationaryGWOExperimentResults = new StationaryGWOExperimentResults();
        stationaryGWOExperimentResults.GWO(stations, GWOPositionsMatrixWithoutCollisions);
        StationaryIGWOExperimentResults stationaryIGWOExperimentResults = new StationaryIGWOExperimentResults();
        stationaryIGWOExperimentResults.IGWO(stations, IGWOPositionsMatrixWithoutCollisions);
        StationaryExGWOExperimentResults stationaryExGWOExperimentResults = new StationaryExGWOExperimentResults();
        stationaryExGWOExperimentResults.ExGWO(stations, ExGWOPositionsMatrixWithoutCollisions);
    }

    public ArrayList<ArrayList<Double>> createRandomStations(int theNumberOfStations, ArrayList<Double> xBoundaries,
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

    public ArrayList<ArrayList<ArrayList<Double>>> createRandomVisitedStations(int population, int dimension, ArrayList<ArrayList<Double>> stations) {
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

    public ArrayList<ArrayList<Double>> createOptimizationMatrix(ArrayList<ArrayList<ArrayList<Double>>> visitedStationsMatrix, ArrayList<Double> sourceStation) {
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

    public ArrayList<Double> findFitnessValues(ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithCollisions, ArrayList<Double> sourceStation, ArrayList<Double> destinationStation) {
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

class StationaryGWOExperimentResults {
    private Random random = null;
    private StationaryExperimentResults stationaryExperimentResults = null;

    public void GWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        stationaryExperimentResults = new StationaryExperimentResults();
        //Boundaries of map
        ArrayList<Double> xBoundaries = new ArrayList<>(), yBoundaries = new ArrayList<>(), zBoundaries = new ArrayList<>();
        //X boundaries
        xBoundaries.add(0.0);
        xBoundaries.add(50.0);
        //Y boundaries
        yBoundaries.add(0.0);
        yBoundaries.add(50.0);
        //Z boundaries
        zBoundaries.add(0.0);
        zBoundaries.add(50.0);
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
        double X, X1, X2, X3;
        double Dalpha, Dbeta, Ddelta;
        double Xalpha, Xbeta, Xdelta;
        int theNumberOfStations = 1000;
        int population = 100, dimension = 5;
        int iteration = 100;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = stationaryExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = stationaryExperimentResults.findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = stationaryExperimentResults.sortFitnessValues(fitnessValues);
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
                visitingStations = new ArrayList<>();
                for (int rdCounter = 0; rdCounter < optimizationMatrix.get(ndCounter).size(); rdCounter = rdCounter + 1) {
                    X = optimizationMatrix.get(ndCounter).get(rdCounter);

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A1 = 2 * a * r1 - a;
                    C1 = 2 * r2;
                    Xalpha = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))).get(rdCounter);
                    Dalpha = C1 * Xalpha - X;
                    if (Dalpha < 0) {
                        Dalpha = Dalpha * -1;
                    }
                    X1 = Xalpha - A1 * Dalpha;

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A2 = 2 * a * r1 - a;
                    C2 = 2 * r2;
                    Xbeta = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(1))).get(rdCounter);
                    Dbeta = C2 * Xbeta - X;
                    if (Dbeta < 0) {
                        Dbeta = Dbeta * -1;
                    }
                    X2 = Xbeta - A2 * Dbeta;

                    r1 = random.nextDouble();
                    r2 = random.nextDouble();
                    A3 = 2 * a * r1 - a;
                    C3 = 2 * r2;
                    Xdelta = optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(2))).get(rdCounter);
                    Ddelta = C3 * Xdelta - X;
                    if (Ddelta < 0) {
                        Ddelta = Ddelta * -1;
                    }
                    X3 = Xdelta - A3 * Ddelta;

                    X = (X1 + X2 + X3) / 3;
                    for (int fourthCounter = 0; fourthCounter < stations.size(); fourthCounter = fourthCounter + 1) {
                        if (!visitingStations.contains(stations.get(fourthCounter))) {
                            if (rdCounter == 0) {
                                distance = stationaryExperimentResults.findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + stationaryExperimentResults.findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = stationaryExperimentResults.findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + stationaryExperimentResults.findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < X && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Update Position
                    if (rdCounter > 0) {
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
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = stationaryExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = stationaryExperimentResults.findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = stationaryExperimentResults.sortFitnessValues(fitnessValues);
        }
        System.out.println("GWO Experiment Results: ");
        System.out.println("After iterations, alpha's fitness value:" + sortedFitnessValues.get(0));
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

class StationaryIGWOExperimentResults {
    private Random random = null;
    private StationaryExperimentResults stationaryExperimentResults = null;

    public void IGWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        stationaryExperimentResults = new StationaryExperimentResults();
        //Boundaries of map
        ArrayList<Double> xBoundaries = new ArrayList<>(), yBoundaries = new ArrayList<>(), zBoundaries = new ArrayList<>();
        //X boundaries
        xBoundaries.add(0.0);
        xBoundaries.add(50.0);
        //Y boundaries
        yBoundaries.add(0.0);
        yBoundaries.add(50.0);
        //Z boundaries
        zBoundaries.add(0.0);
        zBoundaries.add(50.0);
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
        int theNumberOfStations = 1000;
        int population = 100, dimension = 5;
        int iteration = 100;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = stationaryExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = stationaryExperimentResults.findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = stationaryExperimentResults.sortFitnessValues(fitnessValues);
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
                                distance = stationaryExperimentResults.findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + stationaryExperimentResults.findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = stationaryExperimentResults.findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + stationaryExperimentResults.findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Update Position
                    if (rdCounter > 0) {
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
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = stationaryExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = stationaryExperimentResults.findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = stationaryExperimentResults.sortFitnessValues(fitnessValues);
        }
        System.out.println("IGWO Experiment Results: ");
        System.out.println("After iterations, alpha's fitness value:" + sortedFitnessValues.get(0));
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

class StationaryExGWOExperimentResults {
    private Random random = null;
    private StationaryExperimentResults stationaryExperimentResults = null;

    public void ExGWO(ArrayList<ArrayList<Double>> stations, ArrayList<ArrayList<ArrayList<Double>>> positionsMatrixWithoutCollisions) {
        random = new Random();
        stationaryExperimentResults = new StationaryExperimentResults();
        //Boundaries of map
        ArrayList<Double> xBoundaries = new ArrayList<>(), yBoundaries = new ArrayList<>(), zBoundaries = new ArrayList<>();
        //X boundaries
        xBoundaries.add(0.0);
        xBoundaries.add(50.0);
        //Y boundaries
        yBoundaries.add(0.0);
        yBoundaries.add(50.0);
        //Z boundaries
        zBoundaries.add(0.0);
        zBoundaries.add(50.0);
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
        int theNumberOfStations = 1000;
        int population = 100, dimension = 5;
        int iteration = 100;
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<Double>> optimizationMatrix = stationaryExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
        ArrayList<Double> fitnessValues = stationaryExperimentResults.findFitnessValues(positionsMatrixWithoutCollisions, sourceStation, destinationStation);
        ArrayList<Double> sortedFitnessValues = stationaryExperimentResults.sortFitnessValues(fitnessValues);
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
                                distance = stationaryExperimentResults.findEuclideanDistance(sourceStation, stations.get(fourthCounter)) + stationaryExperimentResults.findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            } else {
                                distance = stationaryExperimentResults.findEuclideanDistance(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter)) + stationaryExperimentResults.findEuclideanDistance(destinationStation, stations.get(fourthCounter));
                            }
                            if (distance < x && distance > 0) {
                                distances.add(distance);
                                indexes.add(fourthCounter);
                            }
                        }
                    }
                    if (distances.size() > 0) {
                        positionsMatrixWithoutCollisions.get(ndCounter).set(rdCounter, stations.get(indexes.get(random.nextInt((indexes.size() - 1) + 1))));
                    }
                    //Update Position
                    if (rdCounter > 0) {
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
                    visitingStations.add(positionsMatrixWithoutCollisions.get(ndCounter).get(rdCounter));
                    distances = new ArrayList<>();
                    indexes = new ArrayList<>();
                    optimizationMatrix = stationaryExperimentResults.createOptimizationMatrix(positionsMatrixWithoutCollisions, sourceStation);
                }
                positionsMatrixWithCollisions.add(pathWithCollisiions);
                pathWithCollisiions = new ArrayList<>();
            }
            fitnessValues = stationaryExperimentResults.findFitnessValues(positionsMatrixWithCollisions, sourceStation, destinationStation);
            sortedFitnessValues = stationaryExperimentResults.sortFitnessValues(fitnessValues);
        }
        System.out.println("ExGWO Experiment Results: ");
        System.out.println("After iterations, alpha's fitness value:" + sortedFitnessValues.get(0));
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