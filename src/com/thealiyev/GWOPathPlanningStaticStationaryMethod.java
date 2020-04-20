package com.thealiyev;

import java.util.ArrayList;
import java.util.Random;

public class GWOPathPlanningStaticStationaryMethod {
    private static Random random = null;

    public static void main(String[] args) {
        GWOPathPlanningStaticStationaryMethod gwoPathPlanningStaticStationaryMethod = new GWOPathPlanningStaticStationaryMethod();
        gwoPathPlanningStaticStationaryMethod.GWO();
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
        sourceStation.add(0.0);
        sourceStation.add(0.0);
        sourceStation.add(0.0);
        //Destination station coordinates
        ArrayList<Double> destinatioStation = new ArrayList<>();
        destinatioStation.add(100.0);
        destinatioStation.add(100.0);
        destinatioStation.add(100.0);

        //Gray Wolf Optimization initialization starts here...
        double a;
        double r1, r2;
        double A1, A2, A3;
        double C1, C2, C3;
        double X, X1, X2, X3;
        double Dalpha, Dbeta, Ddelta;
        double Xalpha, Xbeta, Xdelta;
        int theNumberOfStations = 1000;
        int population = 10, dimension = 15;
        int iteration = 100;
        ArrayList<ArrayList<Double>> stations = createRandomStations(theNumberOfStations, xBoundaries, yBoundaries, zBoundaries);
        ArrayList<ArrayList<Double>> visitingStations;
        ArrayList<ArrayList<ArrayList<Double>>> visitedStationsMatrix = createRandomVisitedStations(population, dimension, stations);
        ArrayList<ArrayList<Double>> optimizationMatrix = createOptimizationMatrix(visitedStationsMatrix, sourceStation);
        ArrayList<Double> fitnessValues = findFitnessValues(visitedStationsMatrix, sourceStation, destinatioStation);
        ArrayList<Double> sortedFitnessValues = sortFitnessValues(fitnessValues);
        ArrayList<Double> distances = new ArrayList<>();
        ArrayList<Integer> indexes = new ArrayList<>();
        double distance = 0;

        //Gray Wolf Optimization iterations start here...
        System.out.println("Initialization, alpha's fitness value: " + sortedFitnessValues.get(0));
        System.out.println("Initialization, alpha values: " + optimizationMatrix.get(fitnessValues.indexOf(sortedFitnessValues.get(0))));
        for (int stCounter = 0; stCounter < iteration; stCounter = stCounter + 1) {
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
                                distance = findEuclideanDistance(sourceStation, stations.get(fourthCounter));
                            } else {
                                distance = findEuclideanDistance(visitedStationsMatrix.get(ndCounter).get(rdCounter - 1), stations.get(fourthCounter));
                            }
                            if (distance < X && distance > 0) {
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
            fitnessValues = findFitnessValues(visitedStationsMatrix, sourceStation, destinatioStation);
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