package com.thealiyev;

import java.util.ArrayList;

public class Main {
    public static void main(String[] args) {
        ArrayList<ArrayList<Double>> path = new ArrayList<>();
        ArrayList<Double> coordinate = new ArrayList<>();
        ArrayList<Double> vector = new ArrayList<>();
        coordinate.add(60.0);
        coordinate.add(60.0);
        coordinate.add(87.38138045175235);
        path.add(coordinate);
        coordinate = new ArrayList<>();

        coordinate.add(60.0);
        coordinate.add(60.0);
        coordinate.add(87.38138045175235);
        path.add(coordinate);
        coordinate = new ArrayList<>();

        coordinate.add(60.0);
        coordinate.add(60.0);
        coordinate.add(87.38138045175235);
        path.add(coordinate);
        coordinate = new ArrayList<>();


        for (int counter = 0; counter < path.size(); counter = counter + 1) {
            System.out.print(path.get(counter) + " ");
        }
        System.out.println();

        vector = path.get(0);
        for (int counter = 1; counter < path.size(); counter = counter + 1) {
            if (path.get(counter).get(0).equals(path.get(counter - 1).get(0))) {
                if (path.get(counter).get(1).equals(path.get(counter - 1).get(1))) {
                    if (path.get(counter).get(2).equals(path.get(counter - 1).get(2))) {
                        path.remove(new ArrayList<>(path.get(counter)));
                    }
                }
            }
        }

        for (int counter = 0; counter < path.size(); counter = counter + 1) {
            System.out.print(path.get(counter) + " ");
        }
    }
}
