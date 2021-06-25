package daxing;

import daxing.expression.Expression;

public class Start {

    public static void main(String[] args) {
        String inputFile=args[0];
        double threshold=Double.parseDouble(args[1]);
        String outFile=args[2];
        Expression.calculateConnection(inputFile, threshold, outFile);
    }

}