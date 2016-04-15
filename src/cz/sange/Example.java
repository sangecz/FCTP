package cz.sange;

import gurobi.*;
import java.io.*;

/**
 * INPUT format
 3 8                                                        // numWarehouses numCustomers
 5000000 4000000 4500000                                    // fixed costs for warehouses
 1000000 800000 1250000                                    // capacities for warehouses
 200000 200000 200000 200000 250000 250000 250000 250000    // demands for customers
 4 5 5 4 4 4.2 3.3 5                                        // transport costs from warehouses to customers
 2.5 3.5 4.5 3 2.2 4 2.6 5
 2 4 5 2.5 2.6 3.8 2.9 5.5
 *
 */

public class Example {

    private int [] demands;
    private int [] capacities;
    private int [] fixedCosts;
    private double[][] transportCosts;

    public static void main(String[] args) {

        try {
            System.setIn(new FileInputStream(args[0]));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        Example example = new Example();

        try {
            example.read();
        } catch (IOException e) {
            e.printStackTrace();
        }

        example.solve();
    }

    private void read() throws IOException {
        final StreamTokenizer st = new StreamTokenizer(new BufferedReader(new InputStreamReader(System.in)));

        while (st.nextToken() != StreamTokenizer.TT_EOF) {

            // # warehouses
            int w = (int) st.nval;
            st.nextToken();
            // # customers
            int c = (int) st.nval;
            st.nextToken();

            // fixedCosts costs by warehouse
            fixedCosts = new int[w];
            for (int i = 0; i < w; i++) {
                fixedCosts[i] = (int) st.nval;
                st.nextToken();
            }

            // capacities by warehouse
            capacities = new int[w];
            for (int i = 0; i < w; i++) {
                capacities[i] = (int) st.nval;
                st.nextToken();
            }

            // demands by customer
            demands = new int[c];
            for (int i = 0; i < c; i++) {
                demands[i] = (int) st.nval;
                st.nextToken();
            }

            // trasport costs
            transportCosts = new double[w][c];
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < c; j++) {
                    transportCosts[i][j] = st.nval;
                    st.nextToken();
                }
            }
        }

    }

    private void solve() {
        try {

            int numWarehouses = capacities.length;
            int numCustomers = demands.length;

            //// model
            GRBEnv env = new GRBEnv("example.log") ;
            GRBModel model = new GRBModel(env) ;
            model.set(GRB.StringAttr.ModelName, "FCTP");

            //// vars
            GRBVar [] open = new GRBVar[numWarehouses];
            for (int p = 0; p < numWarehouses; ++p) {
                open[p] = model.addVar(0, 1, fixedCosts[p], GRB.BINARY, "Open" + p);
            }

            GRBVar [][] transport = new GRBVar[numWarehouses][numCustomers];
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    transport[w][c] =
                            model.addVar(0, GRB.INFINITY, transportCosts[w][c], GRB.CONTINUOUS, "trans(" + w + ", " + c + ")");
                }
            }

            //// minimize
            model.set(GRB.IntAttr.ModelSense, 1);

            //// update vars
            model.update();

            //// constraints
            // capacities constraints
            for (int w = 0; w < numWarehouses; ++w) {
                GRBLinExpr ptot = new GRBLinExpr();
                for (int c = 0; c < numCustomers; ++c) {
                    ptot.addTerm(1.0, transport[w][c]);
                }
                GRBLinExpr limit = new GRBLinExpr();
                limit.addTerm(capacities[w], open[w]);
                model.addConstr(ptot, GRB.LESS_EQUAL, limit, "Capacity " + w);
            }
            // demand constraints
            for (int c = 0; c < numCustomers; ++c) {
                GRBLinExpr dtot = new GRBLinExpr();
                for (int w = 0; w < numWarehouses; ++w) {
                    dtot.addTerm(1.0, transport[w][c]);
                }
                model.addConstr(dtot, GRB.GREATER_EQUAL, demands[c], "Demand" + c);
            }

            //// prep
            // open all paths
            for (int p = 0; p < numWarehouses; ++p) {
                open[p].set(GRB.DoubleAttr.Start, 1.0);
            }
            // find max cfixed cost
            System.out.println("Initial guess:");
            double maxFixed = -GRB.INFINITY;
            for (int w = 0; w < numWarehouses; ++w) {
                if (fixedCosts[w] > maxFixed){
                    maxFixed = fixedCosts[w];
                }
            }
            // close most expensive
            for (int w = 0; w < numWarehouses; ++w) {
                if (fixedCosts[w] == maxFixed) {
                    open[w].set(GRB.DoubleAttr.Start, 0.0);
                    System.out.println("Closing w: " + w + "\n");
                    break;
                }
            }

            // relaxation
            model.getEnv().set(GRB.IntParam.Method, GRB.METHOD_BARRIER);

            // solve
            model.optimize();

            // print result
            System.out.println("\nTOTAL COSTS: " + model.get(GRB.DoubleAttr.ObjVal));
            System.out.println("SOLUTION:");
            for (int w = 0; w < numWarehouses; ++w) {
                // if is open
                if(open[w].get(GRB.DoubleAttr.X) == 1.0) {
                    System.out.println("w: " + w + " open:");
                    for (int c = 0; c < numCustomers; ++c) {
                        if(transport[w][c].get(GRB.DoubleAttr.X) > 0.0001){
                            System.out.println("    Transport " + transport[w][c].get(GRB.DoubleAttr.X) + " units to customer " + c);
                        }
                    }
                } else {
                    System.out.println("Warehouse " + w + " closed!");
                }
            }

            // Dispose o f model and environment
            model.dispose();
            env.dispose();
        } catch (GRBException e ) {
            System.out.println("Errorcode: " + e.getErrorCode() + ". " +
                    e.getMessage());
        }
    }
}