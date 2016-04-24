package cz.sange;

import gurobi.*;
import java.io.*;

/**
 * INPUT format
 3 8                                                        // numWarehouses numCustomers
 5000000 4000000 4500000                                    // fixed costs for warehouses
 1000000 800000 1250000                                    // sources for warehouses
 200000 200000 200000 200000 250000 250000 250000 250000    // demands for customers
 4 5 5 4 4 4.2 3.3 5                                        // transport costs from warehouses to customers
 2.5 3.5 4.5 3 2.2 4 2.6 5
 2 4 5 2.5 2.6 3.8 2.9 5.5
 *
 */

public class Example {

    private int [] demands;
    private int [] sources;
//    private int [] fixedCosts;
    private double[][] transportCosts;
    private double[][] fixedTransportCosts;

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

//            // fixedCosts costs by warehouse
//            fixedCosts = new int[w];
//            for (int i = 0; i < w; i++) {
//                fixedCosts[i] = (int) st.nval;
//                st.nextToken();
//            }

            // sources by warehouse
            sources = new int[w];
            for (int i = 0; i < w; i++) {
                sources[i] = (int) st.nval;
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

            // fixed transport costs
            fixedTransportCosts = new double[w][c];
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < c; j++) {
                    fixedTransportCosts[i][j] = st.nval;
                    st.nextToken();
                }
            }
        }

    }

    private void solve() {
        try {

            int numWarehouses = sources.length;
            int numCustomers = demands.length;

            //// model
            GRBEnv env = new GRBEnv("example.log") ;
            GRBModel model = new GRBModel(env) ;
            model.set(GRB.StringAttr.ModelName, "FCTP");

            GRBVar [][] open = new GRBVar[numWarehouses][numCustomers];
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    open[w][c] =
                            model.addVar(0, 1, fixedTransportCosts[w][c], GRB.BINARY, "Open " + w + "->" + c);
                }
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

            /// set objective
//            GRBLinExpr expr = new GRBLinExpr();
//            for (int w = 0; w < numWarehouses; w++) {
//                for (int c = 0; c < numCustomers; c++) {
//                    expr.addTerm(fixedTransportCosts[w][c], open[w][c]);
//                    expr.addTerm(transportCosts[w][c], transport[w][c]);
//                }
//            }
//            model.setObjective(expr, GRB.MINIMIZE);

            //// constraints
            // sources constraints
            for (int w = 0; w < numWarehouses; ++w) {
                GRBLinExpr sourcesSum = new GRBLinExpr();
                for (int c = 0; c < numCustomers; ++c) {
                    sourcesSum.addTerm(1.0, transport[w][c]);
                }
                model.addConstr(sourcesSum, GRB.EQUAL, sources[w], "Capacity " + w);
            }
            // demand constraints
            for (int c = 0; c < numCustomers; ++c) {
                GRBLinExpr demandsSum = new GRBLinExpr();
                for (int w = 0; w < numWarehouses; ++w) {
                    demandsSum.addTerm(1.0, transport[w][c]);
                }
                model.addConstr(demandsSum, GRB.EQUAL, demands[c], "Demand" + c);
            }
            // M constraints
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    double M = Math.min(sources[w], demands[c]);
                    GRBLinExpr expr = new GRBLinExpr();
                    expr.addTerm(M, open[w][c]);
                    model.addConstr(transport[w][c], GRB.LESS_EQUAL, expr, "mins...");
                }
            }

            //// prep
//             open all paths
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    open[w][c].set(GRB.DoubleAttr.Start, 1.0);
                }
            }
            // find max cfixed cost
            System.out.println("Initial guess:");
            double maxFixed = -GRB.INFINITY;
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    if (fixedTransportCosts[w][c] > maxFixed) {
                        maxFixed = fixedTransportCosts[w][c];
                    }
                }
            }
            // close most expensive
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    if (fixedTransportCosts[w][c] == maxFixed) {
                        open[w][c].set(GRB.DoubleAttr.Start, 0.0);
                        System.out.println("Closing w: " + w + "\n");
                        break;
                    }
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
                for (int c = 0; c < numCustomers; ++c) {
                    // if is open
//                if(open[w].get(GRB.DoubleAttr.X) == 1.0) {
                    if (open[w][c].get(GRB.DoubleAttr.X) == 1.0) {
                        System.out.println("Path (" + w + "," + c + ") open, transport:");
//                    for (int c = 0; c < numCustomers; ++c) {
                        if (transport[w][c].get(GRB.DoubleAttr.X) > 0.0001) {
                            System.out.println("\t" + Math.round(transport[w][c].get(GRB.DoubleAttr.X))
                                    + " units for VP: " + transportCosts[w][c]
                                    + " and FP: " + fixedTransportCosts[w][c] + " to customer " + c);
                        }
//                    }
                    } else {
                        System.out.println("Path (" + w + "," + c + ") closed.");
                    }
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