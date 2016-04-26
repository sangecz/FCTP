package cz.sange;

import gurobi.*;
import java.io.*;
import java.util.Arrays;

/**
 * INPUT format
 4 3                                            // numWarehouses numCustomers
 10 30 40 20                                    // supplies for warehouses
 20 50 30                                       // demands for customers
 2.0 3.0 4.0                                    // transport costs matrix numWarehouses x numCustomers
 3.0 2.0 1.0
 1.0 4.0 3.0
 4.0 5.0 2.0
 10.0 30.0 20.0                                 // fixed costs matrix numWarehouses x numCustomers
 10.0 30.0 20.0
 10.0 30.0 20.0
 10.0 30.0 20.0
 *
 */

public class Example {

    private static int MAX_ITERATIONS = 50;
    private static double EPS = 0.0001;
    private static double CONVERGEND_BOUND = 0.1;
    public static final int UNBOUNDED = 1000000;

    private int [] demands;
    private int [] supplies;
    private double[][] transportCosts;
    private double[][] fixedTransportCosts;
    private int[][] xup;
    private int numWarehouses;
    private int numCustomers;


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

        // check supply-demand balance
        example.checkSuppliesDemandsBalance();

        // prepare bbig-M upperbounds on x
        example.prepareUpperbounds();

        try {
//            example.standardMIP();

            example.bendersAlg();
        } catch (GRBException e) {
            e.printStackTrace();
        }

    }

    private void bendersAlg() throws GRBException{

        Solution initialSol = standardMIP(null);
        int [][] y = initialSol.getY();

        int z0; // 'relaxed master objective variable
        boolean[] cutset = new boolean[MAX_ITERATIONS];
        boolean[] unbcutset = new boolean[MAX_ITERATIONS];

        double[] cutconst = new double[MAX_ITERATIONS];
        double[][][] cutcoef = new double[MAX_ITERATIONS][numWarehouses][numCustomers];

        boolean converged = false;
        int iter = 0;
        double bound;
        int [][] ybest = new int[numWarehouses][numCustomers];

        // LB := -inf
        double LB = Double.MIN_VALUE;
        // UB := +inf
        double UB = Double.MAX_VALUE;

        boolean unbounded = false;
        while ((UB - LB) > EPS) {
            // {benders subproblem}

            Solution subProbSol = bendersSubproblem(y);

            if (subProbSol.getStatus() == GRB.Status.UNBOUNDED) {
                unbounded = true;
                // Get unbounded ray u_
//                double u_ = subProbSol.getUnbdRay(); FIXME kde pouzit??

                Solution modifiedSubProbSol = bendersModifiedSubproblem(y);
                if(modifiedSubProbSol.getStatus() != GRB.Status.OPTIMAL){
                    System.err.println("Modified subproblem not solved to optimality");
                    System.exit(7);
                }
                unbcutset[iter] = true;
            } else if (subProbSol.getStatus() == GRB.Status.OPTIMAL){
                unbounded = false;
                double sum = 0;
                for (int i = 0; i < numWarehouses; i++) {
                    for (int j = 0; j < numCustomers; j++) {
                        sum += fixedTransportCosts[i][j] * y[i][j];
                    }
                }
                bound = sum + subProbSol.getObjVal();
                if(bound < UB) {
                    UB = bound;
                    for (int i = 0; i < numWarehouses; i++) {
                        System.arraycopy(y[i], 0, ybest[i], 0, numCustomers);
                    }
                    System.out.println("ybest: " + Arrays.deepToString(ybest));
                }
                cutset[iter] = true;
            }
            // cut data
            double suppliesSum = 0;
            for (int i = 0; i < numWarehouses; i++) {
                suppliesSum += -supplies[i] * subProbSol.getU()[i];
            }
            double demandsSum = 0;
            for (int j = 0; j < numWarehouses; j++) {
                demandsSum += demands[j] * subProbSol.getV()[j];
            }
            cutconst[iter] = suppliesSum + demandsSum;

            for (int i = 0; i < numWarehouses; i++) {
                for (int j = 0; j < numCustomers; j++) {
                    cutcoef[iter][i][j] = -xup[i][j] * subProbSol.getW()[i][j];
                }
            }
            // {benders master problem}
            // benders relaxed master problem:   min_y
            //                                    s.t.    z >= sum{(i,j) in WxC}f_{i,j}*u_{i,j} + sum(i in W)(-s_i)*u__i^(k) + sum{j in C}d_j*v__j^(k) + sum{(i,j) in WxC}(-M_{i,j}*w__{i,j}^(k))*y_{i,j}
            //                                    s.t.    sum(i in W)(-s_i)*u__i^(l) + sum(j in C)d_j*v__j^(l) + sum{(i,j) in WxC}(-M_{i,j}*w__{i,j}^(l))*y_{i,j} <= 0
            Solution masProbSol;
            if (!unbounded) {
                // FIXME cutset[iter] ? 1 : 0 / unbcutset[iter] ? 1 : 0  nedava smysl ??
                masProbSol = bendersMasterProblem(unbounded, cutconst[cutset[iter] ? 1 : 0], cutcoef[cutset[iter] ? 1 : 0]);
            } else {
                masProbSol = bendersMasterProblem(unbounded, cutconst[unbcutset[iter] ? 1 : 0], cutcoef[unbcutset[iter] ? 1 : 0]);
            }

            if(masProbSol.getStatus() == GRB.Status.INFEASIBLE){
                System.err.println("Relaxed Master is infeasible");
                System.exit(3);
            }
            if(masProbSol.getStatus() != GRB.Status.OPTIMAL){
                System.err.println("Masterproblem not solved to optimality");
                System.exit(4);
            }

            LB = masProbSol.getObjVal();
            System.out.println("iter#" + iter + ", LB=" + LB + ", UB=" + UB);
            converged = (UB - LB) < CONVERGEND_BOUND;

            if (iter++ >= MAX_ITERATIONS) {
                break;
            }
        }

        if (!converged) {
            System.err.println("No convergence");
            System.exit(5);
        }

        // recover solution
        for (int i = 0; i < numWarehouses; i++) {
            System.arraycopy(ybest[i], 0, y[i], 0, numCustomers);
        }
        Solution fctpSol = standardMIP(y);
        if (fctpSol.getStatus() != GRB.Status.OPTIMAL) {
            System.err.println("final lp not solved to optimality");
            System.exit(6);
        }

        displaySolution(fctpSol);
    }

    private Solution bendersModifiedSubproblem(int [][] y) throws GRBException {

        //// model
        GRBEnv env = new GRBEnv("bendersModifiedSubproblem.log") ;
        GRBModel model = new GRBModel(env) ;
        // To obtain extreme rays for unbounded models:
        model.getEnv().set(GRB.IntParam.InfUnbdInfo, 1);
        model.getEnv().set(GRB.IntParam.Method, 1);
        model.set(GRB.StringAttr.ModelName, "FCTP");

        //// variables
        GRBVar [] u = new GRBVar[numWarehouses];
        for (int i = 0; i < numWarehouses; ++i) {
            u[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "capacity constraint dual " + i);
        }

        GRBVar [] v = new GRBVar[numCustomers];
        for (int j = 0; j < numCustomers; ++j) {
            v[j] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "demand constraint dual " + j);
        }

        GRBVar [][] w = new GRBVar[numWarehouses][numCustomers];
        for (int i = 0; i < numWarehouses; ++i) {
            for (int j = 0; j < numCustomers; ++j) {
                w[i][j] =
                        model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "trans(" + i + ", " + j + ")");
            }
        }

        //// maximize
        model.set(GRB.IntAttr.ModelSense, -1);

        //// integrate variables
        model.update();

        //// constraints
        GRBLinExpr  supply = new GRBLinExpr();
        for (int i = 0; i < numWarehouses; ++i) {
            supply.addTerm(-supplies[i],  u[i]);  // sum{i in W}(-s_i)*u_i
        }

        GRBLinExpr  demand = new GRBLinExpr();
        for (int j = 0; j < numCustomers; ++j) {
            demand.addTerm(demands[j],  v[j]);  // sum{j in C}(d_j*v_j)
        }

        GRBLinExpr e = new GRBLinExpr();
        for (int i = 0; i < numWarehouses; ++i) {
            for (int j = 0; j < numCustomers; ++j) {
                e.addTerm(-xup[i][j] * y[i][j], w[i][j]);  // sum{(i,j) in WxC)(-xup_{i,j}*y_{i,j})*w_{i,j}
            }
        }

        GRBLinExpr modifiedsubobj = new GRBLinExpr();
        modifiedsubobj.add(supply);
        modifiedsubobj.add(demand);
        modifiedsubobj.add(e);
        model.addConstr(modifiedsubobj, GRB.EQUAL, 1, "modifiedsubobj");

        // modifiedsubconstr
        GRBLinExpr modifiedsubconstr = new GRBLinExpr();
        for (int i = 0; i < numWarehouses; ++i) {
            for (int j = 0; j < numCustomers; ++j) {
                modifiedsubconstr.addTerm(-1.0, u[i]);
                modifiedsubconstr.addTerm(1.0, v[j]);
                modifiedsubconstr.addTerm(-1.0, w[i][j]);
                model.addConstr(modifiedsubconstr, GRB.LESS_EQUAL, 0, "modifiedsubconstr");
            }
        }

        // solve
        model.optimize();

        int status = model.get(GRB.IntAttr.Status);

        // Dispose o f model and environment
        model.dispose();
        env.dispose();

        Solution z = new Solution();
        z.setStatus(status);
        return z;
    }

    private void displaySolution(Solution fctpSol) {
        // print result
        System.out.println("\nTOTAL COSTS: " + fctpSol.getObjVal());
        System.out.println("SOLUTION:");
        for (int w = 0; w < numWarehouses; ++w) {
            for (int c = 0; c < numCustomers; ++c) {
                if (fctpSol.getY()[w][c] == 1) {
                    System.out.println("Path (" + w + "," + c + ") open, transport:");
                    if (fctpSol.getX()[w][c] > EPS) {
                        System.out.println("\t" + Math.round(fctpSol.getX()[w][c])
                                + " units for VP: " + transportCosts[w][c]
                                + " and FP: " + fixedTransportCosts[w][c] + " to customer " + c);
                    }
                } else {
                    System.out.println("Path (" + w + "," + c + ") closed.");
                }
            }
        }
    }

    private Solution bendersMasterProblem(boolean unbounded, double cutconst,
                                          double [][] cutcoef) throws GRBException {
        //// model
        GRBEnv env = new GRBEnv("bendersMasterProblem.log") ;
        GRBModel model = new GRBModel(env) ;
        // To obtain extreme rays for unbounded models:
        model.getEnv().set(GRB.IntParam.InfUnbdInfo, 1);
        model.getEnv().set(GRB.IntParam.Method, 1);
        model.set(GRB.StringAttr.ModelName, "FCTP");

        GRBVar [][] y = new GRBVar[numWarehouses][numCustomers];
        for (int i = 0; i < numWarehouses; ++i) {
            for (int j = 0; j < numCustomers; ++j) {
                y[i][j] =
                        model.addVar(0, 1, 0, GRB.BINARY, "Open " + i + "->" + j);
            }
        }

        //// update vars
        model.update();

        if (!unbounded) {
            //// objective fnc
            /// cut(cutset)
            GRBLinExpr fixedSum = new GRBLinExpr();
            GRBLinExpr cutcoeffSum = new GRBLinExpr();
            for (int i = 0; i < numWarehouses; i++) {
                for (int j = 0; j < numCustomers; j++) {
                    fixedSum.addTerm(fixedTransportCosts[i][j], y[i][j]);
                    cutcoeffSum.addTerm(cutcoef[i][j], y[i][j]);
                }
            }
            GRBLinExpr cut = new GRBLinExpr();
            cut.add(fixedSum);
            cut.addConstant(cutconst);
            cut.add(cutcoeffSum);
            model.setObjective(cut, GRB.MINIMIZE);
        } else {
            model.set(GRB.IntAttr.ModelSense, 1);

            /// unboundedcut(unbcutset)
            GRBLinExpr cutcoeffUnbSum = new GRBLinExpr();
            for (int i = 0; i < numWarehouses; i++) {
                for (int j = 0; j < numCustomers; j++) {
                    cutcoeffUnbSum.addTerm(cutcoef[i][j], y[i][j]);
                }
            }
            GRBLinExpr unboundedcut = new GRBLinExpr();
            unboundedcut.addConstant(cutconst);
            unboundedcut.add(cutcoeffUnbSum);

            model.addConstr(unboundedcut, GRB.LESS_EQUAL, 0, "unboundedcut");
        }

        // solve
        model.optimize();

        int status = model.get(GRB.IntAttr.Status);
        double objVal = model.get(GRB.DoubleAttr.ObjVal);

        // Dispose o f model and environment
        model.dispose();
        env.dispose();

        Solution sol = new Solution();
        sol.setStatus(status);
        sol.setObjVal(objVal);
        return sol;
    }

    private Solution bendersSubproblem(int [][] y) throws GRBException {
        //// model
        GRBEnv env = new GRBEnv("bendersSubproblem.log") ;
        GRBModel model = new GRBModel(env) ;
        // To obtain extreme rays for unbounded models:
        model.getEnv().set(GRB.IntParam.InfUnbdInfo, 1);
        model.getEnv().set(GRB.IntParam.Method, 1);
        model.set(GRB.StringAttr.ModelName, "FCTP");

        //// variables
        GRBVar [] u = new GRBVar[numWarehouses];
        for (int i = 0; i < numWarehouses; ++i) {
                u[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "capacity constraint dual " + i);
        }

        GRBVar [] v = new GRBVar[numCustomers];
        for (int j = 0; j < numCustomers; ++j) {
            v[j] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "demand constraint dual " + j);
        }

        GRBVar [][] w = new GRBVar[numWarehouses][numCustomers];
        for (int i = 0; i < numWarehouses; ++i) {
            for (int j = 0; j < numCustomers; ++j) {
                w[i][j] =
                        model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "trans(" + i + ", " + j + ")");
            }
        }

        //// integrate variables
        model.update();

        //// equations
        // Benders' subproblem:   max_{u,v,w} sum{i in W}(-s_i)*u_i + sum{j in C}(d_j*v_j) + sum{(i,j) in WxC)(-xup_{i,j}*y_{i,j})*w_{i,j}
        //                          s.t.  -u_i + v_j - w_{i,j} <= c_{i,j}
        //                          s.t.  u_i >= 0,   v_j >= 0,   w_{i,j}>= 0

        //// objective function
        GRBLinExpr  supply = new GRBLinExpr();
        for (int i = 0; i < numWarehouses; ++i) {
            supply.addTerm(-supplies[i],  u[i]);  // sum{i in W}(-s_i)*u_i
        }

        GRBLinExpr  demand = new GRBLinExpr();
        for (int j = 0; j < numCustomers; ++j) {
            demand.addTerm(demands[j],  v[j]);  // sum{j in C}(d_j*v_j)
        }

        GRBLinExpr e = new GRBLinExpr();
        for (int i = 0; i < numWarehouses; ++i) {
            for (int j = 0; j < numCustomers; ++j) {
                e.addTerm(-xup[i][j] * y[i][j], w[i][j]);  // sum{(i,j) in WxC)(-xup_{i,j}*y_{i,j})*w_{i,j}
            }
        }

        // subobj
        GRBLinExpr subobj = new GRBLinExpr();
        subobj.add(supply);
        subobj.add(demand);
        subobj.add(e);
        model.setObjective(subobj, GRB.MAXIMIZE);

        //// constraints
        for (int i = 0; i < numWarehouses; ++i) {
            for (int j = 0; j < numCustomers; ++j) {
                GRBLinExpr subconstr = new GRBLinExpr();
                subconstr.addTerm(-1.0, u[i]);
                subconstr.addTerm(1.0, v[j]);
                subconstr.addTerm(-1.0, w[i][j]);
                model.addConstr(subconstr, GRB.LESS_EQUAL, transportCosts[i][j], "subconstr");
            }
        }

        // solve
        model.optimize();

        int status = model.get(GRB.IntAttr.Status);
        double objVal = model.get(GRB.DoubleAttr.ObjVal);
        double unbdRay = model.get(GRB.DoubleAttr.UnbdRay);
        double [][] wSol = new double[numWarehouses][numCustomers];
        double [] uSol = new double[numWarehouses];
        double [] vSol = new double[numCustomers];
        for (int i = 0; i < numWarehouses; i++) {
            uSol[i] = u[i].get(GRB.DoubleAttr.X);
            for (int j = 0; j < numCustomers; j++) {
                vSol[j] = v[j].get(GRB.DoubleAttr.X);
                wSol[i][j] = w[i][j].get(GRB.DoubleAttr.X);
            }
        }

        // Dispose o f model and environment
        model.dispose();
        env.dispose();

        Solution z = new Solution();
        z.setStatus(status);
        z.setObjVal(objVal);
        z.setUnbdRay(unbdRay);
        z.setW(wSol);
        z.setU(uSol);
        z.setV(vSol);
        return z;
    }

    private Solution standardMIP(int[][] y_fx) throws GRBException {
        //// model
        GRBEnv env = new GRBEnv("standardMIP.log") ;
        GRBModel model = new GRBModel(env) ;
        model.set(GRB.StringAttr.ModelName, "FCTP");

        //// variables
        boolean fixedY = y_fx != null;
        GRBVar[][] y = null;
        if (!fixedY) {
            y = new GRBVar[numWarehouses][numCustomers];
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    y[w][c] =
                            model.addVar(0, 1, fixedTransportCosts[w][c], GRB.BINARY, "Open " + w + "->" + c);
                }
            }
        }

        GRBVar [][] x = new GRBVar[numWarehouses][numCustomers];
        for (int w = 0; w < numWarehouses; ++w) {
            for (int c = 0; c < numCustomers; ++c) {
                x[w][c] =
                        model.addVar(0, GRB.INFINITY, transportCosts[w][c], GRB.CONTINUOUS, "trans(" + w + ", " + c + ")");
            }
        }

        //// minimize
        model.set(GRB.IntAttr.ModelSense, 1);

        //// update vars
        model.update();

        //// constraints
        // supplies constraints
        for (int w = 0; w < numWarehouses; ++w) {
            GRBLinExpr sourcesSum = new GRBLinExpr();
            for (int c = 0; c < numCustomers; ++c) {
                sourcesSum.addTerm(1.0, x[w][c]);
            }
            model.addConstr(sourcesSum, GRB.EQUAL, supplies[w], "Capacity " + w);
        }
        // demand constraints
        for (int c = 0; c < numCustomers; ++c) {
            GRBLinExpr demandsSum = new GRBLinExpr();
            for (int w = 0; w < numWarehouses; ++w) {
                demandsSum.addTerm(1.0, x[w][c]);
            }
            model.addConstr(demandsSum, GRB.EQUAL, demands[c], "Demand" + c);
        }
        // M constraints
        for (int w = 0; w < numWarehouses; ++w) {
            for (int c = 0; c < numCustomers; ++c) {
                double M = Math.min(supplies[w], demands[c]);
                GRBLinExpr expr = new GRBLinExpr();
                // upper bound pro transport[w][c] neboli x_{i,j}
                if (!fixedY) {
                    expr.addTerm(M, y[w][c]);
                } else {
                    expr.addConstant(M * y_fx[w][c]);
                }
                model.addConstr(x[w][c], GRB.LESS_EQUAL, expr, "mins...");
            }
        }

        // {initialization}
        //// prep -- initial feasible integer solution
//             open all paths
        if (!fixedY) {
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    y[w][c].set(GRB.DoubleAttr.Start, 1.0);
                }
            }
            // find max fixed cost
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
            branch:
            for (int w = 0; w < numWarehouses; ++w) {
                for (int c = 0; c < numCustomers; ++c) {
                    if (fixedTransportCosts[w][c] == maxFixed) {
                        y[w][c].set(GRB.DoubleAttr.Start, 0.0);
                        System.out.println("Closing path: (" + w + "," + c + ") of FP=" + maxFixed + "\n");
                        break branch;
                    }
                }
            }
        }

        //// solve
        model.optimize();

        Solution sol = new Solution();

        int[][] ySol = new int[numWarehouses][numCustomers];
        for (int i = 0; i < numWarehouses; i++) {
            for (int j = 0; j < numCustomers; j++) {
                if (!fixedY) {
                    ySol[i][j] = (int) y[i][j].get(GRB.DoubleAttr.X);
                } else {
                    ySol[i][j] = y_fx[i][j];
                }
            }
        }
        sol.setY(ySol);

        // Dispose o f model and environment
        model.dispose();
        env.dispose();
        return sol;
    }

    private void prepareUpperbounds() {
        xup = new int[numWarehouses][numCustomers];

        for (int w = 0; w < numWarehouses; w++) {
            for (int c = 0; c < numCustomers; c++) {
                xup[w][c] = Math.min(supplies[w], demands[c]);
            }
        }
    }

    private void checkSuppliesDemandsBalance() {

        int sumDemands = 0;
        for (int i = 0; i < demands.length; i++) {
            sumDemands += demands[i];
        }
        int sumSupplies = 0;
        for (int j = 0; j < supplies.length; j++) {
            sumSupplies += supplies[j];
        }

        if(sumDemands != sumSupplies) {
            System.err.println("Demands - Supplies not in balance");
            System.exit(1);
        }
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

            // supplies by warehouse
            supplies = new int[w];
            for (int i = 0; i < w; i++) {
                supplies[i] = (int) st.nval;
                st.nextToken();
            }

            // demands by customer
            demands = new int[c];
            for (int i = 0; i < c; i++) {
                demands[i] = (int) st.nval;
                st.nextToken();
            }

            // transport costs
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

        numWarehouses = supplies.length;
        numCustomers = demands.length;

    }
}