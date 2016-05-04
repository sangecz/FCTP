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

public class Benders {

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
    private GRBEnv envRMP;
    private GRBModel modelRMP;
    private GRBVar[][] yRMP;
    private boolean firstRunRMP;
    private GRBLinExpr objRMP;


    public static void main(String[] args) {

        try {
            System.setIn(new FileInputStream(args[0]));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        Benders benders = new Benders();

        try {
            benders.read();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // check supply-demand balance
        benders.checkSuppliesDemandsBalance();

        // prepare big-M upperbounds on x
        benders.prepareUpperbounds();

        try {
//            benders.displaySolution(benders.standardMIP(null));

            benders.bendersAlg();
        } catch (GRBException e) {
            e.printStackTrace();
        }

    }

    private void bendersAlg() throws GRBException{

        int[][] y = getInitialY();

        // sumy u, v pro bounded
        double cutConst = 0;
        // koeficienty sumy -xup*w*y
        double[][] cutCoefCommon = null;
        // pro bounded
        double[][] cutCoef = new double[numCustomers][numWarehouses];

        boolean converged = false;
        int iter = 0;
        int [][] ybest = new int[numWarehouses][numCustomers];

        double bound;
        // LB := -inf
        double LB = Double.MIN_VALUE;
        // UB := +inf
        double UB = Double.MAX_VALUE;

        boolean bounded = false;

        bendersRMPInit();

        while (!converged) {
            // {benders subproblem}
            Solution subProbSol = bendersDSP(y);

            if (subProbSol.getStatus() == GRB.Status.UNBOUNDED) {
                bounded = false;
                // add extreme rays
            } else if (subProbSol.getStatus() == GRB.Status.OPTIMAL){
                bounded = true;
                double sum = 0;
                for (int i = 0; i < numWarehouses; i++) {
                    for (int j = 0; j < numCustomers; j++) {
                        sum += fixedTransportCosts[i][j] * y[i][j];
                    }
                }
                // updates UB
                bound = sum + subProbSol.getObjVal();
                if(bound < UB) {
                    UB = bound;
                    for (int i = 0; i < numWarehouses; i++) {
                        System.arraycopy(y[i], 0, ybest[i], 0, numCustomers);
                    }
                    System.out.println("ybest#"+iter+ " arr: " + Arrays.deepToString(ybest));
               }
            } else {
                System.err.println("!UNBOUNDED and !OPTIMAL");
            }

            // cut data
            double suppliesSum = 0;
            for (int i = 0; i < numWarehouses; i++) {
                suppliesSum += -supplies[i] * subProbSol.getU()[i];
            }
            double demandsSum = 0;
            for (int j = 0; j < numCustomers; j++) {
                demandsSum += demands[j] * subProbSol.getV()[j];
            }
            cutConst = suppliesSum + demandsSum;

            cutCoef = new double[numWarehouses][numCustomers];
            for (int i = 0; i < numWarehouses; i++) {
                for (int j = 0; j < numCustomers; j++) {
                    cutCoef[i][j] = -xup[i][j] * subProbSol.getW()[i][j];
                }
            }

            // {benders master problem}
            Solution masProbSol = bendersRMP(bounded, cutConst,  cutCoef);
            // fixing y for DSP
            y = masProbSol.getY();

            if(masProbSol.getStatus() == GRB.Status.INFEASIBLE){
                System.err.println("RMP is infeasible");
                System.exit(3);
            }
            if(masProbSol.getStatus() != GRB.Status.OPTIMAL){
                System.err.println("RMP not solved to optimality");
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

        bendersRMPDestroy();

        // recover solution
        Solution fctpSol = standardMIP(ybest);
        if (fctpSol.getStatus() != GRB.Status.OPTIMAL) {
            System.err.println("Final LP not solved to optimality");
            System.exit(6);
        }

        displaySolution(fctpSol);
    }

    /**
     * Open all paths - feasible
     * @return
     * @throws GRBException
     */
    private int[][] getInitialY() throws GRBException {
        int [][] y = new int[numWarehouses][numCustomers]; // init zeros
        for (int i = 0; i < numWarehouses; i++) {
            for (int j = 0; j < numCustomers; j++) {
                y[i][j] = 1;
            }
        }
        return y;
    }

    private void displaySolution(Solution fctpSol)  {
        // print result
        System.out.println("\nTOTAL COSTS: " + Math.round(fctpSol.getObjVal()));
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

    private void bendersRMPInit()  throws GRBException {
        //// model
        envRMP = new GRBEnv("bendersRMP.log");
        modelRMP = new GRBModel(envRMP);
        // method: https://www.gurobi.com/documentation/6.5/refman/method.html
        modelRMP.getEnv().set(GRB.IntParam.Method, 1);
        modelRMP.set(GRB.StringAttr.ModelName, "FCTP");

        yRMP = new GRBVar[numWarehouses][numCustomers];
        for (int i = 0; i < numWarehouses; ++i) {
            for (int j = 0; j < numCustomers; ++j) {
                yRMP[i][j] = modelRMP.addVar(0, 1, 0, GRB.BINARY, "Open " + i + "->" + j);
            }
        }

        firstRunRMP = true;

        objRMP = new GRBLinExpr();
    }

    private void bendersRMPDestroy()  throws GRBException {
        // Dispose o f model and environment
        modelRMP.dispose();
        envRMP.dispose();
    }

    /**
     * Benders relaxed master problem
     * @throws GRBException
     */
    private Solution bendersRMP(boolean boundedDSP, double cutConst, double [][] cutCoef) throws GRBException {

        if(!firstRunRMP) {
            // vycisti promenne modelu
           modelRMP.reset();
        } else {
            firstRunRMP = false;
        }

        if(!boundedDSP) {
            modelRMP.set(GRB.IntAttr.ModelSense, 1);
        }

        //// update vars
        modelRMP.update();

        //// constraints
        // spolecne pro bounded i unbounded:
        GRBLinExpr cutCoeffSum = new GRBLinExpr();
        for (int i = 0; i < numWarehouses; i++) {
            for (int j = 0; j < numCustomers; j++) {
                cutCoeffSum.addTerm(cutCoef[i][j], yRMP[i][j]);
            }
        }
        GRBLinExpr cut = new GRBLinExpr();
        cut.addConstant(cutConst);
        cut.add(cutCoeffSum);

        if(boundedDSP) {
            GRBLinExpr fixedSum = new GRBLinExpr();
            for (int i = 0; i < numWarehouses; i++) {
                for (int j = 0; j < numCustomers; j++) {
                    fixedSum.addTerm(fixedTransportCosts[i][j], yRMP[i][j]);
                }
            }
            cut.add(fixedSum);
            // cut: (sum fixed*y) + (sum -supplies*u) + (sum demands*v) + (sum -xup*w*y) <= z
//            model.addConstr(z0, GRB.GREATER_EQUAL, cut, "RMP cut");
            objRMP.add(cut);
            // RMP obj
            modelRMP.setObjective(objRMP, GRB.MINIMIZE);
        } else {
            // unbounded
            // cut:  (sum -supplies*u) + (sum demands*v) + (sum -xup*w*y) <= 0
            modelRMP.addConstr(cut, GRB.LESS_EQUAL, 0, "RMP unboundedcut");
        }

        GRBLinExpr ycond1 = new GRBLinExpr();
        for (int j = 0; j < numCustomers; j++) {
            for (int i = 0; i < numWarehouses; i++) {
                ycond1.addTerm(supplies[i], yRMP[i][j]);
            }
            modelRMP.addConstr(ycond1, GRB.GREATER_EQUAL, demands[j], "-");
        }
        GRBLinExpr ycond2 = new GRBLinExpr();
        for (int i = 0; i < numWarehouses; i++) {
            for (int j = 0; j < numCustomers; j++) {
                ycond2.addTerm(demands[j], yRMP[i][j]);
            }
            modelRMP.addConstr(ycond2, GRB.GREATER_EQUAL, supplies[i], "-");
        }

        // solve
        modelRMP.optimize();

        int status = modelRMP.get(GRB.IntAttr.Status);
        double objVal = 0;
        int[][] solY = null;
        if (status == GRB.Status.OPTIMAL) {
            objVal = modelRMP.get(GRB.DoubleAttr.ObjVal);
            solY = new int[numWarehouses][numCustomers];
            for (int i = 0; i < numWarehouses; i++) {
                for (int j = 0; j < numCustomers; j++) {
                    solY[i][j] = (int) yRMP[i][j].get(GRB.DoubleAttr.X);  // 0 or 1
                }
            }
        }

        // retrieve and store solution
        Solution sol = new Solution();
        sol.setStatus(status);
        if (status == GRB.Status.OPTIMAL) {
            sol.setObjVal(objVal);
            sol.setY(solY);
        }
        return sol;
    }

    /**
     * Benders dual subproblem
     * @param y fixnute y
     * @return reseni POJO
     * @throws GRBException
     */
    private Solution bendersDSP(int [][] y) throws GRBException {
        //// model
        GRBEnv env = new GRBEnv("bendersDSP.log") ;
        GRBModel model = new GRBModel(env) ;
        // To obtain extreme rays for unbounded models: https://www.gurobi.com/documentation/6.0/refman/infunbdinfo.html#parameter:InfUnbdInfo
        model.getEnv().set(GRB.IntParam.InfUnbdInfo, 1);
        // method: primal simplex: https://www.gurobi.com/documentation/6.5/refman/method.html
        model.getEnv().set(GRB.IntParam.Method, 0);
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
                GRBLinExpr subConstr = new GRBLinExpr();
                subConstr.addTerm(-1.0, u[i]);
                subConstr.addTerm(1.0, v[j]);
                subConstr.addTerm(-1.0, w[i][j]);
                // -u_i + v_j - w_ij <= c_ij
                model.addConstr(subConstr, GRB.LESS_EQUAL, transportCosts[i][j], "subconstr");
            }
        }

        // solve
        model.optimize();

        // retrieve and store solution
        int status = model.get(GRB.IntAttr.Status);
        boolean unbounded = status == GRB.Status.UNBOUNDED;
        double objVal = model.get(GRB.DoubleAttr.ObjVal);

        // gets (un)bounded values for u,v,w for use in RMP
        // https://www.gurobi.com/documentation/6.5/refman/unbdray.html#attr:UnbdRay
        double [][] wSol = new double[numWarehouses][numCustomers];
        double [] uSol = new double[numWarehouses];
        double [] vSol = new double[numCustomers];
        for (int i = 0; i < numWarehouses; i++) {
            uSol[i] = unbounded ? u[i].get(GRB.DoubleAttr.UnbdRay) : u[i].get(GRB.DoubleAttr.X);
            for (int j = 0; j < numCustomers; j++) {
                vSol[j] = unbounded ? v[j].get(GRB.DoubleAttr.UnbdRay) : v[j].get(GRB.DoubleAttr.X);
                wSol[i][j] = unbounded ? w[i][j].get(GRB.DoubleAttr.UnbdRay) : w[i][j].get(GRB.DoubleAttr.X);
            }
        }

        // Dispose o f model and environment
        model.dispose();
        env.dispose();

        Solution z = new Solution();
        z.setStatus(status);
        z.setObjVal(objVal);
        z.setW(wSol);
        z.setU(uSol);
        z.setV(vSol);
        return z;
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

    private Solution standardMIP(int[][] y_fx) throws GRBException {
        //// model
        GRBEnv env = new GRBEnv("standardMIP.log") ;
        GRBModel model = new GRBModel(env) ;
        model.getEnv().set(GRB.IntParam.Method, 1);
        model.set(GRB.StringAttr.ModelName, "FCTP");

        //// variables
        boolean fixedY = y_fx != null;
        GRBVar[][] y = null;
        y = new GRBVar[numWarehouses][numCustomers];
        for (int w = 0; w < numWarehouses; ++w) {
            for (int c = 0; c < numCustomers; ++c) {
                if (!fixedY) {
                    y[w][c] = model.addVar(0, 1, fixedTransportCosts[w][c], GRB.BINARY, "Open " + w + "->" + c);
                } else {
                    y[w][c] = model.addVar(y_fx[w][c], 1, fixedTransportCosts[w][c], GRB.BINARY, "Open " + w + "->" + c);
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
        //// prep -- initial feasible integer solution: open all paths
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

        // retrieve and store solution
        int[][] ySol = new int[numWarehouses][numCustomers];
        double[][] xSol = new double[numWarehouses][numCustomers];
        for (int i = 0; i < numWarehouses; i++) {
            for (int j = 0; j < numCustomers; j++) {
                if (!fixedY) {
                    ySol[i][j] = (int) y[i][j].get(GRB.DoubleAttr.X);
                    xSol[i][j] =  x[i][j].get(GRB.DoubleAttr.X);
                } else {
                    ySol[i][j] = y_fx[i][j];
                    xSol[i][j] =  x[i][j].get(GRB.DoubleAttr.X);
                }
            }
        }

        Solution sol = new Solution();
        sol.setY(ySol);
        sol.setX(xSol);
        sol.setObjVal(model.get(GRB.DoubleAttr.ObjVal));
        sol.setStatus(model.get(GRB.IntAttr.Status));

        // Dispose o f model and environment
        model.dispose();
        env.dispose();
        return sol;
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