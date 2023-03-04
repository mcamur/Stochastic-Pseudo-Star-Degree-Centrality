import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

public class FullModel {
	public static final double MSTOSEC = 0.001;
	public static final double epsilon =  0.99;
	public static final double FUZZ = 1e-6;
	public static final boolean warmStart = false; 
	public static HashMap<Integer, HashSet<Integer>> adjacencyList;
	public static HashMap<Edge, Double> edgeWeights;
	public static int numNodes;
	public static int numEdges;
	public static float density;
	private static IloCplex masterCplex;
	public static void main(String[] args) throws IOException, IloException {


		adjacencyList = new HashMap<Integer, HashSet<Integer>>();	
	    edgeWeights = new HashMap<Edge, Double>();			
	    
	
		int numNodes =   ...;		
		int parameter = ...;
		double  p = ...;
		

		String DATADIR = ...;
		System.out.println(ConsoleColors.GREEN +   "--------------------------------------------------------------------------------------------------------------"  +ConsoleColors.RESET); 
		System.out.println(numNodes+ "_"+ parameter+"_" +p+"_" +l);
		System.out.println(ConsoleColors.GREEN +   "--------------------------------------------------------------------------------------------------------------"  +ConsoleColors.RESET); 

		
		HashMap<Integer, Double> upperBounds = new HashMap<>();
		Path myFilePath = Paths.get(DATADIR).resolve(Paths.get("upperBound.txt"));    //upper bound for the open neighborhood   	
		try {
			Files.lines(myFilePath).forEach(line -> {
				String[] columnValues = line.split("\\s+");
				upperBounds.put(Integer.parseInt(columnValues[0]), Double.parseDouble(columnValues[1]));
			});
		} catch (IOException e) {
			e.printStackTrace();
		}

		Scanner scanner = new Scanner(new FileReader(DATADIR+"adjacencyList.txt"));
		
		
		int edgeNum=0;
		while (scanner.hasNextLine()) {
			final String[] line = scanner.nextLine().split("\\s+");
			edgeNum++;
			Edge myEdge = new Edge(Integer.valueOf(line[0]), Integer.valueOf(line[1]));
			edgeWeights.put(myEdge, Double.valueOf(line[2]));
			adjacencyList.computeIfAbsent(Integer.valueOf(line[0]),
					k -> new HashSet<>()).add(Integer.valueOf(line[1]));
		}
		scanner.close();

	  
		numEdges = edgeNum/2;
		density = (float) (2*numEdges)/ (float)(numNodes*(numNodes-1));
		numNodes = adjacencyList.size();
	  
	   
		masterCplex = new IloCplex();
		masterCplex.setParam(IloCplex.Param.TimeLimit, 4*1800);
		masterCplex.setParam(IloCplex.Param.Threads, masterCplex.getNumCores());
		masterCplex.setOut(null);
		IloIntVar[] x = new IloIntVar[numNodes]; //center node variables
		IloIntVar[] y = new IloIntVar[numNodes]; //leaf variables	
		HashMap<Edge,IloNumVar> z = new HashMap<Edge,IloNumVar>(); //assignment variables
		HashMap<Edge,IloNumVar> m = new HashMap<Edge,IloNumVar>(); 
		HashMap<Edge,IloNumVar> n = new HashMap<Edge,IloNumVar>(); 
		
		IloLinearNumExpr obJexpr = masterCplex.linearNumExpr(); 
		IloLinearNumExpr expr = masterCplex.linearNumExpr(); 
		
		for(int i=0; i< numNodes; i++) { //Define the variables
			x[i] = masterCplex.boolVar("x_"+ i);
			masterCplex.add(x[i]);
			
			y[i] = masterCplex.boolVar("y_"+ i);
			masterCplex.add(y[i]);
			
			for(int j : adjacencyList.get(i)) {
					IloNumVar zVar =  masterCplex.numVar(0.0, Double.MAX_VALUE,"z_"+ i+ "_"+ j);  
					masterCplex.add(zVar);
					
					Edge myEdge = new Edge(i,j);
					z.put(myEdge, zVar);
					
					double prob = edgeWeights.get(myEdge);
					obJexpr.addTerm(prob, zVar);
			}
			
		}
		
		masterCplex.addMaximize(obJexpr, "Objective_Function");

		obJexpr.clear();
		expr.clear();
	
		for(int i=0; i< numNodes; i++) {
			for(int j : adjacencyList.get(i)) {
				Edge myEdge= new Edge(i,j);
				double edgeProb= edgeWeights.get(myEdge);
				if(j > i) {
					IloNumVar mVar = masterCplex.numVar(0.0, Double.MAX_VALUE,"m_"+ i+ "_"+ j); 
					masterCplex.add(mVar);   
					m.put(myEdge, mVar);
					expr.addTerm(mVar, Math.log(1-edgeProb));
					masterCplex.addGe(mVar, masterCplex.sum(-1,masterCplex.sum(y[i],y[j])) , "mConst_"+i+ "_"+j);

				}

				IloNumVar nVar = masterCplex.numVar(0.0, Double.MAX_VALUE,"n_"+ i+ "_"+ j); 
				masterCplex.add(nVar);
				n.put(myEdge, nVar);
				
				expr.addTerm(nVar, Math.log(edgeProb));
				masterCplex.addGe(nVar, masterCplex.sum(-1,masterCplex.sum(x[i],y[j])) , "nConst_"+i+"_"+j);
			}
		}

		
		double RHS= Math.log(1 - epsilon);
		
		IloRange knapsackConst = masterCplex.addGe(expr, RHS, "Constraint_F");
		masterCplex.add(knapsackConst);
		expr.clear();
		
		for(int i=0; i< numNodes; i++) { 
			for(int j : adjacencyList.get(i)) {
				masterCplex.addLe(z.get(new Edge(i,j)),masterCplex.sum(x[i], y[i]), "C_"+i);
			}
		}
	
		for(int i=0; i< numNodes; i++) { 
			expr.clear();
			expr.addTerm(1, x[i]);
			expr.addTerm(1, y[i]);
			for(int j :adjacencyList.get(i)) {
				expr.addTerm(1, z.get(new Edge(j, i)));
			}
			masterCplex.addLe(expr, 1, "B_"+i);
		}
		
		
		for (int i = 0; i < numNodes; i++) { //to be a leaf node, you should be adjacent to the center
			expr.clear();	
			if(adjacencyList.get(i) != null) {
				adjacencyList.get(i).forEach((val) -> { 
					try { 
						expr.addTerm(x[val], 1);
					}catch (IloException e) {
						e.printStackTrace();
					} 
				});
			}
			masterCplex.addLe(y[i],expr, "D_"+i);		                	
		}
		expr.clear();
		
		for(int i=0; i< numNodes; i++)  //only one node can be center a.k.a a single star is created
			expr.addTerm(x[i], 1);  
		 
		masterCplex.addEq(expr, 1, "E");
		expr.clear(); 
		
		//To warm start, the solution obtained via the greedy heuristic is needed
		if(warmStart) {		        	 
			TreeMap<Integer, Set<Integer>> stars = new TreeMap<Integer, Set<Integer>>();  
			 myFilePath = Paths.get(DATADIR).resolve(Paths.get("star.txt"));      	
			 try {
				 Files.lines(myFilePath).forEach(line -> {
					 String[] columnValues = line.split("\\s+");
					 Set<Integer> leaves = new HashSet<Integer>();
					 for(int i =1 ; i < columnValues.length; i++) {
						leaves.add(Integer.valueOf(columnValues[i]));
					 }
					 
					 TreeSet<Integer> treeSet = new TreeSet<Integer>(leaves); 
					 leaves.clear();    
					 stars.put(Integer.valueOf(columnValues[0]),treeSet );
					 treeSet.clear();
				   
				 });
				 
			 } catch (IOException e) {
				 e.printStackTrace();
			 }
		 
			 for (Map.Entry<Integer,Set<Integer>> e : stars.entrySet()) {
				  IloNumVar[] vars = new IloNumVar[1 + e.getValue().size()];
				  double[] vals = new double[vars.length];
				  int idx = 0;
				  vars[idx] = x[e.getKey()];
				  vals[idx++] = 1;
				  for (Integer yIdx : e.getValue()) {
					  vars[idx] = y[yIdx];
					  vals[idx++] = 1;
				  }
				 
				  masterCplex.addMIPStart(vars, vals);
				  vars =null;
				  vals= null;
				}
		   stars.clear();
			 
		  }
//Solve the model via automated Benders provided by CPLEX with the following three lines
//			        masterCplex.setParam(IloCplex.Param.Benders.Strategy,
//	                         IloCplex.BendersStrategy.Full);
//					masterCplex.exportModel("IPmodel.lp");
		//System.out.println("number of variables : " + masterCplex.getNcols());
		
		long start = System.currentTimeMillis();
	   
		if (masterCplex.solve()) {
			// post process the results
			double time_taken = (System.currentTimeMillis() - start) * MSTOSEC;
			masterCplex.end();
			adjacencyList.clear();	
			edgeWeights.clear();	

		   
		}else {
			System.out.println("Optimal solution is not obtained. Check the instance!");
			masterCplex.end();
		}
			   

	}

}


