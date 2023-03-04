import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
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

public class ClassicalBendersDecomposition {
	public static final double FUZZ = 1e-6;
	public static final double MSTOSEC = 0.001;
	public static final double epsilon =  0.99;
	public static final boolean warmStart = true; 
	public static HashMap<Integer, HashSet<Integer>> adjacencyList;
	public static HashMap<Edge, Double> edgeWeights;
	public static HashMap<Integer, TreeSet<NodeWeight>> sortedWeights;
	public static int numNodes;
	public static int numEdges;
	public static float density;
	public static HashMap<Integer, Double> maxProbabilities;
	public static void main(String[] args) throws IOException, IloException {
		// TODO Auto-generated method stub
	
		adjacencyList = new HashMap<Integer, HashSet<Integer>>();	
	    edgeWeights = new HashMap<Edge, Double>();		
	    maxProbabilities = new HashMap<Integer, Double>(); 
	    sortedWeights  = new HashMap<Integer, TreeSet<NodeWeight>>();
	
		int numNodes =   ...;		
		int parameter = ...;
		double  p = ...;
						
		String DATADIR = "....";
			
						
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
		//create the adjacency list	
		while (scanner.hasNextLine()) {
			final String[] line = scanner.nextLine().split("\\s+"); 
			edgeNum++;
			Edge myEdge = new Edge(Integer.valueOf(line[0]), Integer.valueOf(line[1]));
			sortedWeights.computeIfAbsent(Integer.valueOf(line[0]),
					k -> new TreeSet<NodeWeight>()).add(new NodeWeight(Integer.valueOf(line[1]), Float.valueOf(line[2])));
			edgeWeights.put(myEdge, Double.valueOf(line[2]));
			adjacencyList.computeIfAbsent(Integer.valueOf(line[0]),
					k -> new HashSet<>()).add(Integer.valueOf(line[1]));
		}
	  
	   
		numEdges = edgeNum/2;
		density = (float) (2*numEdges)/ (float)(numNodes*(numNodes-1))   ;
		scanner.close();
		maxProbabilities = maxProb(); 
		
		numNodes = adjacencyList.size();
		
	   

		IloCplex masterCplex = new IloCplex();
		masterCplex.setParam(IloCplex.Param.TimeLimit, 4*1800);
		masterCplex.setParam(IloCplex.Param.Threads, masterCplex.getNumCores());

		
		IloIntVar[] x = new IloIntVar[numNodes]; //center node variables
		IloIntVar[] y = new IloIntVar[numNodes]; //leaf variables	
		IloNumVar[] estObj = new IloNumVar[numNodes];
	 
		HashMap<Edge,IloNumVar> m = new HashMap<Edge,IloNumVar>(); 
		HashMap<Edge,IloNumVar> n = new HashMap<Edge,IloNumVar>(); 
		
  
		IloLinearNumExpr expr = masterCplex.linearNumExpr(); 
		
	  //Define the variables: x,y, estimation of z
		for(int i=0; i< numNodes; i++) { 
			x[i] = masterCplex.boolVar("x_"+ i);
			masterCplex.add(x[i]);
			
			y[i] = masterCplex.boolVar("y_"+ i);
			masterCplex.add(y[i]);
			
			estObj[i]= masterCplex.numVar(0.0, maxProbabilities.get(i),"t_"+ i );  
			masterCplex.add(estObj[i]);
			
		}
		
		masterCplex.addMaximize(masterCplex.sum(estObj), "obj");//objective
		
		 for (int i = 0; i < numNodes; i++) { //upper bound for the size of open neighborhood
				expr.addTerm(x[i], upperBounds.get(i));
		  }
		masterCplex.addLe(masterCplex.sum(estObj) , expr, "upper_Bound");
		
	   
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
		
	   
	   
		
		double LB = 0;
		double UB = Double.MAX_VALUE;
		int iter = 0;
		
		double[] xSol = new double[numNodes];
		double[] ySol = new double[numNodes];
		double[] zMaster =new double[numNodes];
		
		double[] beta =new double[numNodes]; //dual variables beta and gamma
		double[] dualSol = new double[numNodes]; //objective of dual problem

		masterCplex.setOut(null);
		double  time_tracker =0;
		double remaning_time = Double.MAX_VALUE;;
		long start =0;
	  
		double optimality_gap = 1;
  
		while_loop:
		while ( optimality_gap >= FUZZ && time_tracker <= 4*1800 ) {
			
			start = System.currentTimeMillis();
			if (masterCplex.solve()) {	
				iter = iter + 1;
				UB = masterCplex.getObjValue();
				int numLeaves =0;
				int center = -1;
				
				
				Set<Integer> leaves = new HashSet<Integer>();
				for (int i = 0; i < numNodes; ++i) {	            	 
					xSol[i] = masterCplex.getValue(x[i]);
					ySol[i] = masterCplex.getValue(y[i]);
					
					if(ySol[i] > FUZZ) {
						numLeaves +=1;
						leaves.add(i);
					}
						
					if(xSol[i] > 0.5)
						center = i;
					
					zMaster[i] = masterCplex.getValue(estObj[i]);
				 }
				 TreeSet<Integer> treeSet = new TreeSet<Integer>(leaves);	
				 leaves.clear();
				 HashMap<Edge,Double> gamma = new HashMap<Edge,Double>(); 
				  
				 boolean tracker;
				 int counter;
				 int other_counter;
				 double obj1;    //two components of dual objective
				 double obj2; 
				  
				  for(int i =0; i< numNodes; i++) {
					  obj1 = 1- xSol[i]- ySol[i];
					  obj2 = 0;
					  beta[i] =0;
					  tracker =true;
					  counter =0;
					 
					  for(NodeWeight nw : sortedWeights.get(i)) {
						   int k = nw.getNode();
						   if(tracker) {
							   obj2 = obj2 +  xSol[k] + ySol[k];
							 
							   if(obj2 > obj1 + FUZZ) {		    			
								   beta[i] =  nw.getWeight();
								   tracker = false;
							   }else {
								   counter ++;
							   }
							   
						   }else {	   
					
							  gamma.put(new Edge(k, i), 0.0);
						   }
						}
						  other_counter=0;
						  myLoop:
						  for(NodeWeight nw : sortedWeights.get(i)) {
							  if(other_counter <= counter) {
								  int k = nw.getNode(); 
								  double val =   nw.getWeight() - beta[i];
								  gamma.put(new Edge(k, i),val);
								  other_counter ++;
							  }else {
								  break myLoop;
							  }
						  }
					  }
			  
				  for(int i =0; i< numNodes; i++) { //dual objective for each cut
					  dualSol[i] =   (1- xSol[i] - ySol[i])*beta[i];
					  double sum =0;
					  for(int k :adjacencyList.get(i)) {
						  sum  =  sum + (gamma.get(new Edge(k,i))*(xSol[k] + ySol[k]));
						}
						dualSol[i] = dualSol[i] +  sum;	
				  }
				   double tempLB = Arrays.stream(dualSol).sum();
				   expr.clear();
				   for (int i = 0; i < numNodes; i++) {	

							expr.setConstant(beta[i]);	
							expr.addTerm(-beta[i],x[i]);
							expr.addTerm(-beta[i],y[i]);
							for(int k : adjacencyList.get(i)){
								double val = gamma.get(new Edge(k,i));
								expr.addTerm(val, x[k]);
								expr.addTerm(val, y[k]);
							 }
							masterCplex.addLe(estObj[i],expr); //create the optimality cut for node i
							expr.clear();
						  
				  }
				   gamma.clear();
				   if( tempLB  >  LB + FUZZ) { //update the lower bound
					   LB = tempLB;
				   }
				   
				   time_tracker = time_tracker +  ((System.currentTimeMillis() - start) * MSTOSEC);
				   remaning_time = 4*1800 - time_tracker;
				   
				   if(remaning_time< 60) {
					   break while_loop;
				   }
				   
				   if(masterCplex.isMIP()){
						int nb = masterCplex.getNMIPStarts();
						masterCplex.deleteMIPStarts(0, nb);
					}
				   
				   IloNumVar[] vars = new IloNumVar[1 + numLeaves];
				   double[] vals = new double[vars.length];
				   int idx = 0;
				   vars[idx] = x[center];
				   vals[idx++] = 1;
				   
				   for (Integer yIdx : treeSet) {
					  vars[idx] = y[yIdx];
					  vals[idx++] = 1;
					}
					treeSet.clear();
					masterCplex.addMIPStart(vars, vals);
					vars =null;
					vals= null;
				   
				   masterCplex.setParam(IloCplex.Param.TimeLimit, remaning_time );
				   
				   
				   optimality_gap = (UB-LB)/(1e-10+LB);
				   
			}else {
				System.out.println("Master problem is not solved at iteration :"+iter+" . "
						+ "Check the instance!:"+ numNodes+ "--"+parameter+ "--"+p);
				masterCplex.end();
			}
		}
		//Post process the solutions
		masterCplex.end();
		adjacencyList.clear();	
		edgeWeights.clear();	
		maxProbabilities.clear();
		sortedWeights.clear();
		   


	}
	public static HashMap<Integer, Double> maxProb() {
		HashMap<Integer, Double> myMap = new HashMap<Integer, Double>();
		for(int i=0; i< numNodes ; i++) {
			double maxValue= 0;
			for(int j : adjacencyList.get(i)) {
				Edge myEdge = new Edge(j,i);
				double current = edgeWeights.get(myEdge);
				if(current >= maxValue ) {
					maxValue = current;
				}
			}
			myMap.put(i,maxValue);
		}
		
		return myMap;
	}
}