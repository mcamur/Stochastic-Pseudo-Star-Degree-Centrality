import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import java.util.Map.Entry;

import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

public class GreedyHeuristic {

	public static final double eps =0.99;
	public static HashMap<Integer, HashSet<Integer>> adjacencyList;
	public static HashMap<Edge, Float> edgeWeights;
	public static final double FUZZ = 1e-10;
	public static void main(String[] args) throws IOException {

			
							
		int numNodes =   ...;		
		int parameter = ...;
		double  p = ...;
		String DATADIR = "...",

		
		adjacencyList = new HashMap<Integer, HashSet<Integer>>();
		
		edgeWeights = new HashMap<Edge, Float>();
		Path filePath = null;
		try {
		  filePath = Paths.get(DATADIR).resolve(Paths.get("adjacencyList.txt")); 
		} catch (Exception e) {
			System.out.println("Something went wrong.");
		}
			   
		
		try {
			Files.lines(filePath).forEach(line -> {
				String[] columnValues = line.split("\\s+");
				
				Edge myEdge = new Edge(Integer.valueOf(columnValues[0]), Integer.valueOf(columnValues[1]));
					 edgeWeights.put(myEdge,   Float.valueOf(columnValues[2]));
				adjacencyList.computeIfAbsent(Integer.valueOf(columnValues[0]),
						k -> new HashSet<>()).add(Integer.valueOf(columnValues[1]));
			});
		} catch (IOException e) {
			e.printStackTrace();
		}
		Graph<Integer, DefaultEdge> myGraph =  new SimpleGraph<>(DefaultEdge.class);
		
		
		for(int i = 0 ; i <numNodes; i++) //add nodes
			myGraph.addVertex(i);
						
		for(int i = 0 ; i< numNodes; i++) {	//add edges
			if(adjacencyList.get(i) != null) {
				for(Integer j : adjacencyList.get(i)) {
					myGraph.addEdge(i, j);
				}
			}
		}
	
		HashMap<Integer, Set<Integer>> star = new HashMap<Integer,Set<Integer>>();
		HashMap<Integer, HashMap<Integer, Integer>> assignmentStar = new HashMap<Integer, HashMap<Integer, Integer>>(); 
		
		for(int z= 0 ;  z < numNodes ; z++) {
			
			int center = z;
			
			Set<Integer> starElements = new HashSet<Integer>(); //leaves and center
			starElements.add(center);
		
			Set<Integer> candidateLeaves = Graphs.neighborSetOf(myGraph, center); //get the open neighborhood as candidate leaves
			
			HashMap<Integer, Integer> assignment = new HashMap<Integer, Integer>(); //store neighborhood assignments
			
			
			for(int j : adjacencyList.get(center)) { //assign each node in the open neighborhood
				assignment.put(j, center); //center covers each node initially
			}
			
			double totalProb= 1;
			boolean flag = true;	
			
			outerloop:
			while(flag) {
				
				if(candidateLeaves.isEmpty()) //if no more leaf, then break
					break;
				
				HashMap<Integer, Double> weight = calculateWeight(center, starElements, candidateLeaves, assignment); //calculate
				//weight of each leaf
		
				Iterator<Entry<Integer, Double>> iter = weight.entrySet().iterator();
				while (iter.hasNext()) {	//if weight is zero, then remove it from the neighborhood, no need to evaluate
					 Map.Entry<Integer,Double> entry  = iter.next();
					if(entry .getValue() <= 0) {
						candidateLeaves.remove(entry.getKey()); // it is already assigned
						iter.remove();
					}
				}
				boolean secondFlag =true;	
				
				if(!candidateLeaves.isEmpty()) {
				
					while(secondFlag) {
					
		
						double results[] = contribution(center,starElements , weight, candidateLeaves, totalProb);					
						
						int candidateLeaf = (int) results[0];
						
						double candidateTotalProb = results[1] ;
					
						if(candidateTotalProb < 1 -eps) {										
							candidateLeaves.remove(candidateLeaf); //infeasible! that is why try the next node in the open neighborhood
							weight.remove(candidateLeaf);
							if(candidateLeaves.isEmpty()) { //if no more leave, then break
								break outerloop; 
							}else {
								continue;
							}
				
						}else {
					
							totalProb = candidateTotalProb;
							
							assignment = updateAssignment(assignment, starElements,candidateLeaf );
							starElements.add(candidateLeaf); 		// add leaf into star
							candidateLeaves.remove(candidateLeaf); // shrink the set of candidate leaves
							
							assignment.remove(candidateLeaf); 	  // remove the assignment  
						
							secondFlag =false;
					
						}
					}
				
				}else {
					break outerloop; 
				}
			}
		
			star.put(center, starElements);
			assignmentStar.put(center, assignment);
		}
		//System.out.println(star);
	
		//verifyStar(star);

		int bestCenter = bestObjective(assignmentStar);
		Set<Integer> Leaves = star.get(bestCenter);
		Leaves.remove(bestCenter);
						
		}
	
	}
	private static int bestObjective(HashMap<Integer, HashMap<Integer, Integer>> assignmentStar) {
		double bestObjective = 0;
		int bestCenter = -1;
		for (Map.Entry<Integer, HashMap<Integer, Integer>> entry : assignmentStar.entrySet()) {
			Integer center = entry.getKey();
			double temp  =0;
			
			HashMap<Integer, Integer> myMap = entry.getValue();
			for (Map.Entry<Integer, Integer> it : myMap.entrySet()) {
				temp = temp + edgeWeights.get(new Edge(it.getKey(), it.getValue()));
			}

			if(bestObjective < temp) {
				bestObjective = temp;
				bestCenter = center;
			}
		}
		
		return bestCenter;
	}
	private static void verifyStar(HashMap<Integer, Set<Integer>> stars, double eps) {
		
	    Iterator it = stars.entrySet().iterator();
	    while (it.hasNext()) {
	        Map.Entry pair = (Map.Entry)it.next();
	        int center= (int) pair.getKey();
	        Set<Integer> leaves = (Set<Integer>) pair.getValue();
	        double p = 1;
	        //double centerProb =1;
	        for(int leaf : leaves ) {
	        	if(leaf!=center) {
	        		//centerProb = centerProb * edgeWeights.get(new Edge(center, leaf));
		        	p = p * edgeWeights.get(new Edge(center, leaf));
		        	for(int nextLeaf : leaves ) {
		        		if(nextLeaf> leaf & nextLeaf != center) {
		        			if(adjacencyList.get(leaf).contains(nextLeaf)) {
		        				p = p*(1-edgeWeights.get(new Edge(leaf,nextLeaf)));
		        			}
		        		}
		        	}
	        	}
	        }

	        if(p <  1 -eps) {
	        	System.out.println("trouble at center "+ center);
	        }else {
	        	System.out.println("center " + center+ " is safe!");
	        }
	        
	    }

	}
	
	private static HashMap<Integer, Integer> updateAssignment(HashMap<Integer, Integer> assignment,
			Set<Integer> starElements, int candidateLeaf) {
		
		for(int j : adjacencyList.get(candidateLeaf)) {
			if(assignment.containsKey(j)) {
				double prevProb = edgeWeights.get(new Edge(j, assignment.get(j)));
				double newProb =  edgeWeights.get(new Edge(j, candidateLeaf));
				if(prevProb < newProb) {
					assignment.put(j, candidateLeaf);
				}
			}else {
				if(!starElements.contains(j))
					assignment.put(j, candidateLeaf);
			}
		}
		return assignment;
	}
	
	private static double[] contribution(int center,  Set<Integer>  starElements, HashMap<Integer, Double> weight,
			Set<Integer> candidateLeaves, double totalProb) {
	
		double results[] = new double[2];
		double temp= 0;
		int candidateLeaf= -1;
		Set<Integer> newSet = new HashSet<Integer>();
		newSet.addAll(starElements);
		newSet.remove(center); //go with only leaf nodes
		
		if(newSet.isEmpty()) { //star only has the center. Initial step.
		
			for(int leaf: candidateLeaves) {
				double contr =  weight.get(leaf)*edgeWeights.get(new Edge(center, leaf));
				
				if(temp < contr){
					temp = contr;
					candidateLeaf=  leaf;
				}
			}
			totalProb =  totalProb * edgeWeights.get(new Edge(center, candidateLeaf));	
			
		}else{
			
			
			
			for(int leaf: candidateLeaves) {
			 double contr =  weight.get(leaf)*edgeWeights.get(new Edge(center, leaf));
			 
			 for(int realLeaf : newSet) {
				
		    		if(adjacencyList.get(leaf).contains(realLeaf)){
		    
		    			contr = contr*(1-edgeWeights.get(new Edge(leaf,realLeaf)));
		    			
		    		}
				}	

				if(temp < contr -FUZZ ){
					temp = contr;
					candidateLeaf=  leaf;
				}
			}
			
			totalProb =  totalProb * edgeWeights.get(new Edge(center, candidateLeaf));
		
			for(int realLeaf : newSet) {
	    		if(adjacencyList.get(candidateLeaf).contains(realLeaf)){
	    			totalProb = totalProb*(1-edgeWeights.get(new Edge(candidateLeaf,realLeaf)));
	    		}
	    		
			}
		}
	
		results[0] = candidateLeaf;
		results[1] = totalProb;

		return results;
	}
	
	private static HashMap<Integer, Double> calculateWeight(int center,  Set<Integer>  starElements,
			 Set<Integer> candidateLeaves, HashMap<Integer, Integer> assignment) {
		
		HashMap<Integer, Double> weight = new HashMap<Integer, Double>(); //weight of each node in the neighborhood
	
		for(int j : candidateLeaves) { //candidate leaf
			double w =0;
			for(int k : adjacencyList.get(j)) {
				if(!starElements.contains(k) ) {
					if(assignment.containsKey(k)) {
						int starElement  = assignment.get(k); 
						double currentProb = edgeWeights.get(new Edge(starElement,k));
						double newProb = edgeWeights.get(new Edge(j,k));
						if(newProb > currentProb) {
							w = w + newProb - currentProb;
						}
					}else {
						double addProb = edgeWeights.get(new Edge(j,k));
						w = w + addProb;
					}
				}
			}
			weight.put(j, w);
		}
		//System.out.println(weight);
		
		return weight;
	}
	


}
