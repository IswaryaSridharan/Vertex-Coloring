import java.util.HashMap;
import java.util.InputMismatchException;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.Map;
import java.util.ArrayList;
import java.util.Random;
import java.util.*;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Comparator;

/*
* This file implements Dsatur's vertex coloring algorithm.
* We use three different data structures and measure the time taken and colors used.
* Technique 1: Dsatur algorithm that uses adjacency List to store all vertices and Linked List to store its adjacent vertices
* Technique 2: Dsatur algorithm with a max Heap(Priority Queue) to store max saturation degree of vertices
* Technique 3: Dsatur algorithm with a Red-Black tree(TreeSet) to store max saturation degree of vertices
*/

// A class to encapsulate the graph and common graph utility methods
class Graph
{
	Random random = new Random();

	// Use adjacency list for storing the graph
	List<Node> vertices = new ArrayList<>();
	List<LinkedList<Node>> adjList = new ArrayList<>();

	// To create an edge between vertex randomly
	public double randomNumberGenerator()
	{
		return Math.random();
	}


	// Adds an edge between vertices if the random number generator function returns a number less than the specified density
	public void addEdge(int count, double density)
	{
		for(int i=0; i<count; i++)
		{
			for(int j=i+1; j<count; j++)
			{
				double check = randomNumberGenerator();
				if(check < density)
				{
					adjList.get(i).add(vertices.get(j));
					adjList.get(j).add(vertices.get(i));

				}
			}
		}
	}


	// print the graph with vertices and edges connecting to its adjacent vertices
	public void printGraph(int count)
	{
		for(int i=0; i<count; i++)
		{
			System.out.print("Vertex : "+i +"  Edges :");
			for(int j=0; j<adjList.get(i).size(); j++)
			{
				System.out.print(" "+adjList.get(i).get(j).node);
			}
			System.out.println();

		}
	}

	// function to generate random graph with given vertices
	public void generateRandomGraph(int noOfVertex, double density)
	{
		int vertexCount = noOfVertex;
		for(int i=0; i<vertexCount; i++)
		{
			vertices.add(new Node(i, 1, 0, 0));
			adjList.add(new LinkedList<Node>());
		}
		addEdge(vertexCount, density);
		//printGraph(vertexCount);
	}

	// Populates the copy of original graph to be used in three different Dsatur implementation
	public void populateGraph(List<Node> verticesCopy, List<LinkedList<Node>> adjListCopy)
	{
		for(int i=0; i<vertices.size(); i++)
		{
			verticesCopy.add(new Node(vertices.get(i)));
		}
		for(int i=0; i<adjList.size(); i++)
		{
			adjListCopy.add(new LinkedList<Node>());
			for(int j=0; j<adjList.get(i).size(); j++)
			{
				adjListCopy.get(i).add(verticesCopy.get(adjList.get(i).get(j).node));
			}
		}
	}

	// Helper print functions to print the results in a CSV format 
	// so that we can import it to an excel for charting
	public static void printTimeCSV(int verticesCount, double timeDsatur, double timeDsaturHeap, double timeDsaturTree)
	{
		System.out.println(verticesCount + ", " + timeDsatur + ", " + timeDsaturHeap + ", " + timeDsaturTree);
	}
	
	public static void printColorsCSV(int verticesCount, double cDsatur, double cHeapDsatur, double cTreeDsatur)
	{
		System.out.println(verticesCount + ", " + cDsatur + ", " + cHeapDsatur + ", " + cTreeDsatur);
	}
}


class DsaturPlain
{

	Graph graph;
	List<Node> vertices;
	List<LinkedList<Node>> adjList;

	List<Integer> count = new ArrayList<Integer>();
	List<Integer> vertexColor = new ArrayList<Integer>();
	Scanner scanInput = new Scanner(System.in);
	long timeTakenInMs = -1;
	int colorMax = 0;


	public void storeGraphObj(Graph graphObj)
	{
		graph = graphObj;
	}

	public void populateGraph(int vertexCount)
	{
		vertices = new ArrayList<>();
		adjList = new ArrayList<>();
		graph.populateGraph(vertices, adjList); 
	}

	// stores the number of edges for each vertices and stores the labels of the vertices in seperate list
	public void storecount(int vertexCount)
	{
		for(int i=0; i<adjList.size(); i++)
		{
			int x = adjList.get(i).size();
			count.add(x);
		}
	}


	// finds the maximum Adjacent Degree
	public int maxAdjDegree(int vertexCount)
	{
		int max=0; 
		int vertexNumber =0;
		for(int i=0; i<vertexCount; i++)
		{
			if(max < count.get(i))
			{
				max = count.get(i);
				vertexNumber = i;
			}
		}
		//System.out.println();
		//System.out.println("Maximum number of edges is "+max +" in vertex "+vertexNumber);
		return vertexNumber;	

	}

	// Returns the vertex with maximum saturation degree
	// This method has a complexity of O(V * E) as we go through each vertex and 
	// traverse each vertice's neighbour to find unique color count
	public int saturationDegree(int vertexCount)
	{
		int maxSatDeg = -1, vertexNumber = -1; int x;

		int size = -1;
		int adjDegMax = -1;

		for(int i=0; i<vertexCount;i++)
		{
			if(vertexColor.get(i) == 0)
			{
				x = adjList.get(i).size();
				if(x == 0)
				{
					size = 0;

				}
				Set<Integer> satDeg = new HashSet<Integer>();
				for(int j=0; j<x; j++)
				{
					satDeg.add(vertexColor.get(adjList.get(i).get(j).node));
				}
				size = satDeg.size();
				if(maxSatDeg < size)
				{
					maxSatDeg = size;
					vertexNumber = i;
					adjDegMax = x;

				}
				if((maxSatDeg == size) && (x > adjDegMax))
				{
					adjDegMax = x;
					maxSatDeg = size;
					vertexNumber = i;
				}
			}
		}
		//System.out.println();
		//System.out.println("Maximum Saturation Degree is "+maxSatDeg+" of vertex "+vertexNumber);
		return vertexNumber;
	}

	// To color the vertices
	public void colorVertex(int vertexCount)
	{
		// Initialize color 0 to all vertices at the start
		for(int i=0; i<vertexCount; i++)
		{
			vertexColor.add(0);

		}
		int first = maxAdjDegree(vertexCount);
		vertexColor.set(first,1);

		int vertexToBeColored = saturationDegree(vertexCount);

		// After coloring each vertex we call the saturationDegree method to figure out next vertex to be colored
		// This has a complexity of O(V * V * E). Refer to above saturationDegree method for its complexity.
		while(vertexToBeColored != -1)
		{
			addColorToVertex(vertexToBeColored,vertexCount); // calls a function to add color to the other vertices
			vertexToBeColored = saturationDegree(vertexCount);
		}

	}

	// To add the lower numbered color to a vertex.
	// This method has a complexity of O(V + E) as we use a HashSet (with O(1) look up) to store unique colors for 
	// vertexToBeColored's neighbors. Worst case we would be going through V colors to figure out the resulting color.   
	public void addColorToVertex(int vertexToBeColored, int vertexCount)
	{
		int count = adjList.get(vertexToBeColored).size();
		Set<Integer> checkColor = new HashSet<Integer>();


		for(int j=0; j<count; j++)
		{
			checkColor.add(vertexColor.get(adjList.get(vertexToBeColored).get(j).node));
		}

		for(int i=1; i<=vertexCount; i++)
		{
			if(checkColor.contains(i))
			{
				continue;
			}
			else
			{
				vertexColor.set(vertexToBeColored,i);
				if (colorMax < i)colorMax = i;
				break;
			}
		}
		
	}

	//displays the color of each vertex
	public void displayColor(int vertexCount)
	{
		System.out.println();
		System.out.println("Dsatur Adjacency List Algorithm Coloring");
		for(int i=0; i<vertexCount; i++)
		{
			System.out.println("Vertex number "+i +" Color : "+vertexColor.get(i));
		}
	}

	public GraphResult dsatur(int vertexCount, Graph graphObj)
	{
		storeGraphObj(graphObj);
		populateGraph(vertexCount);

		long startTime = System.currentTimeMillis();
		storecount(vertexCount);
		colorVertex(vertexCount);
		long endTime = System.currentTimeMillis();
		//displayColor(vertexCount);
		return new GraphResult(endTime - startTime, colorMax);

	}
}

class DsaturOptimized
{
	Graph graph;
	List<Node> vertices;
	List<LinkedList<Node>> adjList;
	long timeTakenInMs = -1;
	int colorMax = 0;

	// Max priority Heap to store the saturation degree of vertices
	PriorityQueue<Node> satDegQueue = new PriorityQueue<>();

	public void storeGraphObj(Graph graphObj)
	{
		graph = graphObj;
	}

	public void populateGraph(int vertexCount)
	{
		vertices = new ArrayList<>();
		adjList = new ArrayList<>();
		graph.populateGraph(vertices, adjList); 
	}

	// stores the number of edges for each vertices and stores the labels of the vertices in seperate list
	public void storecount(int vertexCount)
	{
		for(int i=0; i<vertexCount; i++)
		{
			int x = adjList.get(i).size();
			vertices.get(i).adjDeg = x;
			satDegQueue.add(vertices.get(i));
		}

	}

	// To udpate the saturation color we keep a hash set to collect unique neighbor color. O(E).
	int getUpdatedSaturationDegree(Node satNode)
	{
		int size = satNode.adjDeg;
		Set<Integer> tempSet = new HashSet<>();
		for (int i = 0; i < size; ++i)
		{
			tempSet.add(adjList.get(satNode.node).get(i).color);
		}
		return tempSet.size();
	}

	
	// To color the vertex
	// This method has a complexity of O(V * (logV +  (V+E) + (E* (E+logV)))) => O(V*logV*E + V*E*E)
	public void colorVertex(int vertexCount)
	{
		while(!satDegQueue.isEmpty()) {
			// Poll from priority queue in O(logV)
			Node nextNode = satDegQueue.poll();

			int nextVertexDeg = nextNode.adjDeg;
			int nextVertex = nextNode.node;

			// Finding the color to be used for nextVertex O(V+E). Refer to addColorToVertex method for complexity
			addColorToVertex(nextVertex, vertexCount);

			// Update the priority queue for next Vertex's adjacent nodes
			// Complexity O(E * logV) as we traverse for each adjacent vertex, update saturation degree and update the heap. 
			for(int i=0; i<nextVertexDeg; i++)
			{
				if(adjList.get(nextVertex).get(i).color == 0)
				{
					Node nextAdj = adjList.get(nextVertex).get(i);
					nextAdj.satDeg = getUpdatedSaturationDegree(nextAdj);
					satDegQueue.remove(nextAdj);
					satDegQueue.add(nextAdj);
				}
			}
		}
	}

	// To add the lower numbered color to a vertex.
	// This method has a complexity of O(V + E) as we use a HashSet (with O(1) look up) to store unique colors for 
	// vertexToBeColored's neighbors. Worst case we would be going through V colors to figure out the resulting color.
	public void addColorToVertex(int vertexToBeColored, int vertexCount)
	{
		int count = adjList.get(vertexToBeColored).size();
		Set<Integer> checkColor = new HashSet<Integer>();

		for(int j=0; j<count; j++)
		{
			checkColor.add(adjList.get(vertexToBeColored).get(j).color);
		}

		for(int i=1; i<=vertexCount; i++)
		{
			if(checkColor.contains(i))
			{
				continue;
			}
			else
			{
				vertices.get(vertexToBeColored).color = i;
				if (i > colorMax)colorMax = i;
				break;
			}
		}
	
	}

	//displays the color of each vertex
	public void displayColor(int vertexCount)
	{
		System.out.println();
		System.out.println("Dsatur PriorityQueue Algorithm Coloring");
		for(int i=0; i<vertexCount; i++)
		{
			System.out.println("Vertex number "+i +" Color : "+vertices.get(i).color);
		}
	}

	public GraphResult dsatur(int vertexCount, Graph graphObj)
	{
		storeGraphObj(graphObj);
		populateGraph(vertexCount);

		long startTime = System.currentTimeMillis();
		storecount(vertexCount);
		colorVertex(vertexCount);
		long endTime = System.currentTimeMillis();
		//displayColor(vertexCount);
		return new GraphResult(endTime - startTime, colorMax);
	}
}

class DsaturOptimizedTree
{

	Graph graph;
	List<Node> vertices;
	List<LinkedList<Node>> adjList;

	Scanner scanInput = new Scanner(System.in);
	
	// Java implements TreeSet internally using a Red-Black tree. We use it to store saturation degree
	// JDK library reference https://docs.oracle.com/javase/7/docs/api/java/util/TreeSet.html
	// https://docs.oracle.com/javase/7/docs/api/java/util/TreeMap.html
	TreeSet<Node> satDegTree = new TreeSet<>();

	long timeTakenInMs = -1;
	int colorMax = 0;

	public void storeGraphObj(Graph graphObj)
	{
		graph = graphObj;
	}

	public void populateGraph(int vertexCount)
	{
		vertices = new ArrayList<>();
		adjList = new ArrayList<>();
		graph.populateGraph(vertices, adjList); 
	}

	// stores the number of edges for each vertices and stores the labels of the vertices in seperate list
	public void storecount(int vertexCount)
	{
		for(int i=0; i<vertexCount; i++)
		{
			int x = adjList.get(i).size();
			vertices.get(i).adjDeg = x;
			satDegTree.add(vertices.get(i));
		}

	}

	int getUpdatedSaturationDegree(Node satNode)
	{
		int size = satNode.adjDeg;
		Set<Integer> tempSet = new HashSet<>();
		for (int i = 0; i < size; ++i)
		{
			tempSet.add(adjList.get(satNode.node).get(i).color);
		}
		return tempSet.size();
	}	

	// To color the vertex
	// This method has a complexity of O(V * (logV +  (V+E) + (E*E +  2 * logV * E))) => O((V*E*E) + (V*E*logV))  
	public void colorVertex(int vertexCount)
	{
		while(!satDegTree.isEmpty()) {

			Node nextNode = satDegTree.pollFirst();
			int nextVertexDeg = nextNode.adjDeg;
			int nextVertex = nextNode.node;

			addColorToVertex(nextVertex, vertexCount);

			// Update the tree set for next Vertex's adjacent nodes
			for(int i=0; i<nextVertexDeg; i++)
			{
				if(adjList.get(nextVertex).get(i).color == 0)
				{
					Node nextAdj = adjList.get(nextVertex).get(i);
					
					// We need to remove and add the updated vertex to the tree.
					// Otherwise the sorted nature of the tree would be messed up. Unlinke Heap, complexity here is O(logV + E)
					satDegTree.remove(nextAdj);
					nextAdj.satDeg = getUpdatedSaturationDegree(nextAdj);
					satDegTree.add(nextAdj);

				}
			}
		}
	}

	//to add the lower numbered color to a vertex
	public void addColorToVertex(int vertexToBeColored, int vertexCount)
	{
		int count = adjList.get(vertexToBeColored).size();
		Set<Integer> checkColor = new HashSet<Integer>();

		for(int j=0; j<count; j++)
		{
			checkColor.add(adjList.get(vertexToBeColored).get(j).color);
		}

		for(int i=1; i<=vertexCount; i++)
		{
			if(checkColor.contains(i))
			{
				continue;
			}
			else
			{
				vertices.get(vertexToBeColored).color = i;
				if (i > colorMax)colorMax = i;
				break;
			}
		}
	
	}

	//displays the color of each vertex
	public void displayColor(int vertexCount)
	{
		System.out.println();
		System.out.println("Dsatur TreeSet Algorithm Coloring");
		for(int i=0; i<vertexCount; i++)
		{
			System.out.println("Vertex number "+i +" Color : "+vertices.get(i).color);
		}
	}

	public GraphResult dsatur(int vertexCount, Graph graphObj)
	{
		storeGraphObj(graphObj);
		populateGraph(vertexCount);

		long startTime = System.currentTimeMillis();
		storecount(vertexCount);
		colorVertex(vertexCount);
		long endTime = System.currentTimeMillis();
		//displayColor(vertexCount);
		return new GraphResult(endTime - startTime, colorMax);
	}
}

public class TCSS543
{
	public static void main(String []args)
	{
		int numberOfVertex;
		int input;
		int totalNumberOfGraphs = 100;
		int maxNumberOfVertices = 200;
			
		//Gather each data structure results for time and colors required in following maps
		Map<Integer, Double> dPlainTimeRequiredResults = new HashMap<>();
		Map<Integer, Double> dHeapTimeRequiredResults = new HashMap<>();
		Map<Integer, Double> dTreeTimeRequiredResults = new HashMap<>();
		Map<Integer, Double> dPlainColorRequiredResults = new HashMap<>();
		Map<Integer, Double> dHeapColorRequiredResults = new HashMap<>();
		Map<Integer, Double> dTreeColorRequiredResults = new HashMap<>();
		
		double[] densityArray = {0.73, 0.61, 0.44, 0.26};

		// For each density value generate 100 graphs for each with vertices 10, 15, 20, 25, ..., 100
		for (int j = 0; j < densityArray.length; j++)
		{
			for (numberOfVertex = 10; numberOfVertex <= maxNumberOfVertices; numberOfVertex +=5)
			{	
			
				int dsaturPlainTimeInMs = 0;
				int dsaturHeapTimeInMs = 0;
				int dsaturTreeTimeInMs = 0;
				int dsaturPlainColorRequired = 0;
				int dsaturHeapColorRequired = 0;
				int dsaturTreeColorRequired = 0;

				for(int i=1; i<=totalNumberOfGraphs; i++)
				{
					TCSS543 mainObj = new TCSS543();
					Graph graph = new Graph();
					DsaturPlain dsaturPlain = new DsaturPlain();
					DsaturOptimized dsaturOpt = new DsaturOptimized();
					DsaturOptimizedTree dsaturOptTree = new DsaturOptimizedTree();
			
					graph.generateRandomGraph(numberOfVertex, densityArray[j]);
					//graph.printGraph(numberOfVertex);
					GraphResult dsaturPlainResult = dsaturPlain.dsatur(numberOfVertex,graph);
					GraphResult dsaturHeapResult = dsaturOpt.dsatur(numberOfVertex,graph);
					GraphResult dsaturTreeResult = dsaturOptTree.dsatur(numberOfVertex,graph);

					dsaturPlainTimeInMs += dsaturPlainResult.timeRequiredToColorInMs;
					dsaturHeapTimeInMs += dsaturHeapResult.timeRequiredToColorInMs;
					dsaturTreeTimeInMs += dsaturTreeResult.timeRequiredToColorInMs;

					dsaturPlainColorRequired += dsaturPlainResult.maxColorRequired;
					dsaturHeapColorRequired += dsaturHeapResult.maxColorRequired;
					dsaturTreeColorRequired += dsaturTreeResult.maxColorRequired;
				}
				dPlainTimeRequiredResults.put(numberOfVertex, dsaturPlainTimeInMs/(double)totalNumberOfGraphs);
				dHeapTimeRequiredResults.put(numberOfVertex, dsaturHeapTimeInMs/(double)totalNumberOfGraphs);
				dTreeTimeRequiredResults.put(numberOfVertex, dsaturTreeTimeInMs/(double)totalNumberOfGraphs);

				dPlainColorRequiredResults.put(numberOfVertex, dsaturPlainColorRequired/(double)totalNumberOfGraphs);
				dHeapColorRequiredResults.put(numberOfVertex, dsaturHeapColorRequired/(double)totalNumberOfGraphs);
				dTreeColorRequiredResults.put(numberOfVertex, dsaturTreeColorRequired/(double)totalNumberOfGraphs);
			}
			System.out.println("Density = " + densityArray[j]);
			System.out.println("VertexCount, DsaturPlainTimeInMs, DsaturHeapTimeInMs, DsaturTreeTimeInMs");
			for (int i=10; i <=maxNumberOfVertices; i+=5)
			{
				Graph.printTimeCSV(i, dPlainTimeRequiredResults.get(i), dHeapTimeRequiredResults.get(i), dTreeTimeRequiredResults.get(i));
			}
			System.out.println();
			System.out.println("VertexCount, DsaturPlainColorsUsed, DsaturHeapColorsUsed, DsaturTreeColorsUsed");
			for (int i=10; i <=maxNumberOfVertices; i+=5)
			{
				Graph.printColorsCSV(i, dPlainColorRequiredResults.get(i), dHeapColorRequiredResults.get(i), dTreeColorRequiredResults.get(i));
			}
		}
	}
}

// A utility class to encapsulate the results of each run
class GraphResult
{
	public int maxColorRequired;
	public long timeRequiredToColorInMs;

	public GraphResult(long time, int colors)
	{
		this.maxColorRequired = colors;
		this.timeRequiredToColorInMs = time;
	}
}

// A Node instance represents the attributes of a vertex in the graph.
class Node implements Comparable <Node>
{
	public int node;
	public int satDeg;
	public int adjDeg;
	public int color;

	public Node(int node, int satDeg, int adjDeg, int color)
	{
		this.node = node;
		this.satDeg = satDeg;
		this.adjDeg = adjDeg;
		this.color = color;
	}

	public Node(Node node2)
	{
		this.node = node2.node;
		this.satDeg = node2.satDeg;
		this.adjDeg = node2.adjDeg;
		this.color = node2.color;
	}
	
	// Comparator used by priority Heap and Red-Black tree to order the nodes in underlying data structures  
	@Override
	public int compareTo(Node other)
	{
		if (this.satDeg < other.satDeg || 
		   (this.satDeg == other.satDeg && this.adjDeg < other.adjDeg) ||
		   (this.satDeg == other.satDeg && this.adjDeg == other.adjDeg && this.node > other.node)) {
			return 1;
		} else if (this.satDeg > other.satDeg || 
			  (this.satDeg == other.satDeg && this.adjDeg > other.adjDeg) ||
                          (this.satDeg == other.satDeg && this.adjDeg == other.adjDeg && this.node < other.node)) {
			return -1;
		} else {
			return 0;
		}
	}
}

