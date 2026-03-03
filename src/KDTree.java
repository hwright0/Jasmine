/*
 * Data structure for fast k-nearest neighbor queries in variant sets
 * For a given query, the k closest points to it in the dataset will be reported,
 * breaking ties by variant ID to ensure deterministic behavior.
 * 
 * We assume variants are 2-D points; nearness is based on Euclidean distance or its generalizations.
 * 
 * Uses algorithm described here: 
 * https://courses.cs.washington.edu/courses/cse599c1/13wi/slides/lsh-hashkernels-annotated.pdf
 */

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Stack;

public class KDTree 
{
	Node root;
	int cnt;
	int K;
	
	int n;
	
	/*
	 * Initializes a KD-tree from a list of variants
	 */
	public KDTree(Variant[] p) 
	{
		n = p.length;
		K = 2;
		LinkedList<Node> list = new LinkedList<Node>();
		for (Variant q : p) list.add(new Node(q));
		root = build(list, 0);//buildNonrecursive(list, 0).get(0);
	}
	
	public KDTree(Variant[] p, boolean recursive) 
	{
		n = p.length;
		K = 2;
		LinkedList<Node> list = new LinkedList<Node>();
		for (Variant q : p) list.add(new Node(q));
		root = recursive ? build(list, 0) : buildNonrecursive(list).get(0);
	}
	
	private Node build(LinkedList<Node> p, int depth) 
	{
		if (p.size() == 0) return null;
		Node pivot = p.remove();
		
		// Sort the points into left and right subtrees based on current split dimension
		LinkedList<Node> left = new LinkedList<Node>();
		LinkedList<Node> right = new LinkedList<Node>();
		while (!p.isEmpty()) 
		{
			if (p.peek().planes[depth % K] < pivot.planes[depth % K])
				left.add(p.remove());
			else
				right.add(p.remove());
		}
		pivot.children[0] = build(left, depth + 1);
		pivot.children[1] = build(right, depth + 1);
		
		return pivot;
	}
	
	/*
	 * Build the data structure from a list of points without recursion
	 * This avoids stack overflow issues caused by larger datasets
	 */
	private ArrayList<Node> buildNonrecursive(LinkedList<Node> p) 
	{
		ArrayList<Node> nodeList = new ArrayList<Node>();
		if(p.size() == 0)
		{
			return null;
		}
		
		// The stack of node lists to process (in place of recursive calls)
		ArrayDeque<LinkedList<Node>> toProcess = new ArrayDeque<LinkedList<Node>>();
		ArrayDeque<Integer> parents = new ArrayDeque<Integer>();
		ArrayDeque<Integer> depths = new ArrayDeque<Integer>();
		ArrayDeque<Integer> parentsides = new ArrayDeque<Integer>();
		
		// Initialize root to null to be filled 
		//nodeList.add(res);
		toProcess.addFirst(p);
		parents.addFirst(-1); // This is not actually the parent of the root, but it will be ignored anyways
		depths.addFirst(0);
		parentsides.addFirst(-1);
		
		while(!toProcess.isEmpty())
		{
			// Get information for processing this node from stacks
			LinkedList<Node> pcur = toProcess.pollFirst();
			int parentcur = parents.pollFirst();
			int depthcur = depths.pollFirst();
			int parentsidecur = parentsides.pollFirst();
			
			// Sort by the current split dimension and pick the median as pivot.
			// This produces a balanced tree (O(log n) depth) instead of the
			// degenerate O(n) depth that results from always picking the first
			// element when variants are pre-sorted by position.
			final int splitDim = depthcur % K;
			ArrayList<Node> sortedNodes = new ArrayList<Node>(pcur);
			sortedNodes.sort((x, y) -> Double.compare(x.planes[splitDim], y.planes[splitDim]));
			int medianIdx = sortedNodes.size() / 2;
			Node pivot = sortedNodes.get(medianIdx);
			
			// Partition directly from the sorted list: indices < median → left, > median → right
			LinkedList<Node> left = new LinkedList<Node>();
			LinkedList<Node> right = new LinkedList<Node>();
			for(int si = 0; si < sortedNodes.size(); si++)
			{
				if(si == medianIdx) continue;
				if(si < medianIdx)
					left.add(sortedNodes.get(si));
				else
					right.add(sortedNodes.get(si));
			}
			sortedNodes = null; // allow GC
			
			//pcur.clear();
			
			// Update this node's parent's child-pointer to this node.
			if(parentsidecur != -1)
			{
				nodeList.get(parentcur).children[parentsidecur] = pivot;
			}
			
			pivot.children[0] = null;
			pivot.children[1] = null;
			nodeList.add(pivot);
			
			// Add right child to processing stack
			if(right.size() > 0)
			{
				toProcess.addFirst(right);
				parents.addFirst(nodeList.size() - 1);
				parentsides.addFirst(1);
				depths.addFirst(depthcur + 1);
			}
			
			// Add left child to processing stack
			if(left.size() > 0)
			{
				toProcess.addFirst(left);
				parents.addFirst(nodeList.size() - 1);
				parentsides.addFirst(0);
				depths.addFirst(depthcur + 1);
			}
		}
		
		return nodeList;
	}
	
	/*
	 * Used to make sure two KD-trees are the same
	 */
	static boolean compare(String pref, Node a, Node b)
	{
		if(a == null && b != null)
		{
			System.out.println(pref + " only a is null");
			return true;
		}
		if(b == null && a != null)
		{
			System.out.println(pref + " only b is null");
			return true;
		}
		if(a == null && b == null)
		{
			return false;
		}
		if(a.planes[0] != b.planes[0] || a.planes[1] != b.planes[1])
		{
			System.out.println(pref + " diff value: " + a.planes[0] + " " + a.planes[1] + " " + b.planes[0] + " " + b.planes[1]);
			return true;
		}
		
		boolean leftDiff = compare(pref + "L", a.children[0], b.children[0]);
		if(leftDiff)
		{
			return true;
		}
		else
		{
			return compare(pref + "R", a.children[1], b.children[1]);
		}
		
	}
	
	/*
	 * Gets the k nearest neighbors for a query variant.
	 * Thread-safe: all query state is local to this call.
	 */
	public Variant[] kNearestNeighbor(Variant p, int k) {
		Node searchNode = new Node(p);
		PriorityQueue<Candidate> bestCandidates = new PriorityQueue<Candidate>();
		search(root, 0, searchNode, bestCandidates, k);
		Variant[] res = new Variant[bestCandidates.size()];
		int idx = res.length - 1;
		while(!bestCandidates.isEmpty())
		{
			res[idx--] = bestCandidates.poll().v;
		}
		return res;
	}
	
	/*
	 * Search the subtree rooted at cur for candidate points in the set of query's k-nearest neighbors.
	 * All state is passed as parameters so multiple threads can call this concurrently.
	 */
	private void search(Node ocur, int odepth, Node search, PriorityQueue<Candidate> best, int querySize) {
		Stack<Node> curs = new Stack<Node>();
		Stack<Integer> depths = new Stack<Integer>();
		Stack<Boolean> processedBest = new Stack<Boolean>();
		curs.add(ocur);
		depths.add(odepth);
		processedBest.push(false);
		while(!curs.isEmpty())
		{
			Node cur = curs.pop();
			int depth = depths.pop();
			boolean bestDone = processedBest.pop();
			
			if(cur == null) continue;
			
			int betterChild = (int) Math.signum(search.planes[depth % K] - cur.planes[depth % K]) < 0 ? 0 : 1;
			
			if(!bestDone)
			{
				curs.add(cur);
				depths.add(depth);
				processedBest.add(true);
				curs.add(cur.children[betterChild]);
				depths.add(depth+1);
				processedBest.add(false);
				continue;
			}
			Candidate toAdd = new Candidate(cur.p, cur.p.distance(search.p));
			if (best == null || best.size() < querySize || toAdd.compareTo(best.peek()) > 0) 
			{
				if(best.size() == querySize)
				{
					best.poll();
				}
				best.add(toAdd);
			}
			if (best.size() < querySize || Math.abs(search.planes[depth % K] - cur.planes[depth % K]) < best.peek().dist)
			{
				curs.add(cur.children[1 - betterChild]);
				depths.add(depth+1);
				processedBest.add(false);
			}
		}
	}
	
	/*
	 * A node of the KD tree
	 * Each node has a variant, storing alongside it its values along the split planes, as well as two (possibly null) children
	 */
	private class Node {
		Node[] children;
		Variant p;
		double[] planes;
		public Node(Variant pp) 
		{
			p = pp;
			planes = new double[K];
			planes[0] = p.start;
			planes[1] = p.end; // add additional dimensions as necessary
			children = new Node[2];
			children[0] = null;
			children[1] = null;
		}
	}
	
	/*
	 * Candidate k-nearest neighbor of the current query point
	 */
	private static class Candidate implements Comparable<Candidate>
	{
		Variant v;
		double dist;
		Candidate(Variant v, double dist)
		{
			this.v = v;
			this.dist = dist;
		}
		public int compareTo(Candidate o)
		{
			if(Math.abs(dist - o.dist) > 1e-9) return Double.compare(o.dist, dist);
			if(v.hash != o.v.hash) return o.v.hash - v.hash;
			return o.v.id.compareTo(v.id);
			
		}
	}
}
