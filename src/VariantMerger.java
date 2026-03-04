/*
 * Main interface for merging variants.  This assumes all variants are on the
 * same chromosome, and have the same type/strand if that separation is desired.
 */

import java.util.ArrayList;
import java.util.PriorityQueue;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class VariantMerger
{
	// An array of all of the variants to be considered
	Variant[] data;
	
	// The number of total variants
	int n;
	
	// A forest in which connected components will represent merged groups
	Forest forest;
	
	// A KD-tree data structure for fast k-nearest-neighbors queries
	KDTree knn;
	
	// Indices of variants in each group, used for more advanced distance checks like clique and centroid
	ArrayList<Integer>[] merged;

	@SuppressWarnings("unchecked")
	public VariantMerger(Variant[] data)
	{
		n = data.length;
		
		forest = new Forest(data);
		knn = new KDTree(data, false);
		
		for(int i = 0; i<n; i++) data[i].index = i;
		this.data = data;
		
		if(Settings.CENTROID_MERGE || Settings.CLIQUE_MERGE)
		{
			merged = new ArrayList[n];
			for(int i = 0; i<n; i++)
			{
				merged[i] = new ArrayList<Integer>();
				merged[i].add(i);
			}
		}
	}
	
	/*
	 * Helper function to convert an ArrayList to an array to make the constructor more flexible
	 */
	static Variant[] listToArray(ArrayList<Variant> data)
	{
		int length = data.size();
		Variant[] asArray = new Variant[length];
		for(int i = 0; i<length; i++)
		{
			asArray[i] = data.get(i);
		}
		return asArray;
	}
	
	/*
	 * Alternate constructor which takes a list instead of an array
	 */
	public VariantMerger(ArrayList<Variant> data)
	{
		this(listToArray(data));
	}
	
	/*
	 * Runs the core algorithm for building the implicit merging graph and
	 * performing merging
	 */
	void runMerging()
	{
		if(n == 1)
		{
			return;
		}

		// Delegate to hierarchical merging for large dense graphs when requested
		if(Settings.HIERARCHICAL_BATCH_SIZE > 0 && n > Settings.HIERARCHICAL_BATCH_SIZE
			&& !Settings.CLIQUE_MERGE && !Settings.CENTROID_MERGE)
		{
			runHierarchicalMerging();
			return;
		}

		long mergeStartTime = System.currentTimeMillis();
		System.out.printf("[Merge] Graph with %,d variants: building kNN graph (%d threads)%n",
			n, Math.max(1, Settings.THREADS));

		// For each variant v, how many of its nearest neighbors have had their edges
		// from v considered already.  
		int[] countEdgesProcessed = new int[n];
		
		// nearestNeighbors will be used as a cache to store the next few nearest neighbors
		// The purpose of this is to prevent performing a new KNN-query every time an edge
		// is considered, but instead a logarithmic number of times.
		Variant[][] nearestNeighbors = new Variant[n][];
		
		// A heap of edges to be processed in non-decreasing order of distance
		PriorityQueue<Edge> toProcess = new PriorityQueue<Edge>();
		
		// Populate nearestNeighbors[] in parallel: each query is independent since
		// KDTree.kNearestNeighbor is now fully thread-safe (no shared mutable state).
		// Each slot nearestNeighbors[i] is written by exactly one thread so the array
		// itself needs no synchronization.
		int numThreads = Math.max(1, Settings.THREADS);
		ExecutorService pool = Executors.newFixedThreadPool(numThreads);
		CountDownLatch latch = new CountDownLatch(n);
		for(int i = 0; i < n; i++)
		{
			final int fi = i;
			pool.submit(() -> {
				try
				{
					nearestNeighbors[fi] = knn.kNearestNeighbor(data[fi], 4);
				}
				finally
				{
					latch.countDown();
				}
			});
		}
		try { latch.await(); } catch(InterruptedException ie) { Thread.currentThread().interrupt(); }
		pool.shutdown();

		System.out.printf("[Merge] kNN done in %.1fs%n",
			(System.currentTimeMillis() - mergeStartTime) / 1000.0);

		// Build the initial edge heap sequentially now that nearestNeighbors[] is fully populated
		for(int i = 0; i<n; i++)
		{
			if(nearestNeighbors[i] == null || nearestNeighbors[i].length == 0) { continue; }
			int maxDistAllowed = Math.max(data[i].maxDist, nearestNeighbors[i][0].maxDist);
			if(Settings.REQUIRE_MUTUAL_DISTANCE)
			{
				maxDistAllowed = Math.min(data[i].maxDist, nearestNeighbors[i][0].maxDist);
			}
			if(data[i].distance(nearestNeighbors[i][0]) < maxDistAllowed + 1e-9)
			{
				toProcess.add(new Edge(i, nearestNeighbors[i][0].index, data[i].distance(nearestNeighbors[i][0])));
			}
			countEdgesProcessed[i]++;
		}
		
		long edgesProcessed = 0;
		long nextReport = 500_000L;
		long edgeLoopStart = System.currentTimeMillis();

		while(!toProcess.isEmpty())
		{
			edgesProcessed++;
			if(edgesProcessed >= nextReport)
			{
				System.out.printf("[Merge] Edge loop: %,d edges, heap=%,d, elapsed=%.1fs%n",
					edgesProcessed, toProcess.size(),
					(System.currentTimeMillis() - edgeLoopStart) / 1000.0);
				nextReport += 500_000L;
			}
			Edge e = toProcess.poll();
			boolean valid = forest.canUnion(e.from, e.to);
			if(valid)
			{
				int fromRoot = forest.find(e.from);
				int toRoot = forest.find(e.to);
				if(Settings.CLIQUE_MERGE || Settings.CENTROID_MERGE)
				{
					// Makes sure all newly merged variant pairs are within the distance threshold
					if(Settings.CLIQUE_MERGE)
					{
						for(int i = 0; i<merged[fromRoot].size() && valid; i++)
						{
							Variant candidateFrom = data[merged[fromRoot].get(i)];
							for(int j = 0; j<merged[toRoot].size() && valid; j++)
							{
								Variant candidateTo = data[merged[toRoot].get(j)];
								int maxDistAllowed = Math.max(candidateFrom.maxDist, candidateTo.maxDist);
								if(Settings.REQUIRE_MUTUAL_DISTANCE)
								{
									maxDistAllowed = Math.min(data[i].maxDist, nearestNeighbors[i][0].maxDist);
								}
								if(candidateFrom.distance(candidateTo) > maxDistAllowed + 1e-9)
								{
									valid = false;
								}
							}
						}
					}
					// Make sure everything being merged can be merged with their overall centroid
					else if(Settings.CENTROID_MERGE)
					{
						double avgStart = 0.0, avgEnd = 0.0;
						for(int i = 0; i<merged[fromRoot].size(); i++)
						{
							Variant v = data[merged[fromRoot].get(i)];
							avgStart += v.start;
							avgEnd += v.end;
						}
						for(int i = 0; i<merged[toRoot].size(); i++)
						{
							Variant v = data[merged[toRoot].get(i)];
							avgStart += v.start;
							avgEnd += v.end;
						}
						avgStart /= merged[fromRoot].size() + merged[toRoot].size();
						avgEnd /= merged[fromRoot].size() + merged[toRoot].size();
						
						for(int i = 0; i<merged[fromRoot].size() && valid; i++)
						{
							Variant v = data[merged[fromRoot].get(i)];
							valid &= v.distFromPoint(avgStart, avgEnd) <= v.maxDist + 1e-9;
						}
						for(int i = 0; i<merged[toRoot].size() && valid; i++)
						{
							Variant v = data[merged[toRoot].get(i)];
							valid &= v.distFromPoint(avgStart, avgEnd) <= v.maxDist + 1e-9;
						}
					}
					if(valid)
					{
						forest.union(fromRoot, toRoot);
						if(forest.map[fromRoot] < 0)
						{
							// The first variant is the new root of the union-find component
							for(int x : merged[toRoot])
							{
								merged[fromRoot].add(x);
							}
						}
						else
						{
							for(int x : merged[fromRoot])
							{
								merged[toRoot].add(x);
							}
						}
					}
				}
				
				// Two variants are being merged here - nothing needs to be done
				else
				{
					forest.union(e.from, e.to);
				}
			}
			
			while(true)
			{
				
				// If we already used the stored neighbors, query again for twice as many
				if(countEdgesProcessed[e.from] >= nearestNeighbors[e.from].length)
				{
					nearestNeighbors[e.from] = knn.kNearestNeighbor(data[e.from], 2 * nearestNeighbors[e.from].length);
				}
				
				// If we tried to get more and didn't find anymore, then we are done with this variant
				if(countEdgesProcessed[e.from] >= nearestNeighbors[e.from].length)
				{
					break;
				}
				Variant candidateTo = nearestNeighbors[e.from][countEdgesProcessed[e.from]];
				
				// This edge was invalid because of distance from the query, so stop looking at any edges 
				// since they'll only get farther away
				int maxDistAllowed = Math.max(data[e.from].maxDist, candidateTo.maxDist);
				if(Settings.REQUIRE_MUTUAL_DISTANCE)
				{
					maxDistAllowed = Math.min(data[e.from].maxDist, candidateTo.maxDist);
				}
				
				if(data[e.from].distance(candidateTo) > data[e.from].maxDist + 1e-9)
				{
					break;
				}
				
				else if(data[e.from].distance(candidateTo) > maxDistAllowed + 1e-9)
				{
					countEdgesProcessed[e.from]++;
					continue;
				}
				
				// If edge was invalid because of coming from the same sample, ignore it and try the next one
				else if(!Settings.ALLOW_INTRASAMPLE && data[e.from].sample == candidateTo.sample)
				{
					toProcess.add(new Edge(e.from, candidateTo.index, data[e.from].distance(candidateTo)));
					countEdgesProcessed[e.from]++;
					break;
				}
				
				// If sequences weren't similar enough for two insertions, ignore and try again
				else if(!data[e.from].passesStringSimilarity(candidateTo))
				{
					countEdgesProcessed[e.from]++;
					continue;
				}
				
				else if(!data[e.from].passesOverlap(candidateTo))
				{
					countEdgesProcessed[e.from]++;
					continue;
				}
				
				// The next edge is something we want to consider since it is close enough and goes to a
				// different sample
				else
				{
					toProcess.add(new Edge(e.from, candidateTo.index, data[e.from].distance(candidateTo)));
					countEdgesProcessed[e.from]++;
					break;
				}
			}
		}
		System.out.printf("[Merge] Complete: %,d edges processed, total=%.1fs%n",
			edgesProcessed, (System.currentTimeMillis() - mergeStartTime) / 1000.0);
	}

	/*
	 * Hierarchical merging: splits variants into batches, merges each batch independently
	 * in parallel (Phase 1), then merges the batch-level representatives in a second pass
	 * (Phase 2), and finally propagates the combined groupings back into this.forest
	 * (Phase 3).
	 *
	 * Phase 1 batches are formed by splitting the position-sorted input array into
	 * consecutive windows of Settings.HIERARCHICAL_BATCH_SIZE variants.  Because the
	 * data is sorted, each window covers a contiguous genomic region, so variants that
	 * should merge almost always land in the same batch.
	 *
	 * Same-sample filtering (ALLOW_INTRASAMPLE=false):
	 *   After Phase 1 each batch representative stands for many samples.  Before Phase 2
	 *   runs, the representative's sampleSets entry in finalVM.forest is replaced with
	 *   the full union of all samples in its batch group.  This ensures that Phase 2
	 *   correctly rejects merges whose underlying groups share a sample.
	 *
	 * Known approximation (inherent to any hierarchical scheme):
	 *   A cluster whose union-find root is not the variant closest to the batch boundary
	 *   may miss a valid cross-batch merge if the root is outside maxDist of the
	 *   adjacent batch's root.  In practice this affects only a small fraction of
	 *   boundary variants and is accepted for the large speedup it enables.
	 *
	 * Complexity: O(batches × B·k·log k·log B  +  R·k'·log k'·log R)
	 *   vs  O(n·k·log k·log n)  for the monolithic edge loop,
	 *   where B = batch size, R = number of phase-2 representatives.
	 *
	 * Phase 1 batches are run in parallel using Settings.THREADS threads.
	 * This method is skipped when CLIQUE_MERGE or CENTROID_MERGE is set because
	 * those modes require global group membership during merging.
	 */
	void runHierarchicalMerging()
	{
		int batchSize = Settings.HIERARCHICAL_BATCH_SIZE;
		int numBatches = (n + batchSize - 1) / batchSize;
		long t0 = System.currentTimeMillis();
		System.out.printf("[Hierarchical] %,d variants → %d batches of ≤%,d%n",
			n, numBatches, batchSize);

		// ----------------------------------------------------------------
		// Phase 1: intra-batch merging (batches are independent → parallel)
		// batchRoot[i] = global index in data[] of the representative for i
		// ----------------------------------------------------------------
		int[] batchRoot = new int[n];
		int numThreads = Math.max(1, Settings.THREADS);
		ExecutorService pool = Executors.newFixedThreadPool(numThreads);
		CountDownLatch latch = new CountDownLatch(numBatches);

		for(int b = 0; b < numBatches; b++)
		{
			final int bStart = b * batchSize;
			final int bEnd   = Math.min(n, (b + 1) * batchSize);
			pool.submit(() -> {
				try
				{
					int len = bEnd - bStart;
					Variant[] batchData = new Variant[len];
					for(int i = 0; i < len; i++) batchData[i] = data[bStart + i];

					// Constructor sets batchData[i].index = i (local); runMerging uses it
					VariantMerger bvm = new VariantMerger(batchData);
					bvm.runMerging();

					// Map local roots → global indices; restore .index to global position
					// Different batches write to non-overlapping ranges so no synchronisation needed
					for(int i = 0; i < len; i++)
					{
						int localRoot = bvm.forest.find(i);
						batchRoot[bStart + i] = bStart + localRoot;
						batchData[i].index = bStart + i;
					}
				}
				finally { latch.countDown(); }
			});
		}
		try { latch.await(); } catch(InterruptedException ie) { Thread.currentThread().interrupt(); }
		pool.shutdown();

		int numReps = 0;
		for(int i = 0; i < n; i++) if(batchRoot[i] == i) numReps++;
		System.out.printf("[Hierarchical] Phase 1 done in %.1fs – %,d batch representatives%n",
			(System.currentTimeMillis() - t0) / 1000.0, numReps);

		// ---------------------------------------------------------------
		// Phase 2: inter-batch merge on the batch representatives
		// globalToRepPos[i] = position of global index i in reps list, or -1
		// repOrigIdx[j]     = global data[] index of the j-th representative
		// ---------------------------------------------------------------
		int[] globalToRepPos = new int[n];
		for(int i = 0; i < n; i++) globalToRepPos[i] = -1;
		int[] repOrigIdx = new int[numReps];
		ArrayList<Variant> reps = new ArrayList<>(numReps);

		for(int i = 0; i < n; i++)
		{
			if(batchRoot[i] == i)
			{
				int pos = reps.size();
				globalToRepPos[i] = pos;
				repOrigIdx[pos]   = i;
				reps.add(data[i]);
			}
		}

		long tP2 = System.currentTimeMillis();
		System.out.printf("[Hierarchical] Phase 2: merging %,d representatives%n", numReps);
		// Constructor overwrites .index for each rep to its position in reps[]
		VariantMerger finalVM = new VariantMerger(reps);

		// ---- Fix Phase-2 same-sample filtering ----
		// Forest(reps) initialises each rep's sampleSets entry with only its own
		// single sample.  But a batch representative may stand for a group that
		// already contains many samples (merged in Phase 1).  If we don't inject
		// the full sample membership here, Phase 2 will allow merges between two
		// representatives whose underlying batch groups share a sample — producing
		// an output variant with duplicate genotypes from that sample.
		//
		// Fix: replace each single-element sampleSets[j] with the union of all
		// samples that belong to that batch group.  batchRoot[i] already maps
		// every original variant i to its batch-group representative (global idx),
		// and globalToRepPos maps that global idx to the Phase-2 position j.
		if(!Settings.ALLOW_INTRASAMPLE && finalVM.forest.sampleSets != null)
		{
			// Wipe the single-element sets the constructor just built
			for(int j = 0; j < numReps; j++)
				finalVM.forest.sampleSets[j] = new java.util.HashSet<Integer>();
			// Accumulate the full sample membership for each batch group
			for(int i = 0; i < n; i++)
			{
				int repPos = globalToRepPos[batchRoot[i]];
				finalVM.forest.sampleSets[repPos].add(data[i].sample);
			}
		}

		finalVM.runMerging();
		// Restore .index for representatives back to global positions
		for(int j = 0; j < numReps; j++) reps.get(j).index = repOrigIdx[j];
		System.out.printf("[Hierarchical] Phase 2 done in %.1fs%n",
			(System.currentTimeMillis() - tP2) / 1000.0);

		// ---------------------------------------------------------------
		// Phase 3: propagate final groupings back into this.forest
		// for each original variant i:
		//   i → batchRoot[i] (batch rep, global) → globalToRepPos → repPos in reps[]
		//     → finalVM.forest.find(repPos) → finalRepPos → repOrigIdx[finalRepPos]
		//     = global root → union i with that root in this.forest
		// ---------------------------------------------------------------
		System.out.println("[Hierarchical] Phase 3: rebuilding union-find");
		for(int i = 0; i < n; i++)
		{
			int repPos      = globalToRepPos[batchRoot[i]];
			int finalRepPos = finalVM.forest.find(repPos);
			int finalGlobal = repOrigIdx[finalRepPos];
			if(i != finalGlobal)
			{
				int rootI     = forest.find(i);
				int rootFinal = forest.find(finalGlobal);
				if(rootI != rootFinal)
					forest.union(rootI, rootFinal);
			}
		}
		System.out.printf("[Hierarchical] Total time: %.1fs%n",
			(System.currentTimeMillis() - t0) / 1000.0);
	}

	/*
	 * Get an array of all of the groups of variants
	 */
	@SuppressWarnings("unchecked")
	ArrayList<Variant>[] getGroups()
	{
		ArrayList<Variant>[] res = new ArrayList[n];
		for(int i = 0; i<n; i++)
		{
			res[i] = new ArrayList<Variant>();
		}
		for(int i = 0; i<n; i++)
		{
			if(forest.map[i] < 0)
			{
				res[i].add(data[i]);
			}
			else
			{
				res[forest.find(i)].add(data[i]);
			}
		}
		return res;
	}
	
	/*
	 * An edge between two variants indicating that they can be merged
	 * Sorting is non-decreasing order of edge weights (ties broken with variant IDs)
	 */
	class Edge implements Comparable<Edge>
	{
		int from, to;
		double dist;
		Edge(int from, int to, double dist)
		{
			this.from = from;
			this.to = to;
			this.dist = dist;
		}
		@Override
		public int compareTo(Edge o) {
			if(Math.abs(dist - o.dist) > 1e-9)
			{
				return Double.compare(dist, o.dist);
			}
			if(data[from].hash != data[o.from].hash) return data[from].hash - (data[o.from].hash);
			if(data[to].hash != data[o.to].hash) return data[to].hash - (data[o.to].hash);
			if(from != o.from) return data[from].id.compareTo(data[o.from].id);
			return data[to].id.compareTo(data[o.to].id);
		}
	}
}
