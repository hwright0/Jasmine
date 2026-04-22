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
		runMerging(Settings.THREADS, true);
	}

	/*
	 * Runs merging with an explicit thread count, skipping hierarchical logic.
	 * Used by Phase 1 batch sub-mergers to avoid recursive re-entry.
	 */
	void runMerging(int threads)
	{
		runMerging(threads, false);
	}

	/*
	 * Runs merging with an explicit thread count.
	 * allowHierarchical controls whether sample-based hierarchical merging
	 * is considered; Phase 1 sub-mergers pass false to prevent infinite recursion.
	 */
	void runMerging(int threads, boolean allowHierarchical)
	{
		if(n == 1)
		{
			return;
		}

		// Delegate to hierarchical merging for large dense graphs when requested.
		if(allowHierarchical && !Settings.CLIQUE_MERGE && !Settings.CENTROID_MERGE)
		{
			if(Settings.HIERARCHICAL_BATCH_SAMPLES > 0)
			{
				runSampleHierarchicalMerging();
				return;
			}
		}

		long mergeStartTime = System.currentTimeMillis();
		int numThreads = Math.max(1, threads);
		System.out.printf("[Merge] Graph with %,d variants: building kNN graph (%d threads)%n",
			n, numThreads);

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
	 * Hierarchical merging partitioned by sample ID.
	 *
	 * Preferred over position-based batching for large cohorts because:
	 *   1. Every batch covers ALL genomic loci, so variants at the same position
	 *      always compete within one batch.  No positional boundary artefacts.
	 *   2. Batches are disjoint sample sets by construction, so Phase 2 can never
	 *      produce a merged variant with two genotypes from the same person — no
	 *      same-sample injection fix is needed.
	 *   3. Phase 2 is tiny: (numBatches × distinctLociPerBatch) representatives,
	 *      e.g. 50 batches × 20 loci = 1000 reps for a 500k-sample cohort.
	 *
	 * Phase 1: each batch of Settings.HIERARCHICAL_BATCH_SAMPLES consecutive sample
	 *          IDs is merged independently, in parallel.
	 * Phase 2: one Phase-2 VariantMerger runs over all batch representatives.
	 * Phase 3: final groupings are propagated back into this.forest.
	 *
	 * Complexity: O(numBatches × B_var·k·log k·log B_var  +  R·k'·log k'·log R)
	 *   where B_var = variants in one sample batch (≈ n/numBatches),
	 *         R     = Phase-2 representatives (≈ numBatches × numDistinctLoci).
	 */
	void runSampleHierarchicalMerging()
	{
		int samplesPerBatch = Settings.HIERARCHICAL_BATCH_SAMPLES;
		long t0 = System.currentTimeMillis();

		// Find the range of sample IDs present in this graph
		int maxSample = 0;
		for(int i = 0; i < n; i++)
			if(data[i].sample > maxSample) maxSample = data[i].sample;

		int numBatches = (maxSample / samplesPerBatch) + 1;

		// Partition variant indices by sample batch.
		// Each bucket b holds all data[] index positions where data[i].sample
		// falls in [ b*samplesPerBatch, (b+1)*samplesPerBatch ).
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] buckets = new ArrayList[numBatches];
		for(int b = 0; b < numBatches; b++) buckets[b] = new ArrayList<>();
		for(int i = 0; i < n; i++) buckets[data[i].sample / samplesPerBatch].add(i);

		// Drop empty buckets so the thread pool isn't flooded with no-ops
		ArrayList<ArrayList<Integer>> activeBuckets = new ArrayList<>();
		for(ArrayList<Integer> bucket : buckets)
			if(!bucket.isEmpty()) activeBuckets.add(bucket);
		int numActive = activeBuckets.size();

		System.out.printf("[SampleHierarchical] %,d variants, %d active sample batches " +
			"(≤%d samples each), %d threads%n",
			n, numActive, samplesPerBatch, Math.max(1, Settings.THREADS));

		// ----------------------------------------------------------------
		// Phase 1: intra-batch merging — each batch is independent → parallel
		// batchRoot[i] = global data[] index of i's batch-group representative
		// ----------------------------------------------------------------
		int[] batchRoot = new int[n];
		for(int i = 0; i < n; i++) batchRoot[i] = i; // default: each variant is its own root

		int numThreads = Math.max(1, Settings.THREADS);
		ExecutorService pool = Executors.newFixedThreadPool(numThreads);
		CountDownLatch latch = new CountDownLatch(numActive);

		for(int b = 0; b < numActive; b++)
		{
			final ArrayList<Integer> bucket = activeBuckets.get(b);
			pool.submit(() -> {
				try
				{
					int len = bucket.size();
					Variant[] batchData = new Variant[len];
					for(int i = 0; i < len; i++) batchData[i] = data[bucket.get(i)];

					VariantMerger bvm = new VariantMerger(batchData);
					bvm.runMerging(1); // single-threaded: batches already run in parallel

					// Map local roots back to global data[] indices;
					// restore .index to global positions.
					// Each bucket covers a non-overlapping range of sample IDs so
					// no two threads write to the same batchRoot[] slot.
					for(int i = 0; i < len; i++)
					{
						int localRoot     = bvm.forest.find(i);
						int globalIdx     = bucket.get(i);
						int globalRootIdx = bucket.get(localRoot);
						batchRoot[globalIdx]    = globalRootIdx;
						batchData[i].index = globalIdx;
					}
				}
				finally { latch.countDown(); }
			});
		}
		try { latch.await(); } catch(InterruptedException ie) { Thread.currentThread().interrupt(); }
		pool.shutdown();

		int numReps = 0;
		for(int i = 0; i < n; i++) if(batchRoot[i] == i) numReps++;
		System.out.printf("[SampleHierarchical] Phase 1 done in %.1fs – %,d representatives%n",
			(System.currentTimeMillis() - t0) / 1000.0, numReps);

		// ---------------------------------------------------------------
		// Phase 2: merge the inter-batch representatives
		// globalToRepPos[i] = position of global index i in reps[], or -1
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
		System.out.printf("[SampleHierarchical] Phase 2: merging %,d representatives%n", numReps);

		VariantMerger finalVM = new VariantMerger(reps);

		// Inject accumulated sample sets from Phase 1 into Phase 2's forest.
		// Each Phase 1 representative may stand for a group spanning many samples.
		// Without this injection, Phase 2's canUnion would only see the single
		// sample of the representative itself and could incorrectly merge groups
		// that share a sample.
		if(!Settings.ALLOW_INTRASAMPLE && finalVM.forest.sampleSets != null)
		{
			for(int j = 0; j < numReps; j++)
				finalVM.forest.sampleSets[j] = new java.util.HashSet<Integer>();
			for(int i = 0; i < n; i++)
			{
				int repPos = globalToRepPos[batchRoot[i]];
				finalVM.forest.sampleSets[repPos].add(data[i].sample);
			}
		}

		finalVM.runMerging(Settings.THREADS);  // allowHierarchical=false to prevent re-entry
		for(int j = 0; j < numReps; j++) reps.get(j).index = repOrigIdx[j];
		System.out.printf("[SampleHierarchical] Phase 2 done in %.1fs%n",
			(System.currentTimeMillis() - tP2) / 1000.0);

		// ---------------------------------------------------------------
		// Phase 3: propagate final groupings back into this.forest
		// ---------------------------------------------------------------
		System.out.println("[SampleHierarchical] Phase 3: rebuilding union-find");
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
		System.out.printf("[SampleHierarchical] Total time: %.1fs%n",
			(System.currentTimeMillis() - t0) / 1000.0);
	}

	/*
	 * Get an array of all of the groups of variants.
	 * Returns a compact array sized to the number of actual merged groups
	 * (not the total variant count), so downstream data structures stay small.
	 */
	@SuppressWarnings("unchecked")
	ArrayList<Variant>[] getGroups()
	{
		// First pass: assign a compact index to each unique root
		int[] rootToCompact = new int[n];
		java.util.Arrays.fill(rootToCompact, -1);
		int numGroups = 0;
		for(int i = 0; i < n; i++)
		{
			int root = forest.find(i);
			if(rootToCompact[root] == -1)
			{
				rootToCompact[root] = numGroups++;
			}
		}
		// Allocate only as many lists as there are groups
		ArrayList<Variant>[] res = new ArrayList[numGroups];
		for(int i = 0; i < numGroups; i++) res[i] = new ArrayList<Variant>();
		// Second pass: place each variant into its compact group
		for(int i = 0; i < n; i++)
		{
			res[rootToCompact[forest.find(i)]].add(data[i]);
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
