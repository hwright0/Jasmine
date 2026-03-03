/*
 * A representation of a forest using a union-find data structure
 * It allows nodes to be merged, checking if their components share
 * any variants from the same sample.
 *
 * Sample membership is tracked with a sparse HashSet<Integer> per root node
 * rather than a dense bitset.  This is critical when the number of samples is
 * large (e.g. 500 000) because the dense bitset would require
 * ceil(samples/63) longs × n variants of storage even though
 * each component typically contains only a tiny fraction of all samples.
 *
 * Union-by-rank (size) is preserved so that the small-to-large addAll trick
 * keeps the total merge cost O(n log n) across all unions.
 */

import java.util.Arrays;
import java.util.HashSet;

public class Forest
{
	int[] map; // map[i] is negative if root, more negative means bigger set; if nonneg, it is the parent index

	// Sparse per-root sample set. sampleSets[i] is non-null only for root nodes.
	// Null when ALLOW_INTRASAMPLE is true (no tracking needed).
	HashSet<Integer>[] sampleSets;

	@SuppressWarnings("unchecked")
	public Forest(Variant[] data)
	{
		int n = data.length;
		map = new int[n];
		Arrays.fill(map, -1);

		if(!Settings.ALLOW_INTRASAMPLE)
		{
			sampleSets = new HashSet[n];
			for(int i = 0; i < n; i++)
			{
				sampleSets[i] = new HashSet<Integer>();
				sampleSets[i].add(data[i].sample);
			}
		}
	}

	/*
	 * Get the root of the component containing a variant
	 */
	public int find(int x)
	{
		if(map[x] < 0)
			return x;
		else
		{
			map[x] = find(map[x]);
			return map[x];
		}
	}

	/*
	 * Check whether the two variants can be unioned without creating an intra-sample merge
	 */
	public boolean canUnion(int a, int b)
	{
		int roota = find(a), rootb = find(b);
		if(roota == rootb)
		{
			return false;
		}
		if(!okayEdge(roota, rootb))
		{
			return false;
		}
		return true;
	}

	public void union(int a, int b)
	{
		int roota = find(a), rootb = find(b);
		if(map[roota] < map[rootb])
		{
			// roota's component is larger (more negative = bigger); make roota the root
			map[roota] += map[rootb];
			map[rootb] = roota;
			if(!Settings.ALLOW_INTRASAMPLE)
			{
				sampleSets[roota].addAll(sampleSets[rootb]);
				sampleSets[rootb] = null; // allow GC to reclaim the non-root's set
			}
		}
		else
		{
			map[rootb] += map[roota];
			map[roota] = rootb;
			if(!Settings.ALLOW_INTRASAMPLE)
			{
				sampleSets[rootb].addAll(sampleSets[roota]);
				sampleSets[roota] = null; // allow GC to reclaim the non-root's set
			}
		}
	}

	/*
	 * Whether adding an edge between two root nodes would create an intra-sample merge.
	 * Iterates the smaller set and probes the larger for O(min(|A|,|B|)) performance.
	 */
	private boolean okayEdge(int rootA, int rootB)
	{
		if(Settings.ALLOW_INTRASAMPLE)
		{
			return true;
		}
		HashSet<Integer> smaller = sampleSets[rootA];
		HashSet<Integer> larger  = sampleSets[rootB];
		if(smaller.size() > larger.size())
		{
			smaller = sampleSets[rootB];
			larger  = sampleSets[rootA];
		}
		for(int s : smaller)
		{
			if(larger.contains(s))
			{
				return false;
			}
		}
		return true;
	}
}

