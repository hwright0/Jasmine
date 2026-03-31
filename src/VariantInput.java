/*
 * Methods for reading VCF entries from VCF files and dividing
 * the entries into separate groups by graph ID
 */

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class VariantInput {
	
	// How many samples were merged to produce each input file
	// ConcurrentHashMap so that parallel file reading in readAllFiles is safe
	static ConcurrentHashMap<Integer, Integer> previouslyMergedSamples = new ConcurrentHashMap<Integer, Integer>();
	
	/*
	 * Count the number of VCF files in a list
	 */
	public static int countFiles(String fileList) throws Exception
	{
		return PipelineManager.getFilesFromList(fileList).size();
	}
	
	/*
	 * Get a list of all variants from a group of files, binning them by graphID
	 */
	@SuppressWarnings("unchecked")
	public static TreeMap<String, ArrayList<Variant>> readAllFiles(String fileList) throws Exception
	{
		ArrayList<String> fileNames = PipelineManager.getFilesFromList(fileList);
		int numFiles = fileNames.size();
		
		TreeMap<String, ArrayList<Variant>>[] variantsPerFile = new TreeMap[numFiles];
		
		// Read all input VCF files in parallel using the configured thread count.
		// Each getSingleList call is independent (distinct sample index, no shared write state).
		int threads = Math.max(1, Math.min(Settings.THREADS, numFiles));
		ExecutorService pool = Executors.newFixedThreadPool(threads);
		ArrayList<Future<?>> futures = new ArrayList<Future<?>>();
		for(int i = 0; i < numFiles; i++)
		{
			final int idx = i;
			futures.add(pool.submit(() -> {
				try
				{
					variantsPerFile[idx] = getSingleList(fileNames.get(idx), idx);
				}
				catch(Exception e)
				{
					throw new RuntimeException(e);
				}
			}));
		}
		pool.shutdown();
		for(Future<?> f : futures)
		{
			f.get(); // re-throws ExecutionException wrapping any per-file error
		}
		
		// Merge per-file maps into the combined result (single-threaded; order preserved)
		TreeMap<String, ArrayList<Variant>> res = new TreeMap<String, ArrayList<Variant>>();
		for(int i = 0; i < numFiles; i++)
		{
			for(String s : variantsPerFile[i].keySet())
			{
				if(!res.containsKey(s))
				{
					res.put(s, new ArrayList<Variant>());
				}
				res.get(s).addAll(variantsPerFile[i].get(s));
			}
		}
		return res;
	}
	
	/*
	 * Get a list of variants binned by graphID for a single VCF file
	 */
	private static TreeMap<String, ArrayList<Variant>> getSingleList(String filename, int sample) throws Exception
	{
		BgzipReader bgzReader = null;
		Scanner input;
		if(filename.endsWith(".gz"))
		{
			bgzReader = new BgzipReader(filename);
			input = bgzReader.getScanner();
		}
		else
		{
			input = new Scanner(new BufferedInputStream(new FileInputStream(new File(filename))));
		}
		ArrayList<Variant> allVariants = new ArrayList<Variant>();
		HashSet<String> ids = new HashSet<String>();
		if(!previouslyMergedSamples.containsKey(sample))
		{
			previouslyMergedSamples.put(sample, 1);
		}
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			if(line.length() >=2 && line.charAt(0) == 31 && (line.charAt(1) == 65533 || line.charAt(1) == 139))
			{
				throw new Exception(filename + " appears to be a gzipped file but does not have a .gz extension");
			}
			VcfEntry entry = VcfEntry.fromLine(line);
			if(!previouslyMergedSamples.containsKey(sample))
			{
				if(entry.getInfo("SUPP_VEC_EXT").length() > 0)
				{
					previouslyMergedSamples.put(sample, entry.getInfo("SUPP_VEC_EXT").length());
				}
				else if(entry.getInfo("SUPP_VEC").length() > 0)
				{
					previouslyMergedSamples.put(sample, entry.getInfo("SUPP_VEC").length());
				}
				else
				{
					previouslyMergedSamples.put(sample, 1);
				}
			}
			if(ids.contains(entry.getId()))
			{
				String oldId = entry.getId();
				int index = 1;
				while(true)
				{
					String newId = oldId + "_duplicate" + index;
					if(!ids.contains(newId))
					{
						entry.setId(newId);
						break;
					}
					else
					{
						index++;
					}
				}
				System.err.println("Warning: Duplicate variant ID " + oldId + " in " + filename + "; Replacing with " + entry.getId());
			}
			ids.add(entry.getId());
			allVariants.add(fromVcfEntry(entry, sample));
			
		}
		
		System.out.println(filename + " has " + allVariants.size() + " variants");
		if(bgzReader != null)
		{
			bgzReader.close();
		}
		else
		{
			input.close();
		}
		
		return divideIntoGraphs(allVariants);
	}
	
	/*
	 * Take a list of variants and bin them by graphID
	 */
	private static TreeMap<String, ArrayList<Variant>> divideIntoGraphs(ArrayList<Variant> data)
	{
		TreeMap<String, ArrayList<Variant>> groups = new TreeMap<String, ArrayList<Variant>>();
		for(Variant v : data)
		{
			String graphID = v.graphID;
			if(!groups.containsKey(graphID))
			{
				groups.put(graphID, new ArrayList<Variant>());
			}
			groups.get(graphID).add(v);
		}
		return groups;
	}
	
	/*
	 * From a line of a VCF file, extract the information needed for merging
	 * and return it as a Variant object
	 */
	public static Variant fromVcfEntry(VcfEntry entry, int sample) throws Exception
	{
		double start = entry.getFirstCoord();
		double end = entry.getSecondCoord();
		
		entry.setId(sample + "_" + entry.getId());
		
		String id = entry.getGraphID();
		
		String seq = null;
		if(entry.getType().equals("INS"))
		{
			String entrySeq = entry.getSeq();
			if(entrySeq.length() > 0)
			{
				seq = entrySeq;
			}
		}
		
		// Default distance threshold model is constant, so set to that first
		int maxDist = Settings.MAX_DIST;
		double minSeqId = Settings.MIN_SEQUENCE_SIMILARITY;
		
		// Then, check if there is a per-variant distance threshold
		String maxDistInfo = entry.getInfo("JASMINE_DIST");
		if(maxDistInfo.length() > 0)
		{
			maxDist = Integer.parseInt(maxDistInfo);
		}
		
		// Next check if a per-sample distance threshold was set
		else if(Settings.PER_SAMPLE_DISTS != null && Settings.PER_SAMPLE_DISTS.length > sample)
		{
			maxDist = Settings.PER_SAMPLE_DISTS[sample];
		}
		
		// Next, check if there is a length-based threshold
		else if(Settings.USE_LINEAR_THRESHOLD && Settings.MAX_DIST_LINEAR > 0)
		{
			maxDist = (int)(Settings.MAX_DIST_LINEAR * Math.abs(entry.getLength()) + 0.5);
			if(Settings.MAX_DIST_SET)
			{
				maxDist = Math.min(maxDist, Settings.MAX_DIST);
			}
			if(Settings.MIN_DIST != -1)
			{
				maxDist = Math.max(maxDist, Settings.MIN_DIST);
			}
		}
		
		// Check for per-variant sequence ID thresholds
		String minIdInfo = entry.getInfo("JASMINE_ID");
		if(minIdInfo.length() > 0)
		{
			minSeqId = Double.parseDouble(minIdInfo);
		}
		
		Variant res = new Variant(sample, entry.getId(), start, end, id, seq, maxDist, minSeqId);
		res.hash = Variant.hash(entry.tabTokens[7]);
		if(Settings.OVERLAP_REQUIRED > 0 && (entry.getType().equals("DEL")) || entry.getType().equals("INV") || entry.getType().equals("DUP"))
		{
			res.interval = new double[] {entry.getPos(), entry.getEnd()};
		}
		return res;
	}
}
