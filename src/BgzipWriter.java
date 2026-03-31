/*
 * Pure-Java multi-threaded BGZF writer.
 * Buffers data into 64 KB blocks and compresses them in parallel using
 * java.util.zip.Deflater, producing spec-compliant BGZF output that is
 * compatible with tabix, samtools, and other htslib-based tools.
 * No external bgzip binary is required.
 */

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

public class BgzipWriter implements Closeable {

	// Maximum uncompressed payload per BGZF block (leaves room for header+trailer)
	private static final int MAX_BLOCK = 65280;

	// 28-byte empty BGZF block that marks end-of-file
	private static final byte[] EOF_BLOCK = {
		31, (byte)139, 8, 4, 0, 0, 0, 0, 0, (byte)255,
		6, 0, 66, 67, 2, 0, 27, 0,
		3, 0,
		0, 0, 0, 0, 0, 0, 0, 0
	};

	private final OutputStream rawOut;
	private final PrintWriter writer;
	private final ExecutorService pool;
	private final LinkedList<Future<byte[]>> pending = new LinkedList<>();
	private byte[] uncompressed = new byte[MAX_BLOCK];
	private int pos = 0;
	private boolean closed = false;
	private final int maxPending;

	BgzipWriter(String outputPath, int threads) throws IOException
	{
		rawOut = new BufferedOutputStream(new FileOutputStream(outputPath));
		int t = Math.max(1, threads);
		pool = Executors.newFixedThreadPool(t);
		maxPending = t * 2;
		writer = new PrintWriter(new BgzfStream(), false);
	}

	PrintWriter getWriter()
	{
		return writer;
	}

	@Override
	public void close() throws IOException
	{
		if(closed) return;
		closed = true;
		writer.flush();
		if(pos > 0) submitBlock();
		drainAll();
		rawOut.write(EOF_BLOCK);
		rawOut.close();
		pool.shutdown();
	}

	/*
	 * Submit the current buffer as a compression task and reset the buffer
	 */
	private void submitBlock() throws IOException
	{
		byte[] data = Arrays.copyOf(uncompressed, pos);
		pos = 0;
		pending.add(pool.submit(() -> compressBgzfBlock(data)));
		// Block if too many pending to bound memory usage
		while(pending.size() > maxPending)
		{
			try { rawOut.write(pending.poll().get()); }
			catch(Exception e) { throw new IOException("BGZF compression error", e); }
		}
		// Drain completed futures from the front (non-blocking)
		while(!pending.isEmpty() && pending.peek().isDone())
		{
			try { rawOut.write(pending.poll().get()); }
			catch(Exception e) { throw new IOException("BGZF compression error", e); }
		}
	}

	/*
	 * Block until all pending compression tasks are written
	 */
	private void drainAll() throws IOException
	{
		while(!pending.isEmpty())
		{
			try { rawOut.write(pending.poll().get()); }
			catch(Exception e) { throw new IOException("BGZF compression error", e); }
		}
	}

	/*
	 * Compress a single uncompressed block into a complete BGZF block
	 */
	static byte[] compressBgzfBlock(byte[] data) throws IOException
	{
		Deflater deflater = new Deflater(Deflater.DEFAULT_COMPRESSION, true); // true = raw deflate (no zlib header)
		deflater.setInput(data);
		deflater.finish();
		byte[] compBuf = new byte[data.length + 1024];
		int compLen = 0;
		while(!deflater.finished())
		{
			compLen += deflater.deflate(compBuf, compLen, compBuf.length - compLen);
		}
		deflater.end();

		CRC32 crc = new CRC32();
		crc.update(data);

		// Block layout: header(18) + compressed data + CRC32(4) + ISIZE(4)
		int blockSize = 18 + compLen + 8;
		byte[] block = new byte[blockSize];

		// Standard gzip header with FEXTRA flag
		block[0] = 31;           // ID1
		block[1] = (byte)139;    // ID2
		block[2] = 8;            // CM = deflate
		block[3] = 4;            // FLG = FEXTRA
		// MTIME(4), XFL, OS
		block[9] = (byte)255;    // OS = unknown
		block[10] = 6;           // XLEN lo
		block[11] = 0;           // XLEN hi
		// BGZF extra subfield: SI1='B', SI2='C', SLEN=2, BSIZE
		block[12] = 66;          // SI1 = 'B'
		block[13] = 67;          // SI2 = 'C'
		block[14] = 2;           // SLEN lo
		block[15] = 0;           // SLEN hi
		int bsize = blockSize - 1;
		block[16] = (byte)(bsize & 0xFF);
		block[17] = (byte)((bsize >> 8) & 0xFF);

		// Compressed payload
		System.arraycopy(compBuf, 0, block, 18, compLen);

		// CRC32 (little-endian)
		int c = (int)crc.getValue();
		int off = 18 + compLen;
		block[off]     = (byte)(c & 0xFF);
		block[off + 1] = (byte)((c >> 8) & 0xFF);
		block[off + 2] = (byte)((c >> 16) & 0xFF);
		block[off + 3] = (byte)((c >> 24) & 0xFF);

		// ISIZE = uncompressed size (little-endian)
		int s = data.length;
		block[off + 4] = (byte)(s & 0xFF);
		block[off + 5] = (byte)((s >> 8) & 0xFF);
		block[off + 6] = (byte)((s >> 16) & 0xFF);
		block[off + 7] = (byte)((s >> 24) & 0xFF);

		return block;
	}

	/*
	 * Internal OutputStream that buffers writes and submits full blocks for compression
	 */
	private class BgzfStream extends OutputStream
	{
		@Override
		public void write(int b) throws IOException
		{
			uncompressed[pos++] = (byte)b;
			if(pos >= MAX_BLOCK) submitBlock();
		}

		@Override
		public void write(byte[] data, int off, int len) throws IOException
		{
			while(len > 0)
			{
				int space = MAX_BLOCK - pos;
				int n = Math.min(space, len);
				System.arraycopy(data, off, uncompressed, pos, n);
				pos += n;
				off += n;
				len -= n;
				if(pos >= MAX_BLOCK) submitBlock();
			}
		}

		@Override
		public void flush() throws IOException
		{
			if(pos > 0) submitBlock();
			drainAll();
			rawOut.flush();
		}

		@Override
		public void close() { /* handled by BgzipWriter.close() */ }
	}
}
