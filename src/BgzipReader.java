/*
 * Pure-Java BGZF reader.
 * Reads BGZF block headers to determine block sizes, decompresses each block
 * using java.util.zip.Inflater, and presents the result as a Scanner.
 * No external bgzip binary is required.
 */

import java.io.BufferedInputStream;
import java.io.Closeable;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Scanner;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

public class BgzipReader implements Closeable {
	private final BgzfInputStream bgzfIn;
	private final Scanner scanner;

	BgzipReader(String inputPath) throws IOException
	{
		bgzfIn = new BgzfInputStream(new BufferedInputStream(new FileInputStream(inputPath)));
		scanner = new Scanner(new BufferedInputStream(bgzfIn));
	}

	Scanner getScanner()
	{
		return scanner;
	}

	@Override
	public void close() throws IOException
	{
		scanner.close();
		bgzfIn.close();
	}

	/*
	 * InputStream that reads and decompresses BGZF blocks one at a time
	 */
	private static class BgzfInputStream extends InputStream
	{
		private final DataInputStream in;
		private byte[] buf = new byte[65536];
		private int pos = 0;
		private int len = 0;
		private boolean eof = false;

		BgzfInputStream(InputStream in)
		{
			this.in = new DataInputStream(in);
		}

		/*
		 * Read the next BGZF block and decompress it into buf.
		 * Returns false if EOF is reached.
		 */
		private boolean nextBlock() throws IOException
		{
			// Read gzip magic bytes
			int id1 = in.read();
			if(id1 == -1) { eof = true; return false; }
			int id2 = in.read();
			if(id2 == -1 || id1 != 31 || id2 != 139)
				throw new IOException("Invalid BGZF block header");

			// CM, FLG, MTIME(4), XFL, OS = 8 bytes
			byte[] hdr = new byte[8];
			in.readFully(hdr);

			// XLEN (2 bytes, little-endian)
			int xlo = in.readUnsignedByte();
			int xhi = in.readUnsignedByte();
			int xlen = xlo | (xhi << 8);

			// Read extra field
			byte[] extra = new byte[xlen];
			in.readFully(extra);

			// Find BC subfield to get BSIZE
			int bsize = -1;
			int i = 0;
			while(i + 5 <= xlen)
			{
				int si1 = extra[i] & 0xFF;
				int si2 = extra[i + 1] & 0xFF;
				int slen = (extra[i + 2] & 0xFF) | ((extra[i + 3] & 0xFF) << 8);
				if(si1 == 66 && si2 == 67 && slen == 2)
				{
					bsize = (extra[i + 4] & 0xFF) | ((extra[i + 5] & 0xFF) << 8);
					break;
				}
				i += 4 + slen;
			}
			if(bsize == -1) throw new IOException("BGZF BC subfield not found");

			// Read the rest of the block: compressed data + CRC32(4) + ISIZE(4)
			int headerRead = 12 + xlen; // 10 (fixed gzip header) + 2 (XLEN field) + xlen
			int remaining = bsize + 1 - headerRead;
			byte[] rest = new byte[remaining];
			in.readFully(rest);

			// ISIZE is the last 4 bytes (little-endian)
			int isize = (rest[remaining - 4] & 0xFF)
			          | ((rest[remaining - 3] & 0xFF) << 8)
			          | ((rest[remaining - 2] & 0xFF) << 16)
			          | ((rest[remaining - 1] & 0xFF) << 24);

			// Empty block = EOF marker
			if(isize == 0) { eof = true; return false; }

			// Decompress using raw inflate (no zlib header)
			int compLen = remaining - 8;
			Inflater inf = new Inflater(true);
			inf.setInput(rest, 0, compLen);
			try
			{
				int total = 0;
				while(!inf.finished())
				{
					total += inf.inflate(buf, total, buf.length - total);
				}
				len = total;
			}
			catch(DataFormatException e)
			{
				throw new IOException("BGZF decompression error", e);
			}
			finally
			{
				inf.end();
			}
			pos = 0;
			return len > 0;
		}

		@Override
		public int read() throws IOException
		{
			if(pos >= len && (eof || !nextBlock())) return -1;
			return buf[pos++] & 0xFF;
		}

		@Override
		public int read(byte[] b, int off, int readLen) throws IOException
		{
			if(pos >= len && (eof || !nextBlock())) return -1;
			int avail = len - pos;
			int n = Math.min(avail, readLen);
			System.arraycopy(buf, pos, b, off, n);
			pos += n;
			return n;
		}

		@Override
		public void close() throws IOException
		{
			in.close();
		}
	}
}
