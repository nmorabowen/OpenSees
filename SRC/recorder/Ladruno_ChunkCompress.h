/* ********************************************************************** **
**  Ladruno recorder — off-thread chunk compression.                        **
**                                                                        **
**  The compressed time-series layout (Ladruno_Hdf5.h, WP-233) stores ONE   **
**  STEP PER CHUNK, so each recorded step is a set of complete chunks.     **
**  A ChunkJob turns one chunk (one id tile of one step) into the exact    **
**  bytes HDF5's shuffle + deflate pipeline would store, on a worker       **
**  thread. The sink then hands the bytes to H5Dwrite_chunk on the solver  **
**  thread. Workers never call HDF5 (the library is not thread-safe        **
**  here); they only run the byte shuffle and zlib compress2.             **
** ********************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

#ifndef Ladruno_ChunkCompress_h
#define Ladruno_ChunkCompress_h

#include <atomic>
#include <cstddef>
#include <memory>
#include <vector>

namespace ladruno {
namespace chunkz {

	// One chunk to compress: rows [id0, id0 + chunk_ids) of a step slab of
	// n_ids x n_comp doubles. Ids past n_ids (the edge tile) are zero-padded:
	// HDF5 stores every chunk at full size.
	struct ChunkJob {
		std::shared_ptr<const std::vector<double> > slab;  // shared by a step's tiles
		size_t id0 = 0;
		size_t chunk_ids = 0;
		size_t n_ids = 0;
		size_t n_comp = 0;
		bool f32 = false;      // on-disk float32 (-precision f32): narrow first
		int level = 1;         // zlib level (the dataset's deflate filter level)

		std::vector<unsigned char> out;   // shuffle + deflate bytes, valid once done
		bool ok = false;
		std::atomic<int> state{0};        // 0 queued, 1 running, 2 done
	};

	// Compress `job` on the calling thread (also the synchronous fallback and the
	// at-exit path). Writes job.out / job.ok; does not touch job.state.
	void compressInto(const ChunkJob& job, std::vector<unsigned char>& out, bool& ok);

	// Queue a job on the process-wide pool. With threads == 0 (the pool could not
	// start), the job is compressed inline before this returns.
	void submit(const std::shared_ptr<ChunkJob>& job);

	// Block until `job` is done; if it is still queued, run it here instead of
	// waiting (keeps a blocking drain from idling behind a busy pool).
	void wait(const std::shared_ptr<ChunkJob>& job);

	// Size the pool (first call wins the minimum; later calls may only grow it).
	// n < 0 = automatic. Returns the number of workers running.
	int ensureThreads(int n);

	// The automatic worker count for this machine.
	int autoThreads();

} // namespace chunkz
} // namespace ladruno

#endif // Ladruno_ChunkCompress_h
