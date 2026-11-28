/* ********************************************************************** **
**  Ladruno recorder — off-thread chunk compression (see the header).      **
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

#include "Ladruno_ChunkCompress.h"

#include <zlib.h>

#include <algorithm>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <mutex>
#include <thread>

namespace ladruno {
namespace chunkz {

	void compressInto(const ChunkJob& job, std::vector<unsigned char>& out, bool& ok)
	{
		ok = false;
		out.clear();
		if (!job.slab || job.chunk_ids == 0 || job.n_comp == 0)
			return;
		const std::vector<double>& slab = *job.slab;
		const size_t esz = job.f32 ? 4 : 8;
		const size_t nelem = job.chunk_ids * job.n_comp;
		const size_t nbytes = nelem * esz;

		// 1. the chunk's raw little-endian bytes, edge ids zero-padded (the fill value)
		std::vector<unsigned char> raw(nbytes, 0);
		const size_t ids_here = (job.id0 < job.n_ids)
			? std::min(job.chunk_ids, job.n_ids - job.id0) : 0;
		const double* src = slab.data() + job.id0 * job.n_comp;
		const size_t nval = ids_here * job.n_comp;
		if (job.f32) {
			float* dst = reinterpret_cast<float*>(raw.data());
			for (size_t i = 0; i < nval; ++i) dst[i] = static_cast<float>(src[i]);
		}
		else if (nval) {
			std::memcpy(raw.data(), src, nval * sizeof(double));
		}

		// 2. HDF5 shuffle filter (H5Z_filter_shuffle): byte j of element i goes to
		//    j * nelem + i. nbytes is a whole number of elements, so no leftover.
		std::vector<unsigned char> shuf(nbytes);
		for (size_t j = 0; j < esz; ++j) {
			unsigned char* d = shuf.data() + j * nelem;
			const unsigned char* s = raw.data() + j;
			for (size_t i = 0; i < nelem; ++i)
				d[i] = s[i * esz];
		}

		// 3. HDF5 deflate filter (H5Z_filter_deflate): a zlib stream from compress2
		//    at the dataset's level; the reader inflates it whatever its size.
		uLongf zlen = compressBound((uLong)nbytes);
		out.resize(zlen);
		const int rc = compress2(out.data(), &zlen, shuf.data(), (uLong)nbytes, job.level);
		if (rc != Z_OK) {
			out.clear();
			return;
		}
		out.resize(zlen);
		ok = true;
	}

	namespace {

		// Process-wide pool. Deliberately leaked and its threads detached: a static
		// pool would be destroyed at exit with joinable threads (std::terminate), and
		// in the .pyd the CRT runs atexit after Windows has already killed them.
		struct Pool {
			std::mutex m;
			std::condition_variable work_cv;
			std::condition_variable done_cv;
			std::deque<std::shared_ptr<ChunkJob> > q;
			int threads = 0;

			void run(const std::shared_ptr<ChunkJob>& job)
			{
				int expect = 0;
				if (!job->state.compare_exchange_strong(expect, 1))
					return;   // a waiter claimed it and runs it itself
				bool ok = false;
				std::vector<unsigned char> out;
				compressInto(*job, out, ok);
				job->out.swap(out);
				job->ok = ok;
				{
					std::lock_guard<std::mutex> lk(m);
					job->state.store(2);
				}
				done_cv.notify_all();
			}

			void worker()
			{
				for (;;) {
					std::shared_ptr<ChunkJob> job;
					{
						std::unique_lock<std::mutex> lk(m);
						work_cv.wait(lk, [this] { return !q.empty(); });
						job = q.front();
						q.pop_front();
					}
					run(job);
				}
			}
		};

		Pool& pool()
		{
			static Pool* p = new Pool();   // never deleted (see above)
			return *p;
		}

	} // namespace

	int autoThreads()
	{
		const unsigned hw = std::thread::hardware_concurrency();
		int n = (hw >= 8) ? (int)(hw / 4) : 1;
		return std::max(1, std::min(n, 4));
	}

	int ensureThreads(int n)
	{
		if (n < 0) n = autoThreads();
		Pool& p = pool();
		std::lock_guard<std::mutex> lk(p.m);
		while (p.threads < n) {
			try {
				std::thread t([&p] { p.worker(); });
				t.detach();
				++p.threads;
			}
			catch (...) {
				break;   // no more threads: the jobs run inline (submit / wait)
			}
		}
		return p.threads;
	}

	void submit(const std::shared_ptr<ChunkJob>& job)
	{
		Pool& p = pool();
		{
			std::lock_guard<std::mutex> lk(p.m);
			if (p.threads > 0) {
				p.q.push_back(job);
				p.work_cv.notify_one();
				return;
			}
		}
		p.run(job);
	}

	void wait(const std::shared_ptr<ChunkJob>& job)
	{
		Pool& p = pool();
		if (job->state.load() == 0)
			p.run(job);   // still queued: claim it (no-op if a worker beat us)
		std::unique_lock<std::mutex> lk(p.m);
		p.done_cv.wait(lk, [&job] { return job->state.load() == 2; });
	}

} // namespace chunkz
} // namespace ladruno
