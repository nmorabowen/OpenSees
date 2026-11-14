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

#ifndef Ladruno_Hdf5_h
#define Ladruno_Hdf5_h

/*************************************************************************************

Ladruno_Hdf5.h

HDF5 interface layer (namespace h5) ported VERBATIM from MPCORecorder.cpp
(frozen). Every function is made `inline` because, unlike the frozen .cpp where
these are plain TU-local functions, here they live in a header included by
multiple TUs and must avoid duplicate-symbol errors. Everything is wrapped in
`namespace ladruno { ... }` to avoid ODR clashes with the frozen file's identical
`namespace h5` symbols.

Context: in this build the macro LADRUNO_HDF5_LOADED_AT_RUNTIME is NEVER defined,
so the real linked HDF5 C API is used directly (H5*). No runtime loader /
ptr_H5* machinery is involved.

**************************************************************************************/

#include "Ladruno_Types.h"

namespace ladruno {

/*ladruno - hdf5 interface utilities*/
namespace h5 {

	namespace attribute {

		// low level functions for c interface

		inline herr_t writei(hid_t obj, const char *attr_name, const int *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t attr = H5Acreate(obj, attr_name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, H5T_NATIVE_INT, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			return status;
		}
		inline herr_t writed(hid_t obj, const char *attr_name, const double *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t attr = H5Acreate(obj, attr_name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, H5T_NATIVE_DOUBLE, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			return status;
		}
		inline herr_t writes(hid_t obj, const char *attr_name, const char *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { 1 };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t atype = H5Tcopy(H5T_C_S1);
			status = H5Tset_size(atype, data_size);
			status = H5Tset_strpad(atype, H5T_STR_NULLTERM);
			hid_t attr = H5Acreate(obj, attr_name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, atype, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			H5Tclose(atype);  // Ladruno WP-163 (ARCH-4): was leaked per call
			return status;
		}

		// higher level c++ utils

		inline herr_t write(hid_t obj, const char *attr_name, int attr_data)
		{
			return writei(obj, attr_name, &attr_data, 1);
		}
		inline herr_t write(hid_t obj, const char *attr_name, const std::vector<int> &attr_data)
		{
			if (attr_data.size() > 0) {
				return writei(obj, attr_name, &attr_data[0], attr_data.size());
			}
			return 0;
		}
		inline herr_t write(hid_t obj, const char *attr_name, double attr_data)
		{
			return writed(obj, attr_name, &attr_data, 1);
		}
		inline herr_t write(hid_t obj, const char *attr_name, std::vector<double> &attr_data)
		{
			if (attr_data.size() > 0) {
				return writed(obj, attr_name, &attr_data[0], attr_data.size());
			}
			return 0;
		}
		inline herr_t write(hid_t obj, const char *attr_name, const std::string &attr_data)
		{
			if (!attr_data.empty()) {
				return writes(obj, attr_name, attr_data.c_str(), attr_data.size());
			}
			return 0;
		}

	}

	namespace group {

		// low level functions for c interface

		inline hid_t create(hid_t loc_id, const char *name, hid_t lcpl_id, hid_t gcpl_id, hid_t gapl_id) {
			return H5Gcreate(loc_id, name, lcpl_id, gcpl_id, gapl_id);
		}
		inline herr_t close(hid_t group_id) {
			if (group_id == HID_INVALID) return -1;
			return H5Gclose(group_id);
		}

		// higher level c++ utils

		inline hid_t createResultGroup(hid_t obj, hid_t gplist, const std::string &result_name,
			const std::string &disp_name, const std::string &components, int num_components,
			const std::string &dimension, const std::string &description,
			int result_type, int data_type)
		{
			herr_t status;
			hid_t h_gp_result = h5::group::create(obj, result_name.c_str(), H5P_DEFAULT, gplist, H5P_DEFAULT);
			status = h5::attribute::write(h_gp_result, "DISPLAY_NAME", disp_name);
			status = h5::attribute::write(h_gp_result, "COMPONENTS", components);
			status = h5::attribute::write(h_gp_result, "DIMENSION", dimension);
			status = h5::attribute::write(h_gp_result, "DESCRIPTION", description);
			status = h5::attribute::write(h_gp_result, "TYPE", result_type);
			status = h5::attribute::write(h_gp_result, "DATA_TYPE", data_type);
			return h_gp_result;
		}

	}

	namespace dataset {

		// low level functions for c interface

		inline herr_t close(hid_t dataset_id) {
			if (dataset_id == HID_INVALID) return -1;
			return H5Dclose(dataset_id);
		}
		inline hid_t createAndWrited1(hid_t obj, const char *name, const double *data, hsize_t data_size) {
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		inline hid_t createAndWrited2(hid_t obj, const char *name, const double *data, hsize_t rows, hsize_t cols)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[2] = { rows, cols };
			hid_t space = H5Screate_simple(2, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		inline hid_t createAndWritei1(hid_t obj, const char *name, const int *data, hsize_t data_size)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		inline hid_t createAndWritei2(hid_t obj, const char *name, const int *data, hsize_t rows, hsize_t cols)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[2] = { rows, cols };
			hid_t space = H5Screate_simple(2, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		inline hid_t createAndWrites(hid_t obj, const char *name, const char *data, hsize_t data_size)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { 1 };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t atype = H5Tcopy(H5T_C_S1);
			status = H5Tset_size(atype, data_size);
			status = H5Tset_strpad(atype, H5T_STR_NULLTERM);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			H5Tclose(atype);  // Ladruno WP-163 (ARCH-4): was leaked per call
			return dset;
		}

		// higher level c++ utils

		inline hid_t createAndWrite(hid_t obj, const char *name, double data)
		{
			return createAndWrited1(obj, name, &data, 1);
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<double> &data)
		{
			if (data.size() > 0) {
				return createAndWrited1(obj, name, &data[0], data.size());
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<double> &data, size_t rows, size_t cols)
		{
			if (data.size() > 0 && data.size() == rows*cols) {
				return createAndWrited2(obj, name, &data[0], rows, cols);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, int data)
		{
			return createAndWritei1(obj, name, &data, 1);
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<int> &data)
		{
			if (data.size() > 0) {
				return createAndWritei1(obj, name, &data[0], data.size());
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<int> &data, size_t rows, size_t cols)
		{
			if (data.size() > 0 && data.size() == rows*cols) {
				return createAndWritei2(obj, name, &data[0], rows, cols);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::string &data)
		{
			if (!data.empty()) {
				return createAndWrites(obj, name, data.c_str(), data.size());
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<utils::locax::quaternion> &data)
		{
			if (data.size() > 0) {
				return createAndWrited2(obj, name, data[0].data(), data.size(), 4);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<detail::element::FiberData> &data)
		{
			if (data.size() > 0) {
				return createAndWrited2(obj, name, data[0].data(), data.size(), 3);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<detail::element::ElemGaussPair> &data)
		{
			if (data.size() > 0) {
				return createAndWritei2(obj, name, data[0].data(), data.size(), 2);
			}
			return HID_INVALID;
		}

		// ---- extensible (chunked) time-series datasets (schema D3) -----------
		// One [T x nIds x nComp] dataset per result grows along the unlimited T
		// (time/mode) axis; one slab is appended per committed step. Chunked +
		// shuffle + deflate so the (compressible) time axis is compressed and a
		// single time-history is an O(1) dataset open instead of O(T) per-step
		// dataset opens. TIME[T] (double) / STEP[T] (int) are the matching
		// extensible 1-D axes. Replaces the per-step DATA/STEP_<k> layout.

		// choose how many T-slabs per chunk to land near ~256 KiB/chunk.
		inline hsize_t chunkSlabsPerBlock(hsize_t slab_elems, hsize_t elem_bytes) {
			const hsize_t target = 262144; // 256 KiB
			hsize_t per = slab_elems * elem_bytes;
			if (per == 0) return 1;
			hsize_t n = target / per;
			if (n < 1)    n = 1;
			if (n > 1024) n = 1024;
			return n;
		}

		// create a [0 x nIds x nComp] dataset, maxdims [UNLIM, nIds, nComp].
		// disk_type selects the ON-DISK element type (default f64 = lossless parity
		// path). Pass H5T_IEEE_F32LE for the opt-in `-precision f32` lossy mode:
		// appendSlab3d still writes from an in-memory double array and HDF5 narrows
		// f64->f32 on write, halving the payload at ~7 significant digits. The chunk
		// byte budget uses the on-disk element size so chunks stay ~256 KiB either way.
		inline hid_t createTimeSeries3d(hid_t obj, const char *name, hsize_t n_ids, hsize_t n_comp,
		                                hid_t disk_type = H5T_IEEE_F64LE) {
			if (n_ids < 1 || n_comp < 1) return HID_INVALID;  // WP-163 R5: unrepresentable
			hsize_t dims[3]    = { 0, n_ids, n_comp };
			hsize_t maxdims[3] = { H5S_UNLIMITED, n_ids, n_comp };
			hid_t space = H5Screate_simple(3, dims, maxdims);
			if (space < 0) return HID_INVALID;
			hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
			if (dcpl < 0) { H5Sclose(space); return HID_INVALID; }
			hsize_t elem_bytes = (H5Tequal(disk_type, H5T_IEEE_F32LE) > 0) ? 4 : 8;
			hsize_t ct = chunkSlabsPerBlock(n_ids * n_comp, elem_bytes);
			// WP-163 R5: HDF5 refuses a chunk of >= 4 GiB, and the chunk spans every
			// id, so a slab above 2^29 f64 values made H5Dcreate fail (silently, the
			// result vanished). Above a 1 GiB slab, tile the id axis instead. Below
			// it the layout is unchanged (one chunk = whole slab x ct steps).
			hsize_t chunk_ids = n_ids;
			const hsize_t max_chunk_bytes = (hsize_t)1 << 30;  // 1 GiB
			const hsize_t row_bytes = n_comp * elem_bytes;
			if (n_ids * row_bytes > max_chunk_bytes) {
				chunk_ids = max_chunk_bytes / row_bytes;
				if (chunk_ids < 1) chunk_ids = 1;
			}
			hsize_t chunk[3] = { ct, chunk_ids, n_comp };
			herr_t st = H5Pset_chunk(dcpl, 3, chunk);
			if (st >= 0) st = H5Pset_shuffle(dcpl);
			if (st >= 0) st = H5Pset_deflate(dcpl, 4);
			hid_t dset = (st >= 0)
				? H5Dcreate(obj, name, disk_type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT)
				: HID_INVALID;
			H5Pclose(dcpl);
			H5Sclose(space);
			return (dset < 0) ? HID_INVALID : dset;
		}

		// ---- WP-164 / WP-233: write-path layout for datasets that stay OPEN -----
		//
		// WP-164 chunked DATA at ~1 MiB, stacking up to 1024 steps per chunk, and
		// deflated a partially filled chunk at every flush. Under SWMR (WP-233,
		// #1082) that layout leaks: the recompressed chunk grows, HDF5 moves it,
		// and SWMR never reuses the space it frees (a 300-step run at -flush 0
		// grew 2.6 -> 206 MB). Two layouts survive SWMR; both write a chunk's
		// final bytes in place or once:
		//
		//  * deflate > 0 (the default, -compress 1) and a slab of >= 4 KiB: ONE STEP
		//    PER CHUNK, ids tiled to <= 1 MiB, shuffle + deflate. A step's write
		//    completes every chunk it touches, so each chunk is compressed and
		//    written exactly once. About WP-164's size; reading one id's history
		//    touches T chunks. Below 4 KiB a per-step chunk compresses worse than
		//    raw (per-chunk overhead) and costs time: such a slab takes the
		//    unfiltered layout whatever -compress says (measured, WP-233 bench).
		//  * deflate == 0 (-compress 0): UNFILTERED, ~64 KiB chunks stacking up to
		//    1024 steps (ct <= 16 for a large slab). An unfiltered chunk keeps its
		//    size, so a partial chunk is rewritten in place (old rows byte-equal).
		//    It is allocated in full, so the 64 KiB target bounds the cost of the
		//    last, partial chunk of every result. Fastest writes and history reads.
		//
		// Measured: Ladruno_scripts/ladruno_recorder_tests/recorder_swmr_layout_bench.py.
		// The dataset access list sizes the chunk cache to one row of chunks (every
		// id tile x ct steps) so a chunk is completed in memory.
		struct TimeSeriesPlan {
			hsize_t ct;          // steps per chunk
			hsize_t chunk_ids;   // ids per chunk
			size_t cache_bytes;  // chunk-cache size for the dataset access list
			size_t cache_slots;  // hash slots for the chunk cache
		};

		// WP-233: does this slab get the compressed (one step per chunk) layout?
		inline bool timeSeriesCompressed(hsize_t n_ids, hsize_t n_comp, hsize_t elem_bytes,
		                                 int deflate) {
			return deflate > 0 && n_ids * n_comp * elem_bytes >= ((hsize_t)4 << 10);
		}

		inline TimeSeriesPlan planTimeSeries(hsize_t n_ids, hsize_t n_comp, hsize_t elem_bytes,
		                                     bool compressed) {
			const hsize_t row_bytes = n_comp * elem_bytes;
			const hsize_t slab = n_ids * row_bytes;
			TimeSeriesPlan p;
			if (compressed) {
				const hsize_t target = (hsize_t)1 << 20;      // <= 1 MiB per chunk
				p.ct = 1;
				p.chunk_ids = row_bytes ? target / row_bytes : 1;
				if (p.chunk_ids < 1) p.chunk_ids = 1;
				if (p.chunk_ids > n_ids) p.chunk_ids = n_ids;
			}
			else {
				const hsize_t target = (hsize_t)64 << 10;     // ~64 KiB per chunk
				const hsize_t row_budget = (hsize_t)2 << 20;  // ct x slab held in cache
				if (slab <= target) {
					p.chunk_ids = n_ids;
					p.ct = slab ? target / slab : 1;
					if (p.ct < 1) p.ct = 1;
					if (p.ct > 1024) p.ct = 1024;
				}
				else {
					p.ct = row_budget / slab;
					if (p.ct < 1) p.ct = 1;
					if (p.ct > 16) p.ct = 16;
					p.chunk_ids = target / (p.ct * row_bytes);
					if (p.chunk_ids < 1) p.chunk_ids = 1;
					if (p.chunk_ids > n_ids) p.chunk_ids = n_ids;
				}
			}
			const hsize_t tiles = (n_ids + p.chunk_ids - 1) / p.chunk_ids;
			const hsize_t chunk_bytes = p.ct * p.chunk_ids * row_bytes;
			// A row of chunks (every id tile x ct steps) plus one spare chunk. With
			// ct == 1 every chunk is completed by a single write, so a cache that
			// cannot hold the row costs nothing; cap it at 256 MiB.
			hsize_t cache = tiles * chunk_bytes + chunk_bytes;
			const hsize_t cap = (hsize_t)256 << 20;
			if (cache > cap) cache = cap;
			if (cache < ((hsize_t)1 << 20)) cache = (hsize_t)1 << 20;
			p.cache_bytes = (size_t)cache;
			p.cache_slots = (size_t)(tiles * 10 + 521) | 1;
			return p;
		}

		// The dataset access list (sized chunk cache, w0 = 1: evict fully written
		// chunks first) for a time-series DATA, at create and at every reopen.
		inline hid_t timeSeriesAccessList(hsize_t n_ids, hsize_t n_comp, hid_t disk_type, int deflate) {
			hsize_t elem_bytes = (H5Tequal(disk_type, H5T_IEEE_F32LE) > 0) ? 4 : 8;
			TimeSeriesPlan p = planTimeSeries(n_ids, n_comp, elem_bytes,
				timeSeriesCompressed(n_ids, n_comp, elem_bytes, deflate));
			hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
			if (dapl >= 0 && H5Pset_chunk_cache(dapl, p.cache_slots, p.cache_bytes, 1.0) < 0) {
				H5Pclose(dapl);
				return HID_INVALID;
			}
			return dapl;
		}

		// Create a [0 x nIds x nComp] extensible dataset with the plan above;
		// deflate 0 = unfiltered 64 KiB layout, deflate > 0 = one step per chunk.
		inline hid_t createTimeSeries3dOpen(hid_t obj, const char *name, hsize_t n_ids,
		                                    hsize_t n_comp, hid_t disk_type, int deflate) {
			if (n_ids < 1 || n_comp < 1) return HID_INVALID;
			hsize_t elem_bytes = (H5Tequal(disk_type, H5T_IEEE_F32LE) > 0) ? 4 : 8;
			const bool compressed = timeSeriesCompressed(n_ids, n_comp, elem_bytes, deflate);
			TimeSeriesPlan p = planTimeSeries(n_ids, n_comp, elem_bytes, compressed);
			hsize_t dims[3]    = { 0, n_ids, n_comp };
			hsize_t maxdims[3] = { H5S_UNLIMITED, n_ids, n_comp };
			hsize_t chunk[3]   = { p.ct, p.chunk_ids, n_comp };
			hid_t space = H5Screate_simple(3, dims, maxdims);
			if (space < 0) return HID_INVALID;
			hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
			hid_t dapl = timeSeriesAccessList(n_ids, n_comp, disk_type, deflate);
			herr_t st = (dcpl >= 0 && dapl >= 0) ? H5Pset_chunk(dcpl, 3, chunk) : -1;
			if (st >= 0 && compressed) {
				st = H5Pset_shuffle(dcpl);
				if (st >= 0) st = H5Pset_deflate(dcpl, (unsigned)(deflate > 9 ? 9 : deflate));
			}
			hid_t dset = (st >= 0)
				? H5Dcreate(obj, name, disk_type, space, H5P_DEFAULT, dcpl, dapl)
				: HID_INVALID;
			if (dapl >= 0) H5Pclose(dapl);
			if (dcpl >= 0) H5Pclose(dcpl);
			H5Sclose(space);
			return (dset < 0) ? HID_INVALID : dset;
		}

		// Write slab t of an open [T x nIds x nComp] dataset (extend to t+1 first).
		// `mspace` is a caller-owned [1 x nIds x nComp] memory space reused per step.
		inline herr_t writeSlab3dAt(hid_t dset, hsize_t t, const double *data,
		                            hsize_t n_ids, hsize_t n_comp, hid_t mspace) {
			hsize_t newdims[3] = { t + 1, n_ids, n_comp };
			if (H5Dset_extent(dset, newdims) < 0) return -1;
			hid_t fspace = H5Dget_space(dset);
			if (fspace < 0) return -1;
			hsize_t start[3] = { t, 0, 0 };
			hsize_t count[3] = { 1, n_ids, n_comp };
			herr_t status = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			if (status >= 0)
				status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, data);
			H5Sclose(fspace);
			return status;
		}

		// Write element t of an open [T] axis (extend to t+1 first).
		inline herr_t writeScalar1dAt(hid_t dset, hsize_t t, hid_t mem_type, const void *value,
		                              hid_t mspace1) {
			hsize_t newdims[1] = { t + 1 };
			if (H5Dset_extent(dset, newdims) < 0) return -1;
			hid_t fspace = H5Dget_space(dset);
			if (fspace < 0) return -1;
			hsize_t start[1] = { t }, count[1] = { 1 };
			herr_t status = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			if (status >= 0)
				status = H5Dwrite(dset, mem_type, mspace1, fspace, H5P_DEFAULT, value);
			H5Sclose(fspace);
			return status;
		}

		// WP-163 R5: current length of the unlimited (first) axis, 0 on error.
		inline hsize_t extent0(hid_t dset) {
			hid_t fspace = H5Dget_space(dset);
			if (fspace < 0) return 0;
			hsize_t cur[3] = { 0, 0, 0 };
			int nd = H5Sget_simple_extent_ndims(fspace);
			if (nd >= 1 && nd <= 3)
				H5Sget_simple_extent_dims(fspace, cur, NULL);
			H5Sclose(fspace);
			return cur[0];
		}

		// append one [1 x nIds x nComp] slab to a [T x nIds x nComp] dataset.
		inline herr_t appendSlab3d(hid_t dset, const double *data, hsize_t n_ids, hsize_t n_comp) {
			hid_t fspace = H5Dget_space(dset);
			hsize_t cur[3];
			H5Sget_simple_extent_dims(fspace, cur, NULL);
			H5Sclose(fspace);
			hsize_t t = cur[0];
			hsize_t newdims[3] = { t + 1, n_ids, n_comp };
			if (H5Dset_extent(dset, newdims) < 0) return -1;
			fspace = H5Dget_space(dset);
			hsize_t start[3] = { t, 0, 0 };
			hsize_t count[3] = { 1, n_ids, n_comp };
			H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			hsize_t mdims[3] = { 1, n_ids, n_comp };
			hid_t mspace = H5Screate_simple(3, mdims, NULL);
			herr_t status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, data);
			H5Sclose(mspace);
			H5Sclose(fspace);
			return status;
		}

		// create a [0] extensible 1-D axis of the given HDF5 file type.
		inline hid_t createTimeAxis1d(hid_t obj, const char *name, hid_t h5type) {
			hsize_t dims[1] = { 0 };
			hsize_t maxdims[1] = { H5S_UNLIMITED };
			hid_t space = H5Screate_simple(1, dims, maxdims);
			hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
			hsize_t chunk[1] = { 1024 };
			H5Pset_chunk(dcpl, 1, chunk);
			hid_t dset = H5Dcreate(obj, name, h5type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
			H5Pclose(dcpl);
			H5Sclose(space);
			return dset;
		}

		// append one scalar to a [T] extensible axis (mem type from the caller).
		inline herr_t appendScalar1d(hid_t dset, hid_t mem_type, const void *value) {
			hid_t fspace = H5Dget_space(dset);
			hsize_t cur[1];
			H5Sget_simple_extent_dims(fspace, cur, NULL);
			H5Sclose(fspace);
			hsize_t t = cur[0];
			hsize_t newdims[1] = { t + 1 };
			if (H5Dset_extent(dset, newdims) < 0) return -1;
			fspace = H5Dget_space(dset);
			hsize_t start[1] = { t }, count[1] = { 1 };
			H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			hsize_t mdims[1] = { 1 };
			hid_t mspace = H5Screate_simple(1, mdims, NULL);
			herr_t status = H5Dwrite(dset, mem_type, mspace, fspace, H5P_DEFAULT, value);
			H5Sclose(mspace);
			H5Sclose(fspace);
			return status;
		}
		inline herr_t appendDouble1d(hid_t dset, double value) {
			return appendScalar1d(dset, H5T_NATIVE_DOUBLE, &value);
		}
		inline herr_t appendInt1d(hid_t dset, int value) {
			return appendScalar1d(dset, H5T_NATIVE_INT, &value);
		}

		// ---- WP-232 (MPM particle group, 232a contract 2) ---------------------
		// The writers above cannot express three things the particle group
		// needs (quirk WP-232a-the-ladruno-hdf5-attribute-writer-skips-empty-
		// strings): an int64 dataset, a string-array dataset, and an EMPTY
		// dataset (createAndWrite skips an empty vector, returning HID_INVALID).
		// These write what they are given, empty included.

		// A [n] array of fixed-length NULLTERM strings, the width the longest
		// string plus its terminator (empty strings are written as "").
		inline hid_t createAndWriteStrings(hid_t obj, const char *name, const std::vector<std::string> &data)
		{
			size_t width = 1;
			for (size_t i = 0; i < data.size(); ++i)
				width = std::max(width, data[i].size() + 1);
			std::vector<char> buf(width * data.size(), '\0');
			for (size_t i = 0; i < data.size(); ++i)
				std::copy(data[i].begin(), data[i].end(), buf.begin() + i * width);
			hsize_t dim[1] = { (hsize_t)data.size() };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t atype = H5Tcopy(H5T_C_S1);
			H5Tset_size(atype, width);
			H5Tset_strpad(atype, H5T_STR_NULLTERM);
			hid_t dset = H5Dcreate(obj, name, atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if (dset >= 0 && !data.empty())
				H5Dwrite(dset, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
			H5Sclose(space);
			H5Tclose(atype);
			return dset;
		}

		// An empty dataset of any type and rank (dims may hold zeros), e.g. the
		// int32 [0 x 1] MODEL/NODES/ID of an empty stage.
		inline hid_t createEmpty(hid_t obj, const char *name, hid_t file_type, int rank, const hsize_t *dims)
		{
			hid_t space = H5Screate_simple(rank, dims, NULL);
			hid_t dset = H5Dcreate(obj, name, file_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			H5Sclose(space);
			return dset;
		}

		// A row-extensible dataset [0 x cols] (cols = 0: 1-D [0]), maxshape
		// [UNLIM (x cols)], chunked rows_per_chunk rows, shuffle + deflate at
		// `compress` (0 = no filter). Rows are appended by appendRows.
		inline hid_t createRows(hid_t obj, const char *name, hid_t file_type, hsize_t cols,
		                        hsize_t rows_per_chunk, int compress)
		{
			const int rank = cols > 0 ? 2 : 1;
			hsize_t dims[2] = { 0, cols };
			hsize_t maxd[2] = { H5S_UNLIMITED, cols };
			hsize_t chunk[2] = { rows_per_chunk > 0 ? rows_per_chunk : 1, cols };
			hid_t space = H5Screate_simple(rank, dims, maxd);
			hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
			H5Pset_chunk(dcpl, rank, chunk);
			if (compress > 0) {
				H5Pset_shuffle(dcpl);
				H5Pset_deflate(dcpl, (unsigned)compress);
			}
			hid_t dset = H5Dcreate(obj, name, file_type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
			H5Pclose(dcpl);
			H5Sclose(space);
			return dset;
		}

		// Append nrows rows (cols per row, 0 = 1-D) of mem_type data to a
		// dataset made by createRows. nrows = 0 is a no-op.
		inline herr_t appendRows(hid_t dset, hid_t mem_type, const void *data, hsize_t nrows, hsize_t cols)
		{
			if (nrows == 0) return 0;
			const int rank = cols > 0 ? 2 : 1;
			hid_t fspace = H5Dget_space(dset);
			hsize_t cur[2] = { 0, 0 };
			H5Sget_simple_extent_dims(fspace, cur, NULL);
			H5Sclose(fspace);
			hsize_t newd[2] = { cur[0] + nrows, cols };
			if (H5Dset_extent(dset, newd) < 0) return -1;
			fspace = H5Dget_space(dset);
			hsize_t start[2] = { cur[0], 0 };
			hsize_t count[2] = { nrows, cols };
			H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			hid_t mspace = H5Screate_simple(rank, count, NULL);
			herr_t status = H5Dwrite(dset, mem_type, mspace, fspace, H5P_DEFAULT, data);
			H5Sclose(mspace);
			H5Sclose(fspace);
			return status;
		}

	}

	namespace file {

		// low level functions for c interface

		inline herr_t close(hid_t file_id) {
			return H5Fclose(file_id);
		}
		inline herr_t flush(hid_t file_id) {
			return H5Fflush(file_id, H5F_SCOPE_LOCAL);
		}
		inline hid_t create(const char *filename, hid_t create_plist, hid_t acc_plist) {
			return H5Fcreate(filename, H5F_ACC_TRUNC, create_plist, acc_plist);
		}
		// WP-233: reopen read-write to lay out objects SWMR may not create.
		inline hid_t openRW(const char *filename, hid_t acc_plist) {
			return H5Fopen(filename, H5F_ACC_RDWR, acc_plist);
		}
		// WP-233: arm SWMR writing (latest format required; only appends after).
		inline herr_t startSWMR(hid_t file_id) {
			return H5Fstart_swmr_write(file_id);
		}
		// WP-233: the recorder's file access list: the latest format (SWMR needs
		// it; extensible-array chunk indexes). The close degree stays the default
		// (weak): HDF5's atexit close of a deck that never wipes must still flush
		// and close the file; openLayout() checks openObjects() before closing.
		inline hid_t accessList() {
			hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
			if (fapl < 0) return HID_INVALID;
			if (H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST) < 0) {
				H5Pclose(fapl);
				return HID_INVALID;
			}
			return fapl;
		}
		// WP-233: objects (datasets, groups, named types, attributes) still open
		// in the file. With a weak close degree, H5Fclose with any of them open
		// would leave the file open, still in SWMR mode, and the reopen would
		// then hand back that same open file.
		inline ssize_t openObjects(hid_t file_id) {
			return H5Fget_obj_count(file_id, H5F_OBJ_DATASET | H5F_OBJ_GROUP |
			                                 H5F_OBJ_DATATYPE | H5F_OBJ_ATTR | H5F_OBJ_LOCAL);
		}
	}

	namespace plist {

		enum CreateOptions {
			FileCreate,
			FileAccess,
			GroupCreate
		};

		// low level functions for c interface

		inline herr_t close(hid_t plist_id) {
			if (plist_id == HID_INVALID) return -1;
			return H5Pclose(plist_id);
		}
		inline hid_t crate(int opt) {
			switch (opt) {
			case FileCreate:
				return H5Pcreate(H5P_FILE_CREATE);
			case FileAccess:
				return H5Pcreate(H5P_FILE_ACCESS);
			case GroupCreate:
				return H5Pcreate(H5P_GROUP_CREATE);
			default:
				return HID_INVALID;
			}
		}
		inline herr_t setLinkCreationOrder(hid_t plist_id, unsigned int crt_order_flags) {
			return H5Pset_link_creation_order(plist_id, crt_order_flags);
		}
		inline herr_t setLibVerBounds(hid_t plist_id, unsigned int minor, unsigned int major) {
			return H5Pset_libver_bounds(plist_id, (H5F_libver_t)minor, (H5F_libver_t)major);
		}
		// WP-233 (#1082): a group that gains links during a run (the root gains a
		// MODEL_STAGE per stage; a modes result's DATA gains a STEP_<k> per
		// snapshot) keeps its links COMPACT for good and gets an object header
		// sized for `n_links` up front (at most ~1,500: one 64 KiB header chunk; past
		// that HDF5 ignores the estimate), so adding one link is ONE in-place write of
		// that header. By default a group goes dense after 8 links (a fractal heap
		// and v2 B-trees rewritten in place) or grows a continuation chunk, and a
		// kill inside the layout flush (which runs without SWMR) then left the
		// group unreadable: measured on the root, every stage lost
		// (recorder_crash_replay.py, --no-presize).
		inline herr_t presizeLinks(hid_t gcpl, unsigned n_links, unsigned name_len) {
			if (H5Pset_link_phase_change(gcpl, 65535, 65534) < 0) return -1;
			return H5Pset_est_link_info(gcpl, n_links, name_len);
		}
	}

}

} // namespace ladruno

#endif // Ladruno_Hdf5_h
