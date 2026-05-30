#ifndef MPCOL_Hdf5_h
#define MPCOL_Hdf5_h

/*************************************************************************************

MPCOL_Hdf5.h

HDF5 interface layer (namespace h5) ported VERBATIM from MPCORecorder.cpp
(frozen). Every function is made `inline` because, unlike the frozen .cpp where
these are plain TU-local functions, here they live in a header included by
multiple TUs and must avoid duplicate-symbol errors. Everything is wrapped in
`namespace mpcol { ... }` to avoid ODR clashes with the frozen file's identical
`namespace h5` symbols.

Context: in this build the macro MPCO_HDF5_LOADED_AT_RUNTIME is NEVER defined,
so the real linked HDF5 C API is used directly (H5*). No runtime loader /
ptr_H5* machinery is involved.

**************************************************************************************/

#include "MPCOL_Types.h"

namespace mpcol {

/*mpco - hdf5 interface utilities*/
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
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<mpco::element::FiberData> &data)
		{
			if (data.size() > 0) {
				return createAndWrited2(obj, name, data[0].data(), data.size(), 3);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<mpco::element::ElemGaussPair> &data)
		{
			if (data.size() > 0) {
				return createAndWritei2(obj, name, data[0].data(), data.size(), 2);
			}
			return HID_INVALID;
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
#ifdef MPCO_USE_SWMR
		inline herr_t startSWMR(hid_t file_id) {
			return H5Fstart_swmr_write(file_id);
		}
#endif // MPCO_USE_SWMR
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
	}

}

} // namespace mpcol

#endif // MPCOL_Hdf5_h
