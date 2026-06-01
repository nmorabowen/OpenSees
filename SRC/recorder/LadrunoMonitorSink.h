/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno: live analysis-monitor SWMR-HDF5 sink.
//
// Written for the Ladruno fork (N. Mora-Bowen). See
// Ladruno_implementation/08_analysis_monitor.md.
//
// A LadrunoMonitorSink owns one HDF5 file written in SWMR
// (Single-Writer / Multiple-Reader) mode so a separate viewer process can
// tail it *while the analysis is still running*. The on-disk layout is:
//
//   /              attrs: FORMAT, FORMAT_VERSION, GENERATOR, units(optional)
//     COLUMNS      [nCols]   variable-length strings (self-describing channels)
//     STEP         [T]       int32, chunked/unlimited  (commitTag per frame)
//     TIME         [T]       float64, chunked/unlimited (pseudo-time per frame)
//     FRAMES       [T,nCols] float64, chunked/unlimited (one row appended/frame)
//
// SWMR rules honoured here:
//   * file created with H5F_LIBVER_LATEST bounds,
//   * every dataset is chunked (required for extensible + SWMR),
//   * the FULL object layout (incl. COLUMNS) is created BEFORE
//     H5Fstart_swmr_write -- no new objects may appear afterwards, which is
//     why the channel set must be fixed at recorder construction,
//   * each append() H5Dflush()es so readers see the new frame.
//
// The same file is a valid at-rest results file once the run ends, so the
// existing profiler-viewer stack can open it unchanged (one file, two
// lifetimes -- see 08_analysis_monitor.md, resolved decision #1).

#ifndef LadrunoMonitorSink_h
#define LadrunoMonitorSink_h

#include <string>
#include <vector>

#include "hdf5.h"

namespace ladruno {

class MonitorSink {
public:
    // columns: one self-describing label per scalar channel, e.g.
    //   "node5.disp.dof1". units: optional free-form string ("" = omit).
    MonitorSink(const std::string &filename,
                const std::vector<std::string> &columns,
                const std::string &units = "");
    ~MonitorSink();

    // Create the file + datasets and arm SWMR. Returns 0 on success, <0 on
    // failure (after which isOpen() is false and append() is a no-op).
    int open();

    // Append one frame. values must have columns().size() entries. Never
    // blocks on a reader. Returns 0 on success, <0 on a write error.
    int append(int step, double time, const std::vector<double> &values);

    // Flush + close. Safe to call more than once. Called from the dtor; we do
    // NOT close in any static/atexit path (heeds the MPCO exit(-1) incident /
    // the profiler "flushing is the command layer's job" rule).
    void close();

    bool isOpen() const { return fileId >= 0; }
    size_t numColumns() const { return columns.size(); }
    long long numFramesWritten() const { return nFrames; }

private:
    MonitorSink(const MonitorSink &) = delete;
    MonitorSink &operator=(const MonitorSink &) = delete;

    int extendAndWrite(hid_t dset, hid_t memType, hsize_t ncols,
                       const void *rowData);

    std::string filename;
    std::vector<std::string> columns;
    std::string units;

    hid_t fileId;     // H5I_INVALID_HID until open()
    hid_t dsStep;
    hid_t dsTime;
    hid_t dsFrames;
    long long nFrames;
    bool writeFailed; // latched on first write error -> further appends no-op
};

} // namespace ladruno

#endif
