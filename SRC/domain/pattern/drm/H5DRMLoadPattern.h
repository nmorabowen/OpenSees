/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */


// ============================================================================
// 2019 By Jose Abell @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// 
// H5DRM provides a load pattern that implements reading H5DRM datasets
// (based on HDF5) containing ground motions useful in the context of
// soil-structure interaction modeling using the Domain Reduction method.
//
// Important assumptions and information:
//
//  + A brute force method for matching FEM domiain nodes to those of the 
//    DRM dataset is used. There is space for improvements here (A scale 
//    factor can be applied to the DRM dataset point coordinates to deal 
//    with non-matching units, DRM datasets should be consistently 
//    created using base SI units and the coordinate system convention).
//    Nodal matching is sufficiently fast for problems with less than
//    1M nodes, it only occurs once during analysis.
//  + Coordinate system for H5DRM datasets coming from ShakerMaker: Z is up, 
//    X is east and Y is
//    north. Units should be SI base units always. The coordinate system
//    of the FEM model should match that of the DRM except, maybe, in
//    in the units of the point coordinates. If the systems dont match, there
//    are option to transform these dataset coordinates into a different 
//    coordinate. 
//  + This H5DRM pattern currently only works for 3-DOF elements of maxumim 8 nodes.  
//  + H5DRM dataset should have at least the acceleration and displacement
//    datasets available for this to work. 
//  + Linear interpolation is done on displacements and accelerations if
//    current analysis timestep does not "fall" in any of the motion time 
//    samples. This can be improved since a linear interpolation in 
//    accelerations is a cubic interpolation in displacements.  
//  + Based on the original papers by Jacobo Bielak [1] and [2]. 
//  + The previous remark means that velocities are not directly used. 
//    This might be important if damping in the DRM layer is thought to be 
//    needed. Please let me (jaabell) know if you think this is needed in
//    some case, I would be happy to discuss and implement this and also
//    know why it might be needed. (coupled poroelasticity?)
//  + It has been tested to work in parallel. 
//  + The H5DRM data format specification is documented in [3].
//  + Based on and and extends the work done for my PhD thesis [4].
// 
//  References: 
//   [1] Bielak, J., Loukakis, K., Hisada, Y., & Yoshimura, C. (2003). 
//       Domain Reduction Method for Three--Dimensional Earthquake 
//       Modeling in Localized Regions. Part {I}: Theory. Bulletin of the 
//       Seismological Society of America, 93(2), 817–824.
//   [2] Yoshimura, C., Bielak, J., & Hisada, Y. (2003). Domain Reduction 
//       Method for Three--Dimensional Earthquake Modeling in Localized 
//       Regions. Part {II}: Verif\/ication and Examples. Bulletin of the 
//       Seismological Society of America, 93(2), 825–840.
//   [3] José A. Abell, Jorge G.F. Crempien, Matías Recabarren, ShakerMaker: A 
//       framework that simplifies the simulation of seismic ground-motions, 
//       SoftwareX, Volume 17, 2022
//   [4] Abell, J. A. (2016). Earthquake-Soil-Structure Interaction Modeling 
//       of Nuclear Power Plants for Near-Field Events. 
//       PhD Thesis. University of California, Davis.
//
// Citation: please cite [3] if you find H5DRM useful. 
// 
// ============================================================================

#ifndef H5DRM_h
#define H5DRM_h

#ifdef _H5DRM

#include <map>
#include <vector>
#include <utility>    // Ladruno (ADR-88): std::pair for interpNodeParents
#include <algorithm>
#include <string>

#include <hdf5.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <ErrorHandler.h>

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>

#include <LoadPattern.h>

#define H5DRM_PREALLOC_TSTEPS 10
#define H5DRM_MAX_RETURN_OPEN_OBJS 100
#define H5DRM_MAX_FILENAME 200
#define H5DRM_MAX_STRINGSIZE 80
#define H5DRMerror std::cerr  << "H5DRM(" << MPI_local_rank << ") - ERROR: "
#define H5DRMwarning std::cout  <<  "H5DRM(" << MPI_local_rank << ") - Warning: "
#define H5DRMout std::cout  << "H5DRM(" << MPI_local_rank << ") : "
// Ladruno WP-187: per-step diagnostics, printed only with the pattern's -verbose flag
// (the if/else form keeps a following `else` from binding to the macro's `if`).
#define H5DRMstep if (!opt_verbose_steps) {} else H5DRMout


class H5DRMLoadPattern : public LoadPattern
{
public:
    H5DRMLoadPattern();
    H5DRMLoadPattern(int tag, std::string dataset_fname_, double cFactor_ = 1.0, 
          double crd_scale_ = 1, 
          double distance_tolerance_ = 1e-3, 
          bool do_coordinate_transformation = true,
        double T00 = 1.0, double T01 = 0.0, double T02 = 0.0,
        double T10 = 0.0, double T11 = 1.0, double T12 = 0.0,
        double T20 = 0.0, double T21 = 0.0, double T22 = 1.0,
        double x00 = 0.0, double x01 = 0.0, double x02 = 0.0);
    ~H5DRMLoadPattern();
    void cleanup(); 

    void setDomain(Domain *theDomain);
    void applyLoad(double time);
    void Print(ostream &s, int flag = 0);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                    FEM_ObjectBroker &theBroker);
    LoadPattern *getCopy(void);

    // Ladruno WP-187: read options (trailing pattern flags, see OPS_H5DRMLoadPattern).
    //   verbose        : print the per-step "H5DRM(r) : t = ..." diagnostics (default off)
    //   prefetchSteps  : time steps per in-memory block (-1 = the dataset's chunk length,
    //                    0 = off: one 3x1 hyperslab read per node and step, as before)
    //   chunkCacheMB   : HDF5 chunk cache per dataset (-1 = sized from the chunk dims and
    //                    the rows this rank reads, 0 = HDF5's default 1 MiB)
    void setReadOptions(bool verbose, int prefetchSteps, double chunkCacheMB);

    virtual std::string getName() const
    {
        return "H5DRM: " + dataset_fname;
    }

protected:

    bool  CalculateBoundaryForces(double t);
    bool  ComputeDRMMotions(double next_integration_time);
    bool  drm_differentiate_displacements(double next_integration_time);
    bool  drm_integrate_velocity(double next_integration_time);
    bool  drm_direct_read(double next_integration_time);
    Vector *getNodalLoad(int node, double time);

    void do_intitialization();

    void node_matching_BruteForce(double d_tol, const ID& internal, const Matrix& xyz, const Vector& drmbox_x0, double& d_err, int & n_nodes_found);

    // Ladruno (ADR-88): reconstructs DRM_D/DRM_A for every node resolved via
    // element-provided interpolation (see interpNodeParents / do_intitialization)
    // from its parents' already-computed values. Called once per step, right
    // after drm_direct_read/drm_differentiate_displacements.
    void resolveInterpolatedNodeMotions();

    // Ladruno WP-187: chunk-cache-aware dataset open and the per-rank block prefetch.
    hid_t open_drm_dataset(const char *path);
    void  collect_rank_rows();
    bool  ensure_block(int i_lo, int i_hi);
    bool  block_value(hid_t dset, int data_pos, int i, double out[3]);

private:

    std::string dataset_fname;   // Name of the HDF5 dataset containing the DRM motions  
                
    ID DRM_Elements;                
    ID DRM_Nodes;                  
    ID DRM_Boundary_Flag;              

    Vector DRM_F;    // vector containing the DRM equivalent-loads
    Vector DRM_D;    // vector containing the displacements at DRM boundary nodes
    Vector DRM_D0;   // vector containing the initial displacements at DRM boundary nodes
    Vector DRM_A;    // vector containing the accelerations at DRM boundary nodes
    
    double last_integration_time;
    double t1, t2, tstart, tend, dt;    // specifies the time increment used in load path vector
    double cFactor;                     // additional factor on the returned load factor
    double crd_scale;                   // Scaling for the point coordinates of the DRM dataset
    double distance_tolerance;          // Distance tolerance for node-matching algorithm

    double drmbox_xmax, drmbox_xmin, drmbox_ymax, drmbox_ymin, drmbox_zmax, drmbox_zmin;


    int N_timesteps;
    bool flag_initialized;
    int N_local_elements;

    std::map<int, int> nodetag2station_id;
    std::map<int, int> nodetag2local_pos;
    ID station_id2data_pos;

    // Ladruno (ADR-88): node tag -> list of (parent node tag, weight) pairs,
    // for nodes resolved via DRMHigherOrderNode interpolation instead of a
    // direct H5DRM station match (e.g. BezierTet10's 6 mid-edge nodes). Absent
    // entry = node matched a real station directly (pre-ADR-88 behavior,
    // unchanged). Populated in do_intitialization, consumed in
    // resolveInterpolatedNodeMotions.
    std::map<int, std::vector<std::pair<int, double> > > interpNodeParents;

    //HDF5 handles
    hid_t ih5_fname;
    hid_t ih5_vel;
    hid_t ih5_vel_ds;
    hid_t ih5_dis;
    hid_t ih5_dis_ds;
    hid_t ih5_acc;
    hid_t ih5_acc_ds;
    hid_t ih5_one_node_ms;
    hid_t ih5_xfer_plist;

    // Ladruno WP-187: read options and the per-rank time-block cache. The block holds,
    // for the rows this rank reads (3 per station, sorted, unique), time steps
    // [blk_t0, blk_t0 + blk_len) of displacement and acceleration, row-major
    // (row r at blk_dis[r * blk_len]). Values are copied bytes of the file: the
    // interpolation that follows is unchanged, so results are bit-identical.
    bool   opt_verbose_steps = false;
    int    opt_prefetch_steps = -1;
    double opt_chunk_cache_mb = -1.0;
    bool   rows_collected = false;
    std::vector<hsize_t> rank_rows;          // sorted unique dataset rows this rank reads
    std::map<int, int>   datapos2blkrow;     // data_pos -> block row of its first component
    std::vector<std::pair<hsize_t, hsize_t> > row_runs;   // (first row, n rows) contiguous runs
    std::vector<double> blk_dis;
    std::vector<double> blk_acc;
    int    blk_t0 = -1;
    int    blk_len = 0;
    int    blk_step = 0;                     // steps per refill (resolved prefetch length)
    hsize_t ds_ntime = 0;                    // time extent of the displacement dataset
    hsize_t ds_chunk[2] = {0, 0};            // chunk dims of the displacement dataset (0 = not chunked)
    long long n_block_refills = 0;

    int MPI_local_rank;        // MPI Process-id (rank) in the case of parallel processing

    bool do_coordinate_transformation;

    //Coordinate transformation... xyz_new = T xyz_old + x0
    Matrix T;
    Vector x0;
};

#endif // _H5DRMLoadPattern

#endif
