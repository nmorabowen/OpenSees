/* ********************************************************************** **
**  Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Ladruno_ResultIO.h — the Source / Sink contract.                        **
**                                                                        **
**  The core abstraction of the Ladruno recorder (ADR 03_ladruno_recorder): **
**  split *what to compute* (ResultSource) from *how to persist*          **
**  (ResultSink). Node, element, and domain results all implement         **
**  ResultSource; StreamingSink and EnvelopeSink are the two persistence  **
**  strategies. Envelopes and global scalars then compose for free across **
**  all three families instead of being bolted on per-family.            **
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

#ifndef Ladruno_ResultIO_h
#define Ladruno_ResultIO_h

#include "Ladruno_Types.h"

#include <string>
#include <vector>

namespace ladruno {

	// Self-describing column metadata for one result (schema v1 §7.1/§7.2).
	struct ResultSchema {
		std::string name;            // on-disk group name (e.g. "DISPLACEMENT", "stress")
		std::string display_name;    // human label
		std::string components_csv;  // "Ux,Uy,Uz"
		std::string dimension;       // "L", "F", "F*L", ...
		std::string description;
		int num_components = 0;
		detail::ResultType::Enum result_type = detail::ResultType::Generic;
		detail::ResultDataType::Enum data_type = detail::ResultDataType::Scalar;
	};

	/*
	ResultSource — WHAT to compute. A node / element / domain result implements
	this. The recorder pairs each source with a sink. A source yields, for the
	current step, a flat row-major buffer of shape [ ids().size() * schema.num_components ]
	(domain sources use a single synthetic id {0}). It is stateless across steps;
	cross-step behavior (streaming vs envelope) is the sink's job.
	*/
	class ResultSource {
	public:
		virtual ~ResultSource() = default;

		// Constant for the lifetime of the source.
		virtual const ResultSchema& schema() const = 0;

		// Row identity: node tags / element-gp ids / {0} for domain scalars.
		virtual const std::vector<int>& ids() const = 0;

		// Fill `buffer` (resized by the source) with this step's values,
		// row-major over ids() x num_components.
		virtual void evaluate(const detail::ProcessInfo& info,
		                      std::vector<double>& buffer) = 0;

		// Parallel routing (ADR D6): true for additive/global quantities that need a
		// per-step partition reduction before a sink may accumulate them (e.g.
		// boundary-node reactions, base shear). Consistent quantities (kinematics,
		// element results) return false and stream/envelope locally with zero comm.
		virtual bool requiresPartitionReduction() const { return false; }
	};

	/*
	ResultSink — HOW to persist. Wraps any source.
	  - StreamingSink writes RESULTS/.../DATA/STEP_<k> per step (parity behavior).
	  - EnvelopeSink keeps running min/max/absmax + arg-step and writes once at
	    finalize (or periodically), per MODEL_STAGE (schema §7.4).
	The sink owns the HDF5 group lifetime for its result.
	*/
	class ResultSink {
	public:
		virtual ~ResultSink() = default;

		// Create groups + the ID dataset for this source's result. Called once when
		// the recorder first records, and again on a new MODEL_STAGE.
		virtual void begin(detail::ProcessInfo& info, const ResultSource& src) = 0;

		// Persist one step's buffer (already partition-reduced if required).
		virtual void accept(detail::ProcessInfo& info, const ResultSource& src,
		                    const std::vector<double>& buffer) = 0;

		// Flush any deferred datasets (envelope writes its accumulators here).
		// Called on stage change, periodically, and at recorder teardown.
		virtual void finalize(detail::ProcessInfo& info) = 0;
	};

} // namespace ladruno

#endif // Ladruno_ResultIO_h
