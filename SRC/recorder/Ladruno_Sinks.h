/* ********************************************************************** **
**  Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Ladruno_Sinks.h — the two persistence sinks (Phase 2).                  **
**                                                                        **
**  StreamingSink — parity behavior: writes RESULTS/<family>/<name>/      **
**    DATA/STEP_<k> per step under the active MODEL_STAGE. Mirrors the    **
**    frozen ResultRecorder::record() group/ID/STEP structure exactly so  **
**    it feeds the 1e-12 parity gate.                                     **
**  EnvelopeSink — wraps any source; keeps componentwise MIN/MAX/ABSMAX + **
**    ARG_STEP accumulators (ADR D7) and writes ENVELOPES/<family>/<name> **
**    at finalize (and on periodic flush), per MODEL_STAGE (schema §7.4). **
**                                                                        **
**  Both write ONLY through ladruno::h5::* (group navigation — opening an    **
**  existing parent — uses the linked HDF5 C API directly; all *writes*   **
**  go through the wrapper).                                              **
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

#ifndef Ladruno_Sinks_h
#define Ladruno_Sinks_h

#include "Ladruno_ResultIO.h"
#include "Ladruno_Types.h"

#include <string>
#include <vector>

namespace ladruno {

	/*
	Which RESULTS sub-tree a sink writes into. Mirrors the schema families
	(§7.1 ON_NODES, §7.2 ON_ELEMENTS, §7.3 ON_DOMAIN) so one sink class serves
	all three source families. Also selects the ENVELOPES sub-tree (§7.4).
	*/
	struct ResultFamily {
		enum Enum {
			OnNodes = 0,
			OnElements,
			OnDomain,
			OnRegions
		};
		// schema group name: "ON_NODES" / "ON_ELEMENTS" / "ON_DOMAIN" /
		// "ON_REGIONS"
		static const char* groupName(Enum f);
	};

	/*
	StreamingSink — the parity sink.

	  begin():  create RESULTS/<family>/<name> result group (DISPLAY_NAME /
	            COMPONENTS / DIMENSION / DESCRIPTION / TYPE / DATA_TYPE attrs),
	            the ID dataset [nIds x 1] from src.ids(), and the empty DATA group.
	  accept(): create <name>/DATA/STEP_<commitTag> [nIds x nComp] from the
	            buffer, with STEP + TIME attrs.
	  finalize(): nothing deferred (handles are opened/closed per call).

	The on-disk shape/attr names match the frozen recorder verbatim; only the
	parent location differs (under MODEL_STAGE[..]/RESULTS/<family> instead of
	the file root).
	*/
	class StreamingSink : public ResultSink {
	public:
		explicit StreamingSink(ResultFamily::Enum family)
			: m_family(family), m_initialized(false) {}
		~StreamingSink() override = default;

		void begin(detail::ProcessInfo& info, const ResultSource& src) override;
		void accept(detail::ProcessInfo& info, const ResultSource& src,
		            const std::vector<double>& buffer) override;
		void finalize(detail::ProcessInfo& info) override;

		// New MODEL_STAGE / re-begin: forget that the group was created.
		void reset() { m_initialized = false; }

	private:
		ResultFamily::Enum m_family;
		bool m_initialized; // group + ID written for the current stage
	};

	/*
	EnvelopeSink — componentwise running extremes (ADR D7), per MODEL_STAGE.

	Accumulators are all shape [nIds x nComp], row-major:
	  m_min     : running minimum per component
	  m_max     : running maximum per component
	  m_absmax  : running max(|value|) per component
	  m_arg_step: commitTag at which m_absmax for that component was last set
	Init happens on the FIRST accept() (we don't know nIds/nComp until then; a
	source's ids()/schema() are stable for the source lifetime but we seed from
	the first buffer so MIN/MAX start at the real first sample, not +/-inf).

	  accept():   element-wise update of all four accumulators from the step
	              buffer; first call seeds them.
	  flush():    rewrite the ENVELOPES datasets in place (delete-and-recreate)
	              for crash safety; callable any time after the first accept.
	  finalize(): same as flush() — the final authoritative write for the stage.
	  reset():    drop accumulators (cheap) so the sink can be reused on the next
	              MODEL_STAGE.
	*/
	class EnvelopeSink : public ResultSink {
	public:
		explicit EnvelopeSink(ResultFamily::Enum family)
			: m_family(family), m_seeded(false), m_n_ids(0), m_n_comp(0),
			  m_result_type(0), m_data_type(0) {}
		~EnvelopeSink() override = default;

		void begin(detail::ProcessInfo& info, const ResultSource& src) override;
		void accept(detail::ProcessInfo& info, const ResultSource& src,
		            const std::vector<double>& buffer) override;
		void finalize(detail::ProcessInfo& info) override;

		// Periodic crash-safety rewrite of the envelope datasets (schema §7.4 /
		// §7b). Identical payload to finalize(); separated so the recorder can
		// call it on its flush cadence without semantic meaning of "stage done".
		void flush(detail::ProcessInfo& info);

		// Drop all accumulators for a new MODEL_STAGE. Cheap (vector clears).
		void reset();

	private:
		// shared writer for flush()/finalize()
		void writeEnvelope(detail::ProcessInfo& info);

	private:
		ResultFamily::Enum m_family;
		bool m_seeded;            // accumulators initialized from first accept
		size_t m_n_ids;
		size_t m_n_comp;
		std::string m_name;       // result name (schema().name); cached at begin/first accept
		// Self-describing result-group metadata (schema §7.1), cached so the
		// ENVELOPES/<name> group carries the same COMPONENTS/DISPLAY_NAME/DIMENSION
		// as the time-series result group (finding (h): authoritative names in-file).
		std::string m_display_name, m_components_csv, m_dimension, m_description;
		int m_result_type, m_data_type;
		std::vector<int> m_ids;   // cached ids() snapshot for the ID dataset
		std::vector<double> m_min;
		std::vector<double> m_max;
		std::vector<double> m_absmax;
		std::vector<int> m_arg_step; // per-component step index of the abs-extreme
	};

} // namespace ladruno

#endif // Ladruno_Sinks_h
