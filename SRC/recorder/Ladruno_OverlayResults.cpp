/* ********************************************************************** **
**  Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Ladruno_OverlayResults.cpp — OVERLAY pore-pressure ResultSource impl.    **
**                                                                        **
**  OverlayPressureSource reads the overlay's committed p through the      **
**  ADR-73 P4 const accessors (getRegionNodeTags / getCommittedP). It does **
**  NOT touch the HDF5 layer; persistence is the generic sink's job.       **
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

#include "Ladruno_OverlayResults.h"

#include "LadrunoPorousOverlay.h"

#include <Domain.h>
#include <LoadPattern.h>
#include <classTags.h>

#include <cstddef>
#include <sstream>

namespace ladruno {

	/* ===================================================================== */
	/* OverlayPressureSource                                                 */
	/* ===================================================================== */

	OverlayPressureSource::OverlayPressureSource(const detail::ProcessInfo& /*info*/,
	                                             LadrunoPorousOverlay* overlay)
		: m_patternTag(0)
		, m_ids()
		, m_schema()
	{
		if (overlay != 0) {
			// Snapshot the region-node id order (the canonical order for this
			// channel — matches the -record CSV header exactly).
			m_ids = overlay->getRegionNodeTags();
			m_patternTag = overlay->getTag();
		}
		buildSchema(m_patternTag);
	}

	void OverlayPressureSource::buildSchema(int patternTag)
	{
		std::stringstream ss;
		ss << "overlayPressure_" << patternTag;
		m_schema.name = ss.str();
		std::stringstream ssd;
		ssd << "Overlay pore pressure " << patternTag;
		m_schema.display_name = ssd.str();
		m_schema.components_csv = "p";
		m_schema.num_components = 1;
		m_schema.dimension = "F/L^2";
		m_schema.description = "LadrunoPorousOverlay committed pore-fluid pressure "
			"field (per region node)";
		m_schema.data_type = detail::ResultDataType::Scalar;
		m_schema.result_type = detail::ResultType::Generic;
	}

	void OverlayPressureSource::evaluate(const detail::ProcessInfo& info,
	                                     std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size(), 0.0);
		// Re-resolve BY TAG every step (F-1): `remove loadPattern` deletes the
		// pattern without a domainChange (no owned SPs), so a cached pointer
		// would dangle. Absent/replaced overlay => zero rows, never a deref.
		if (info.domain == 0)
			return;
		LoadPattern* lp = info.domain->getLoadPattern(m_patternTag);
		if (lp == 0 || lp->getClassTag() != PATTERN_TAG_LadrunoPorousOverlay)
			return;
		LadrunoPorousOverlay* ov = (LadrunoPorousOverlay*)lp;
		// Committed p aligned to getRegionNodeTags() order. Copy up to the id
		// count; a size mismatch (a re-added same-tag overlay with a different
		// region) leaves the remainder at 0.0 rather than reading out of bounds.
		const std::vector<double>& p = ov->getCommittedP();
		const size_t n = (p.size() < m_ids.size()) ? p.size() : m_ids.size();
		for (size_t i = 0; i < n; ++i)
			buffer[i] = p[i];
	}

} // namespace ladruno
