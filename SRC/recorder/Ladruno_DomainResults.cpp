/* ********************************************************************** **
**  Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Ladruno_DomainResults.cpp — DOMAIN / REGION ResultSource impls (ADR D8). **
**                                                                        **
**  EnergyBalanceSource computes the energy balance through the shared     **
**  ebkernel math (EnergyBalanceKernel.h), so its numerical output matches **
**  EnergyBalanceRecorder's whole-model path exactly. It does NOT touch    **
**  the HDF5 layer; persistence is the sink's job (StreamingSink).         **
** ********************************************************************** */

#include "Ladruno_DomainResults.h"

#include "EnergyBalanceKernel.h"

#include "Domain.h"
#include "Element.h"
#include "ElementIter.h"
#include "MeshRegion.h"

#include <cstddef>

namespace ladruno {

	// Size a scratch vector to the largest element DOF count in the domain so
	// the per-entity sweep performs no per-call allocation (mirrors the
	// recorder's buildCache() velScratch sizing).
	static void sizeVelScratch(Domain* dom, Vector& velScratch)
	{
		if (dom == 0)
			return;
		int maxNumDOF = 0;
		Element* ele;
		ElementIter& elements = dom->getElements();
		while ((ele = elements()) != 0) {
			const int n = ele->getNumDOF();
			if (n > maxNumDOF)
				maxNumDOF = n;
		}
		if (velScratch.Size() < maxNumDOF)
			velScratch.resize(maxNumDOF);
	}

	/* ===================================================================== */
	/* EnergyBalanceSource                                                   */
	/* ===================================================================== */

	EnergyBalanceSource::EnergyBalanceSource(const detail::ProcessInfo& /*info*/)
		: m_per_region(false)
		, m_ids()
		, m_schema()
		, m_model_acc()
		, m_region_acc()
		, m_prev_time(0.0)
		, m_first(true)
		, m_velScratch()
	{
		m_ids.push_back(0); // synthetic domain id {0}
		buildSchema();
	}

	EnergyBalanceSource::EnergyBalanceSource(const detail::ProcessInfo& /*info*/,
	                                         const std::vector<int>& region_tags)
		: m_per_region(true)
		, m_ids(region_tags)
		, m_schema()
		, m_model_acc()
		, m_region_acc(region_tags.size())
		, m_prev_time(0.0)
		, m_first(true)
		, m_velScratch()
	{
		buildSchema();
	}

	void EnergyBalanceSource::buildSchema()
	{
		m_schema.name = "energyBalance";
		m_schema.display_name = "Energy Balance";
		m_schema.components_csv = "KE,IE,DW,ULW,RES,ERR";
		m_schema.num_components = ebkernel::NUM_ENERGY_COMPONENTS; // 6
		m_schema.dimension = "F*L";
		m_schema.description = "Structural-dynamics energy balance "
			"(kinetic/internal/damping/load-work + closure residual)";
		m_schema.data_type = detail::ResultDataType::Vectorial;
		m_schema.result_type = detail::ResultType::Generic;
	}

	void EnergyBalanceSource::evaluate(const detail::ProcessInfo& info,
	                                   std::vector<double>& buffer)
	{
		const size_t ncomp = (size_t)m_schema.num_components;
		buffer.assign(m_ids.size() * ncomp, 0.0);

		Domain* dom = info.domain;
		if (dom == 0)
			return;

		// keep the element-DOF scratch sized to the current domain
		sizeVelScratch(dom, m_velScratch);

		// dT from the previous recorded time (guard the first step like the
		// recorder: the first evaluate() seeds the rates and takes no spurious
		// increment over the initial time jump).
		const double dT = info.current_time_step - m_prev_time;

		if (!m_per_region) {
			double ke = 0., ir = 0., dr = 0., ur = 0.;
			ebkernel::sweepDomain(dom, m_velScratch, ke, ir, dr, ur);
			double out[ebkernel::NUM_ENERGY_COMPONENTS];
			m_model_acc.step(dT, m_first, ke, ir, dr, ur, out);
			for (size_t c = 0; c < ncomp; ++c)
				buffer[c] = out[c];
		}
		else {
			for (size_t i = 0; i < m_ids.size(); ++i) {
				MeshRegion* reg = dom->getRegion(m_ids[i]);
				double ke = 0., ir = 0., dr = 0., ur = 0.;
				ebkernel::sweepRegion(dom, reg, m_velScratch, ke, ir, dr, ur);
				double out[ebkernel::NUM_ENERGY_COMPONENTS];
				m_region_acc[i].step(dT, m_first, ke, ir, dr, ur, out);
				const size_t base = i * ncomp;
				for (size_t c = 0; c < ncomp; ++c)
					buffer[base + c] = out[c];
			}
		}

		m_first = false;
		m_prev_time = info.current_time_step;
	}

} // namespace ladruno
