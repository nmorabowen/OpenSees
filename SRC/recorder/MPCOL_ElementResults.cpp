/* ********************************************************************** **
**  MPCO_Ladruno recorder — MPCOL_ElementResults.cpp                      **
**  ElementResultSource adapter implementation (Phase 3, Step 2/3).       **
**                                                                        **
**  Wraps one already-built (header, response-collection) bucket — the    **
**  output of the orchestrator's initElementRecorders port — as a         **
**  mpcol::ResultSource. evaluate() reproduces the frozen                  **
**  recordResultsOnElements per-step packing [nElems x NUM_COLUMNS]       **
**  byte-for-byte (1e-12 parity gate). buildSchema() turns the            **
**  OutputDescriptorHeader into the ResultSchema + carries the            **
**  COLUMN_MAP-equivalent metadata (schema §7.2) for the recorder to      **
**  persist after the sink creates the result group.                       **
** ********************************************************************** */

#include "MPCOL_ElementResults.h"

#include <sstream>

namespace mpcol {

	ElementResultSource::ElementResultSource(const std::vector<std::string>& request,
		const mpco::element::OutputDescriptorHeader& header,
		mpco::element::OutputResponseCollection& bucket)
		: m_header(header)
		, m_bucket(bucket)
	{
		// Row identity = the element tag of each response in the bucket (matches
		// the frozen recordResultsOnElements ID dataset order).
		m_ids.resize(m_bucket.items.size());
		for (size_t i = 0; i < m_bucket.items.size(); ++i)
			m_ids[i] = m_bucket.items[i].element->getTag();
		buildSchema(request);
	}

	void ElementResultSource::buildSchema(const std::vector<std::string>& request)
	{
		// The on-disk group name is the per-bucket dir_name produced by the
		// orchestrator (e.g. "204-FourNodeQuad[201:0:0]"), nested under the
		// per-result display name (e.g. "force" / "stress"). The recorder pairs
		// this source with a StreamingSink(OnElements) whose family root is
		// ON_ELEMENTS; the source name therefore embeds "<display>/<bucket>".
		std::stringstream ss_disp;
		for (size_t i = 0; i < request.size(); ++i) {
			if (i > 0) ss_disp << ".";
			ss_disp << request[i];
		}
		m_schema.display_name = ss_disp.str();
		// name = "<result_display>/<bucket_dir>" so the sink lands it at
		// ON_ELEMENTS/<result_display>/<bucket_dir>.
		std::stringstream ss_name;
		ss_name << m_schema.display_name << "/" << m_bucket.dir_name;
		m_schema.name = ss_name.str();

		// Element results: components/dimension are not a single CSV — the column
		// structure lives in the COLUMN_MAP (written separately by the recorder).
		// num_components is the flat NUM_COLUMNS so the sink writes [nElem x NCOL].
		m_schema.components_csv = "";
		m_schema.dimension = "";
		m_schema.description = "";
		m_schema.num_components = m_header.num_columns;
		m_schema.result_type = mpco::ResultType::Generic;
		m_schema.data_type = mpco::ResultDataType::Scalar;
	}

	void ElementResultSource::evaluate(const mpco::ProcessInfo& /*info*/, std::vector<double>& buffer)
	{
		// Byte-faithful port of the frozen recordResultsOnElements inner packing:
		//   for each element response in the bucket: getResponse(); copy its
		//   getInformation().getData() row into [i * NUM_COLUMNS .. +NUM_COLUMNS).
		const int num_columns = m_header.num_columns;
		const size_t num_rows = m_bucket.items.size();
		buffer.assign(num_rows * (size_t)num_columns, 0.0);
		for (size_t i = 0; i < num_rows; ++i) {
			mpco::element::OutputResponse& current_response = m_bucket.items[i];
			current_response.response->getResponse();
			const Vector& current_data = current_response.response->getInformation().getData();
			if (current_data.Size() != num_columns) {
				opserr << "MPCORecorderLadruno Error: invalid response data size\n";
				continue;
			}
			size_t offset = i * (size_t)num_columns;
			for (int j = 0; j < num_columns; ++j)
				buffer[offset + (size_t)j] = current_data[j];
		}
	}

} // namespace mpcol
