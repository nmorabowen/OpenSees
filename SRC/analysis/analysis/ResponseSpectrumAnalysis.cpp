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

// $Revision: 1.0 $
// $Date: 2020-11-27 18:07:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/ResponseSpectrumAnalysis.h,v $

// Written: Massimo Petracca (ASDEA Software) 
// Created: Fri Nov 27 18:07:11: 2020
// Revision: A
//
// Description: This file contains the class definition for ResponseSpectrumAnalysis.
//
// What: "@(#) ResponseSpectrumAnalysis.h, revA"

#include <ResponseSpectrumAnalysis.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <DomainModalProperties.h>
#include <TimeSeries.h>
#include <elementAPI.h>
#include <Node.h>
#include <NodeIter.h>
#include <LadrunoModalCombination.h> // Ladruno ADR44 P1b
#include <LadrunoRitzKernel.h>       // Ladruno WP-213: Sa-period override
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <cstdio>  // Ladruno WP-212
#include <cstdlib> // Ladruno WP-212

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

namespace {

	bool string_to_list_of_doubles(const std::string& text, char sep, std::vector<double>& out) {
		auto to_double = [](const std::string& text, double& num) -> bool {
			num = 0.0;
			try {
				num = std::stod(text);
				return true;
			}
			catch (...) {
				return false;
			}
		};
		if (out.size() > 0) out.clear();
		std::size_t start = 0, end = 0;
		double value;
		while (true) {
			end = text.find(sep, start);
			if (end == std::string::npos) {
				if (start < text.size()) {
					if (!to_double(text.substr(start), value))
						return false;
					out.push_back(value);
				}
				break;
			}
			std::string subs = text.substr(start, end - start);
			if (subs.size() > 0) {
				if (!to_double(subs, value))
					return false;
				out.push_back(value);
			}
			start = end + 1;
		}
		return true;
	}

}

int
OPS_ResponseSpectrumAnalysis(void)
{
	// responseSpectrum $tsTag $dir <-scale $scale>

	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ResponseSpectrumAnalysis - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// get analysis model
	AnalysisModel* theAnalysisModel = *OPS_GetAnalysisModel();
	if (theAnalysisModel == nullptr) {
		opserr << "ResponseSpectrumAnalysis Error: no AnalysisModel available.\n";
		return -1;
	}
	if (theAnalysisModel->getDomainPtr() == nullptr) {
		opserr << "ResponseSpectrumAnalysis Error: no Domain available.\n";
		return -1;
	}

	// Ladruno WP-207: `responseSpectrumAnalysis -revert` on its own drops the
	// trial state an -increment call left behind (back to the committed state).
	if (OPS_GetNumRemainingInputArgs() == 1) {
		const char* first = OPS_GetString();
		if (strcmp(first, "-revert") == 0) {
			if (theAnalysisModel->revertDomainToLastCommit() < 0) {
				opserr << "ResponseSpectrumAnalysis Error: -revert failed to revert the domain.\n";
				return -1;
			}
			return 0;
		}
		OPS_ResetCurrentInputArg(-1);
	}

	// init default arguments
	TimeSeries* ts = nullptr;
	int dir = 1;
	double scale = 1.0;
	std::vector<double> Tn;
	std::vector<double> Sa;
	int mode_id = 0;
	bool single_mode = false;
	// Ladruno ADR44 P1b: opt-in modal combination
	int combine_rule = LadrunoModalCombination::RULE_NONE;
	std::vector<double> xi_list;
	// Ladruno WP-207: opt-in incremental mode about the committed state
	bool increment = false;
	bool sign_given = false;
	bool revert_at_end = false;
	double sign = 1.0;
	// Ladruno WP-212: required cumulative mass ratio (0 => warn only)
	double require_mass = 0.0;

	// make sure eigenvalue and modal properties have been called before
	DomainModalProperties modal_props;
	if (theAnalysisModel->getDomainPtr()->getModalProperties(modal_props) < 0) {
		opserr << "ResponseSpectrumAnalysis - eigen and modalProperties have not been called" << endln;
		return -1;
	}
	int ndf = modal_props.totalMass().Size();

	// parse
	int nargs = OPS_GetNumRemainingInputArgs();
	if (nargs < 2) {
		opserr << "ResponseSpectrumAnalysis $tsTag $dir <-scale $scale> <-damp $damp>\n"
			<< "or\n"
			<< "ResponseSpectrumAnalysis $dir -Tn $TnValues -fn $fnValues -Sa $SaValues <-scale $scale> <-damp $damp>\n"
			// Ladruno WP-207
			<< "with <-increment <-sign $sign> <-revert>> to impose each mode about the committed state,\n"
			<< "and ResponseSpectrumAnalysis -revert to return to it.\n"
			// Ladruno WP-212
			<< "<-requireMass $ratio> errors when the modes carry less than $ratio of the free mass.\n"
			"Error: at least 2 arguments should be provided.\n";
		return -1;
	}

	// search for -Tn -Sa, if both found we use them as lists
	// otherwise we fallback to the old implementation of the timeSeries
	bool found_Tn = false;
	bool found_Sa = false;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* arg = OPS_GetString();
		if (!found_Tn && (strcmp(arg, "-Tn") == 0 || strcmp(arg, "-fn") == 0))
			found_Tn = true;
		if (!found_Sa && strcmp(arg, "-Sa") == 0)
			found_Sa = true;
	}
	bool use_lists = found_Tn && found_Sa;
	if (found_Tn && !found_Sa) {
		opserr << "ResponseSpectrumAnalysis Error: found -Tn or -fn but not -Sa, please specify both of them or use a timeSeries\n";
		return -1;
	}
	if (found_Sa && !found_Tn) {
		opserr << "ResponseSpectrumAnalysis Error: found -Sa but not -Tn or -fn, please specify both of them or use a timeSeries\n";
		return -1;
	}
	OPS_ResetCurrentInputArg(-nargs);

	// num data
	int numData = 1;

	// get time series
	if (!use_lists) {
		int tstag;
		if (OPS_GetInt(&numData, &tstag) < 0) {
			opserr << "ResponseSpectrumAnalysis Error: Failed to get timeSeries tag.\n";
			return -1;
		}
		ts = OPS_getTimeSeries(tstag);
		if (ts == nullptr) {
			opserr << "ResponseSpectrumAnalysis Error: Failed to get timeSeries with tag = " << tstag << ".\n";
			return -1;
		}
	}

	// get direction
	if (OPS_GetInt(&numData, &dir) < 0) {
		opserr << "ResponseSpectrumAnalysis Error: Failed to get direction.\n";
		return -1;
	}
	if (dir < 1 || dir > ndf) {
		opserr << "ResponseSpectrumAnalysis Error: provided direction (" << dir << ") should be in the range 1-" << ndf << ".\n";
		return -1;
	}

	// parse optional data
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-scale") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDouble(&numData, &scale) < 0) {
					opserr << "ResponseSpectrumAnalysis Error: Failed to get scale factor.\n";
					return -1;
				}
			}
			else {
				opserr << "ResponseSpectrumAnalysis Error: scale factor requested but not provided.\n";
				return -1;
			}
		}
		else if (strcmp(value, "-mode") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetInt(&numData, &mode_id) < 0) {
					opserr << "ResponseSpectrumAnalysis Error: Failed to get the mode_id.\n";
					return -1;
				}
				--mode_id; // make it 0-based
				single_mode = true;
			}
			else {
				opserr << "ResponseSpectrumAnalysis Error: mode_id requested but not provided.\n";
				return -1;
			}
		}
		// Ladruno ADR44 P1b: opt-in modal combination stage
		else if (strcmp(value, "-combine") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "ResponseSpectrumAnalysis Error: -combine requires a rule "
				          "(SRSS | CQC | ABS | TenPercent).\n";
				return -1;
			}
			const char* rname = OPS_GetString();
			combine_rule = LadrunoModalCombination::parseRule(rname);
			if (combine_rule == LadrunoModalCombination::RULE_NONE) {
				opserr << "ResponseSpectrumAnalysis Error: unknown -combine rule '"
				       << rname << "' (use SRSS | CQC | ABS | TenPercent).\n";
				return -1;
			}
		}
		else if (strcmp(value, "-damp") == 0) {
			double xi_one;
			if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &xi_one) < 0) {
				opserr << "ResponseSpectrumAnalysis Error: -damp requires a value.\n";
				return -1;
			}
			xi_list.assign(1, xi_one); // single value broadcast to all modes
		}
		// Ladruno WP-207: -increment [-sign +1|-1] [-revert]
		else if (strcmp(value, "-increment") == 0) {
			increment = true;
		}
		else if (strcmp(value, "-sign") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &sign) < 0) {
				opserr << "ResponseSpectrumAnalysis Error: -sign requires a value (+1 or -1).\n";
				return -1;
			}
			if (sign != 1.0 && sign != -1.0) {
				opserr << "ResponseSpectrumAnalysis Error: -sign must be +1 or -1 (got " << sign << ").\n";
				return -1;
			}
			sign_given = true;
		}
		else if (strcmp(value, "-revert") == 0) {
			revert_at_end = true;
		}
		// Ladruno WP-212: -requireMass $ratio, in (0, 1]
		else if (strcmp(value, "-requireMass") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &require_mass) < 0) {
				opserr << "ResponseSpectrumAnalysis Error: -requireMass requires a ratio in (0, 1].\n";
				return -1;
			}
			if (!(require_mass > 0.0 && require_mass <= 1.0)) {
				opserr << "ResponseSpectrumAnalysis Error: -requireMass must be in (0, 1] (got " << require_mass << ").\n";
				return -1;
			}
		}
		else if (strcmp(value, "-modalDamp") == 0 || strcmp(value, "-modalDamping") == 0) {
			xi_list.clear();
			while (OPS_GetNumRemainingInputArgs() > 0) {
				double item; int one = 1;
				int rem_before = OPS_GetNumRemainingInputArgs();
				if (OPS_GetDoubleInput(&one, &item) < 0) {
					if (OPS_GetNumRemainingInputArgs() < rem_before)
						OPS_ResetCurrentInputArg(-1);
					break;
				}
				xi_list.push_back(item);
			}
			if (xi_list.empty()) {
				opserr << "ResponseSpectrumAnalysis Error: -modalDamp requires >=1 ratio.\n";
				return -1;
			}
		}
		else if (strcmp(value, "-Tn") == 0 || strcmp(value, "-fn") == 0) {
			// first try expanded list like {*}$the_list,
			// also used in python like *the_list
			Tn.clear();
			while (OPS_GetNumRemainingInputArgs() > 0) {
				double item;
				auto old_num_rem = OPS_GetNumRemainingInputArgs();
				if (OPS_GetDoubleInput(&numData, &item) < 0) {
					auto new_num_rem = OPS_GetNumRemainingInputArgs();
					if (new_num_rem < old_num_rem)
						OPS_ResetCurrentInputArg(-1);
					break;
				}
				if (strcmp(value, "-fn") == 0 && item != 0) {
					item = 1.0 / item;
				}
				Tn.push_back(item);
			}
			// try Tcl list (it's a string after all...)
			if (Tn.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
				std::string list_string = OPS_GetString();
				if (!string_to_list_of_doubles(list_string, ' ', Tn)) {
					opserr << "ResponseSpectrumAnalysis Error: cannot parse the Tn list.\n";
					return -1;
				}
			}
		}
		else if (strcmp(value, "-Sa") == 0) {
			// first try expanded list like {*}$the_list,
			// also used in python like *the_list
			Sa.clear();
			while (OPS_GetNumRemainingInputArgs() > 0) {
				double item;
				auto old_num_rem = OPS_GetNumRemainingInputArgs();
				if (OPS_GetDoubleInput(&numData, &item) < 0) {
					auto new_num_rem = OPS_GetNumRemainingInputArgs();
					if (new_num_rem < old_num_rem)
						OPS_ResetCurrentInputArg(-1);
					break;
				}
				Sa.push_back(item);
			}
			// try Tcl list (it's a string after all...)
			if (Sa.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
				std::string list_string = OPS_GetString();
				if (!string_to_list_of_doubles(list_string, ' ', Sa)) {
					opserr << "ResponseSpectrumAnalysis Error: cannot parse the Sa list.\n";
					return -1;
				}
			}
		}
	}

	// check Tn and Sa vectors
	if (use_lists) {
		if (Tn.size() != Sa.size()) {
			opserr << "ResponseSpectrumAnalysis Error: Sa and Tn lists must have the same length\n";
			opserr << (int)Tn.size() << " != " << (int)Sa.size() << "\n";
			return -1;
		}
		if (Tn.size() == 0) {
			opserr << "ResponseSpectrumAnalysis Error: Sa and Tn lists cannot be empty\n";
			return -1;
		}
		for (std::size_t i = 0; i < Tn.size(); ++i) {
			if (Tn[i] < 0.0) {
				opserr << "ResponseSpectrumAnalysis Error: Tn values must be positive (found " << Tn[i] << ")\n";
				return -1;
			}
			if (i > 0 && Tn[i] <= Tn[i - 1]) {
				opserr << "ResponseSpectrumAnalysis Error: Tn values must be monotonically increasing (found " << Tn[i] << " after " << Tn[i - 1] << ")\n";
				return -1;
			}
		}
		for (double item : Sa) {
			if (item < 0.0) {
				opserr << "ResponseSpectrumAnalysis Error: Sa values must be positive (found " << item << ")\n";
				return -1;
			}
		}
	}

	// ok, create the response spectrum analysis and run it here...
	// no need to store it
	// Ladruno WP-207: -sign / -revert only make sense about a committed state
	if (!increment && (sign_given || revert_at_end)) {
		opserr << "ResponseSpectrumAnalysis Error: -sign and -revert require -increment.\n";
		return -1;
	}
	if (increment && combine_rule != LadrunoModalCombination::RULE_NONE) {
		opserr << "ResponseSpectrumAnalysis Error: -increment and -combine are "
		          "mutually exclusive (-increment imposes one mode at a time).\n";
		return -1;
	}

	ResponseSpectrumAnalysis rsa(theAnalysisModel, ts, Tn, Sa, dir, scale);
	if (increment)
		rsa.setIncrement(sign, revert_at_end); // Ladruno WP-207
	rsa.setRequiredMassRatio(require_mass); // Ladruno WP-212
	int result;
	// Ladruno ADR44 P1b: if -combine given, run the combination stage (one
	// committed field); otherwise the stock per-mode behavior is byte-identical.
	if (combine_rule != LadrunoModalCombination::RULE_NONE) {
		if (single_mode) {
			opserr << "ResponseSpectrumAnalysis Error: -combine and -mode are "
			          "mutually exclusive (combination needs all modes).\n";
			return -1;
		}
		rsa.setCombination(combine_rule, xi_list);
		result = rsa.analyzeCombined();
	}
	else if (single_mode)
		result = rsa.analyze(mode_id);
	else
		result = rsa.analyze();

	return result;
}

ResponseSpectrumAnalysis::ResponseSpectrumAnalysis(
	AnalysisModel* theModel,
	TimeSeries* theFunction,
	const std::vector<double>& Tn,
	const std::vector<double>& Sa,
	int theDirection,
	double scale
)
	: m_model(theModel)
	, m_function(theFunction)
	, m_Tn(Tn)
	, m_Sa(Sa)
	, m_direction(theDirection)
	, m_scale(scale)
	, m_current_mode(0)
{

}

ResponseSpectrumAnalysis::~ResponseSpectrumAnalysis()
{
}

int ResponseSpectrumAnalysis::analyze()
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - failed to get modal properties" << endln;
		return -1;
	}

	// size info
	int num_eigen = domain->getEigenvalues().Size();

	// check consistency
	int error_code;
	error_code = check();
	if (error_code < 0) return error_code;

	// Ladruno WP-212: mass-shortfall report and guard
	error_code = checkMass(mp);
	if (error_code < 0) return error_code;

	// loop over all required eigen-modes, compute the modal displacement
	// and save the results.
	// we just compute modal displacements without doing any (SRSS, CQC, etc..)
	// modal combination otherwise derived results cannot be computed.
	// for each mode, this analysis produces a new analysis step.
	// modal combination of displacements (or any derived results)
	// it's up to the user.
	for (m_current_mode = 0; m_current_mode < num_eigen; ++m_current_mode)
	{
		// init the new step
		error_code = beginMode();
		if (error_code < 0) return error_code;

		// compute modal acceleration for this mode using the 
		// provided response spectrum function (time series)
		error_code = solveMode();
		if (error_code < 0) return error_code;

		// done with the current step.
		// here the modal displacements will be recorded (and all other results
		// if requested via recorders...)
		error_code = endMode();
		if (error_code < 0) return error_code;
	}

	return 0;
}

int ResponseSpectrumAnalysis::analyze(int mode_id)
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	//const DomainModalProperties& mp = domain->getModalProperties();
	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze(mode_id) - failed to get modal properties" << endln;
		return -1;
	}

	// size info
	int num_eigen = domain->getEigenvalues().Size();
	if (mode_id < 0 || mode_id >= num_eigen) {
		opserr << "ResponseSpectrumAnalysis::analyze(mode_id) - The provided mode_id (" << mode_id + 1 << ") is out of range (1, " << num_eigen << ")\n";
		return -1;
	}
	m_current_mode = mode_id;

	// check consistency
	int error_code;
	error_code = check();
	if (error_code < 0) return error_code;

	// Ladruno WP-212: mass-shortfall report and guard
	error_code = checkMass(mp);
	if (error_code < 0) return error_code;

	// init the new step
	error_code = beginMode();
	if (error_code < 0) return error_code;

	// compute modal acceleration for this mode using the 
	// provided response spectrum function (time series)
	error_code = solveMode();
	if (error_code < 0) return error_code;

	// done with the current step.
	// here the modal displacements will be recorded (and all other results
	// if requested via recorders...)
	error_code = endMode();
	if (error_code < 0) return error_code;

	return 0;
}

int ResponseSpectrumAnalysis::check()
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	//const DomainModalProperties& mp = domain->getModalProperties();
	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::check() - failed to get modal properties" << endln;
		return -1;
	}

	// number of eigen-modes
	int num_eigen = domain->getEigenvalues().Size();
	if (num_eigen < 1) {
		opserr << "ResponseSpectrumAnalysis::check() - No Eigenvalue provided.\n";
		return -1;
	}

	// check consistency
	auto check_eigen = [&mp, domain]() -> bool {
		const Vector& ev = domain->getEigenvalues();
		if (ev.Size() != mp.eigenvalues().Size())
			return false;
		double tol = std::max(1.0e-15, 1.0e-12 * ev.Norm());
		for (int i = 0; i < ev.Size(); ++i) {
			double a = ev(i);
			double b = mp.eigenvalues()(i);
			if (std::abs(a - b) > tol)
				return false;
		}
		return true;
	};
	if (!check_eigen()) {
		opserr << "ResponseSpectrumAnalysis::check() - Eigenvalues stored in DomainModalProperties are not equal to the eigenvalues in the model.\n"
			"Make sure to call the 'modalProperties' command\n"
			"after the 'eigen' command, and right before the 'responseSpectrum' command.\n";
		return -1;
	}

	return 0;
}

int ResponseSpectrumAnalysis::beginMode()
{
	// Ladruno WP-207: -increment starts every mode from the committed state
	// (drops the trial state a previous -increment mode left in place)
	if (m_increment && m_model->revertDomainToLastCommit() < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed in revertDomainToLastCommit"
			" at mode " << m_current_mode << "\n";
		return -1;
	}

	// new step... (do nothing)
	if (m_model->analysisStep() < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed"
			" at mode " << m_current_mode << "\n";
		return -1;
	}

	return 0;
}

int ResponseSpectrumAnalysis::endMode()
{
	// update domain
	if (m_model->updateDomain() < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed in updateDomain"
			" at mode " << m_current_mode << "\n";
		return -1;
	}

	// Ladruno WP-207: -increment records the trial state but never commits it,
	// so the committed (e.g. gravity + loadConst) state is preserved
	if (m_increment) {
		m_model->getDomainPtr()->record();
		if (m_revert_at_end && m_model->revertDomainToLastCommit() < 0) {
			opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed in revertDomainToLastCommit"
				" at mode " << m_current_mode << "\n";
			return -1;
		}
		return 0;
	}

	// commit domain
	if (m_model->commitDomain() < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed in commitDomain"
			" at mode " << m_current_mode << "\n";
		return -1;
	}

	return 0;
}

int ResponseSpectrumAnalysis::solveMode()
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	//const DomainModalProperties& mp = domain->getModalProperties();
	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::solveMode() - failed to get modal properties" << endln;
		return -1;
	}

	// size info
	int ndf = mp.totalMass().Size();

	// excited DOF.
	// Note: now we assume that a RS acts along one of the global directions, so we need
	// to consider only the associated column of the modal participation factors.
	// in future versions we can implement a general direction vector.
	int exdof = m_direction - 1; // make it 0-based

	// compute modal acceleration for this mode using the 
	// provided response spectrum function (time series)
	double lambda = mp.eigenvalues()(m_current_mode);
	double omega = std::sqrt(lambda);
	double freq = omega / 2.0 / M_PI;
	double period = 1.0 / freq;
	// Ladruno WP-213: a residualModes -period T pseudo-mode reads Sa(T); lambda
	// (and so u = phi Gamma Sa / lambda) stays its own. No override => stock.
	const double sa_period = LadrunoModeKinds::saPeriod(domain, m_current_mode);
	double mga = getSa(sa_period > 0.0 ? sa_period : period);
	double Vscale = mp.eigenVectorScaleFactors()(m_current_mode);
	double MPF = mp.modalParticipationFactors()(m_current_mode, exdof);

	// loop over all nodes and compute the modal displacements
	Node* node;
	NodeIter& theNodes = domain->getNodes();
	while ((node = theNodes()) != 0) {

		// get the nodal eigenvector, according to the ndf of modal properties
		// and its scaling
		const Matrix& node_evec = node->getEigenvectors();
		int node_ndf = node_evec.noRows();
		// Ladruno WP-207: -increment imposes u_k about the committed displacement
		const Vector* u_commit = m_increment ? &node->getDisp() : nullptr;

		// for each DOF...
		for (int i = 0; i < std::min(node_ndf, ndf); ++i) {
			// exclude any dof >= node_ndf (if ndf > node_ndf... in case node is U-only)
			// exclude any dof >= ndf (if node_ndf > ndf... in case the node has more DOFs than U and R)

			// we also need to exclude pressure dofs... easy in 3D because node_ndf is 4,
			// but in 2D it is 3, as in U-R...
			if (ndf == 6 && node_ndf == 4 && i == 3)
				continue;

			// eigenvector
			double V = node_evec(i, m_current_mode) * Vscale;

			// compute modal displacements for the i-th DOF
			double u_modal = V * MPF * mga / lambda;

			// Ladruno WP-207: -increment => trial = committed + sign * u_modal
			if (u_commit)
				u_modal = (*u_commit)(i) + m_sign * u_modal;

			// save this displacement at the i-th dof as new trial displacement
			node->setTrialDisp(u_modal, i);
		}

	}

	return 0;
}

// Ladruno ADR44 P1b: opt-in modal combination.
void ResponseSpectrumAnalysis::setCombination(int rule, const std::vector<double>& xi)
{
	m_combine = rule;
	m_xi = xi;
}

// Ladruno WP-207: opt-in incremental mode about the committed state.
void ResponseSpectrumAnalysis::setIncrement(double sign, bool revertAtEnd)
{
	m_increment = true;
	m_sign = sign;
	m_revert_at_end = revertAtEnd;
}

// Ladruno WP-212: -requireMass.
void ResponseSpectrumAnalysis::setRequiredMassRatio(double ratio)
{
	m_require_mass = ratio;
}

// Ladruno WP-212: ASCE 7-22 12.9.1.1 wants the modes to carry 100% of the mass
// (90% per horizontal direction as the exception); NEC-SE-DS 6.2.2e wants 90%.
// Report the cumulative effective-mass ratio of the excited direction, warn
// below 90%, and refuse a direction with no free or no participating mass (the
// signature of the WP-211 diaphragm-mass bug: every Gamma zero, RSA field zero).
// Only stderr is touched: the computed field is the stock one.
int ResponseSpectrumAnalysis::checkMass(const DomainModalProperties& mp) const
{
	const int exdof = m_direction - 1;
	const int num_eigen = mp.eigenvalues().Size();
	if (num_eigen < 1 || exdof < 0 || exdof >= mp.totalFreeMass().Size())
		return 0; // check() and the parser report these

	const double free_mass = mp.totalFreeMass()(exdof);
	if (!(free_mass > 0.0)) {
		opserr << "ResponseSpectrumAnalysis Error: no free mass in direction " << m_direction
			<< ": the spectrum excites nothing. Check the mass assignment, the constraints and the direction.\n";
		return -1;
	}
	const double ratio = mp.modalParticipationMassRatiosCumulative()(num_eigen - 1, exdof);

	char pct[32];
	std::snprintf(pct, sizeof(pct), "%.2f%%", 100.0 * ratio);

	if (!(ratio > 1.0e-12)) {
		opserr << "ResponseSpectrumAnalysis Error: no participating mass in direction " << m_direction
			<< " (every participation factor of the " << num_eigen << " modes is zero, ratio " << pct
			<< "): the response would be zero. Check that the modes see the mass (constraints, diaphragms).\n";
		return -1;
	}

	if (std::getenv("LADRUNO_OPENSEES_QUIET") == nullptr)
		opserr << "ResponseSpectrumAnalysis: cumulative effective mass in direction " << m_direction
			<< " = " << pct << " of the free mass (" << num_eigen << " modes).\n";

	// relative 1e-12 slack: a complete basis sums to 1 only to round-off
	if (m_require_mass > 0.0 && ratio < m_require_mass * (1.0 - 1.0e-12)) {
		char req[32];
		std::snprintf(req, sizeof(req), "%.2f%%", 100.0 * m_require_mass);
		opserr << "ResponseSpectrumAnalysis Error: the modes carry " << pct << " of the free mass in direction "
			<< m_direction << ", below -requireMass " << req << ". Extract more modes.\n";
		return -1;
	}
	if (ratio < 0.9)
		opserr << "ResponseSpectrumAnalysis WARNING: the modes carry only " << pct << " of the free mass in direction "
			<< m_direction << ", below 90% (ASCE 7-22 12.9.1.1, NEC-SE-DS 6.2.2e). Extract more modes.\n";

	return 0;
}

// Ladruno ADR44 P1b: combine the per-mode peak nodal DISPLACEMENTS into one
// design field (CQC/SRSS/ABS/TenPercent) and commit a single state, instead of
// the per-mode states the stock analyze() produces. Combination is per-quantity
// and nonlinear: this combines displacements only (see LadrunoModalCombination.h).
int ResponseSpectrumAnalysis::analyzeCombined()
{
	Domain* domain = m_model->getDomainPtr();

	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::analyzeCombined() - failed to get modal properties\n";
		return -1;
	}
	int error_code = check();
	if (error_code < 0) return error_code;

	// Ladruno WP-212: mass-shortfall report and guard
	error_code = checkMass(mp);
	if (error_code < 0) return error_code;

	// Ladruno WP-212: what is committed below is the combined DISPLACEMENT field;
	// any force, reaction or drift recovered from it is not the combination of
	// that quantity (ASCE 7-22 12.9.1.3). Say so once per call.
	opserr << "ResponseSpectrumAnalysis WARNING: -combine combines nodal displacements only; "
	          "forces, reactions and drifts recovered from them are not modally combined. "
	          "Use -mode k per mode and combine each quantity.\n";

	const int num_eigen = domain->getEigenvalues().Size();
	const int ndf = mp.totalMass().Size();
	const int exdof = m_direction - 1;

	// per-mode scalars over ALL extracted modes
	std::vector<double> w(num_eigen), T(num_eigen), MPF(num_eigen),
	                    Vscale(num_eigen), Sa(num_eigen), lambda(num_eigen),
	                    xi(num_eigen);
	// Ladruno WP-213: per-mode Sa-period overrides (residualModes -period T);
	// T[a] and w[a] (TenPercent closeness, CQC rho) stay the modes' own
	const std::vector<double> sa_period = LadrunoModeKinds::saPeriods(domain);
	for (int a = 0; a < num_eigen; ++a) {
		lambda[a] = mp.eigenvalues()(a);
		const double omega = (lambda[a] > 0.0) ? std::sqrt(lambda[a]) : 0.0;
		w[a] = omega;
		const double period = (omega > 0.0) ? (2.0 * M_PI / omega) : 0.0;
		T[a] = period;
		Sa[a] = getSa((a < (int)sa_period.size() && sa_period[a] > 0.0)
		              ? sa_period[a] : period); // Ladruno WP-213
		Vscale[a] = mp.eigenVectorScaleFactors()(a);
		MPF[a] = mp.modalParticipationFactors()(a, exdof);
		if (m_xi.size() == 1) xi[a] = m_xi[0];               // broadcast
		else if (a < static_cast<int>(m_xi.size())) xi[a] = m_xi[a];
		else xi[a] = 0.0;
	}

	// Non-oscillatory modes (lambda <= 0: rigid, or spurious negative from a
	// near-singular / geometric-stiffness eigen solve) do not contribute to a
	// response-spectrum combination and would divide by ~0 in R_a — exclude them
	// (R_a = 0, w_a = 0 already) and warn once.
	for (int a = 0; a < num_eigen; ++a)
		if (lambda[a] <= 0.0) {
			opserr << "ResponseSpectrumAnalysis - mode " << (a + 1)
			       << " has a non-positive eigenvalue (" << lambda[a]
			       << "); excluded from the -combine field.\n";
		}
	if (m_combine == LadrunoModalCombination::RULE_CQC) {
		// CQC needs modal damping to be meaningful (with zero damping every
		// rho_ij collapses to 0 and CQC degenerates to SRSS). An undamped SUBSET
		// is still well defined — those modes just get rho=0 (uncorrelated) and
		// the rhoCQC 0/0 edge (coincident undamped pair) is guarded in the kernel.
		bool any_damp = false;
		for (int a = 0; a < num_eigen; ++a)
			if (w[a] > 0.0 && xi[a] > 0.0) { any_damp = true; break; }
		if (!any_damp) {
			opserr << "ResponseSpectrumAnalysis - CQC combination requires modal "
			          "damping (-damp or -modalDamp); none was given.\n";
			return -1;
		}
	}

	error_code = beginMode();
	if (error_code < 0) return error_code;

	Node* node;
	NodeIter& theNodes = domain->getNodes();
	std::vector<double> R(num_eigen);
	while ((node = theNodes()) != 0) {
		const Matrix& node_evec = node->getEigenvectors();
		const int node_ndf = node_evec.noRows();
		if (node_ndf < 1) continue;
		const int nd = std::min(node_ndf, ndf);
		for (int i = 0; i < nd; ++i) {
			if (ndf == 6 && node_ndf == 4 && i == 3) // 3D U-P pressure DOF
				continue;
			for (int a = 0; a < num_eigen; ++a) {
				const double V = node_evec(i, a) * Vscale[a];
				// exclude non-oscillatory modes (matches the w[a] clamp); lambda>0
				// avoids the div-by-~0 blow-up for a rigid/spurious mode
				R[a] = (lambda[a] > 0.0) ? V * MPF[a] * Sa[a] / lambda[a] : 0.0;
			}
			const double u_comb = LadrunoModalCombination::combine(
				static_cast<LadrunoModalCombination::Rule>(m_combine),
				R, w, xi, T);
			node->setTrialDisp(u_comb, i);
		}
	}

	return endMode();
}

// Ladruno WP-202: -scale was parsed and stored in m_scale but never applied,
// so every spectral ordinate ignored it. Apply it once, here, so the per-mode
// (analyze) and the combined (analyzeCombined) paths both honour it.
double ResponseSpectrumAnalysis::getSa(double T) const
{
	return m_scale * getSaUnscaled(T);
}

double ResponseSpectrumAnalysis::getSaUnscaled(double T) const
{
	// use the time series if provided
	if (m_function)
		return m_function->getFactor(T);
	// otherwise use the vectors
	std::size_t n = m_Tn.size();
	// check empty
	if (n < 1)
		return 0.0;
	// check constant
	if (n == 1)
		return m_Sa[0];
	// check bounds
	if (T <= m_Tn.front())
		return m_Sa.front();
	if (T >= m_Tn.back())
		return m_Sa.back();
	// interpolate
	for (std::size_t i = 1; i < n; ++i) {
		double t1 = m_Tn[i];
		if (T <= t1) {
			double t0 = m_Tn[i - 1];
			double s1 = m_Sa[i];
			double s0 = m_Sa[i - 1];
			double dT = t1 - t0;
			double dS = s1 - s0;
			double scale = dT > 0.0 ? (T - t0) / dT : 0.5;
			return s0 + scale * dS;
		}
	}
	return m_Sa.back();
}

