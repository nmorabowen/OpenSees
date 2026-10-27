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

#ifndef ResponseSpectrumAnalysis_h
#define ResponseSpectrumAnalysis_h

#include <vector>
class AnalysisModel;
class TimeSeries;
class DomainModalProperties; // Ladruno WP-212

class ResponseSpectrumAnalysis
{
public:
	ResponseSpectrumAnalysis(
		AnalysisModel* theModel,
		TimeSeries* theFunction,
		const std::vector<double>& Tn,
		const std::vector<double>& Sa,
		int theDirection,
		double scale
	);
	~ResponseSpectrumAnalysis();

public:
	int analyze();
	int analyze(int mode_id);

	// Ladruno ADR44 P1b: opt-in modal combination (CQC/SRSS/ABS/TenPercent).
	// rule = LadrunoModalCombination::Rule; xi = per-mode damping ratios (used by
	// CQC; a single value is broadcast to all modes). Committing ONE combined
	// nodal-displacement field instead of the per-mode states.
	void setCombination(int rule, const std::vector<double>& xi);
	int analyzeCombined();

	// Ladruno WP-207: opt-in incremental mode (-increment). Each mode imposes
	// trial = committed + sign*u_k, updates the domain and records, but does NOT
	// commit, so a preloaded (gravity + loadConst) state survives. The trial
	// state of the last mode is left in place for eleResponse; revertAtEnd
	// (-revert) returns to the committed state before the call returns.
	void setIncrement(double sign, bool revertAtEnd);

	// Ladruno WP-212: -requireMass r, in (0, 1]. Below r the mass-shortfall
	// warning of checkMass() becomes an error. 0 => warn only (default).
	void setRequiredMassRatio(double ratio);

	// Ladruno WP-216: public so responseSpectrumPerMode applies the same
	// modalProperties consistency check and the same spectral ordinate
	// (-scale included) as the per-mode path.
	int check();
	double getSa(double T) const;

private:
	// Ladruno WP-212: report the cumulative effective-mass ratio of the excited
	// direction once per call; warn below 90%, error with no participating mass
	// or below -requireMass.
	int checkMass(const DomainModalProperties& mp) const;
	int beginMode();
	int endMode();
	int solveMode();
	double getSaUnscaled(double T) const; // Ladruno WP-202

private:
	// the model
	AnalysisModel* m_model;
	// the response spectrum function
	TimeSeries* m_function;
	std::vector<double> m_Tn;
	std::vector<double> m_Sa;
	// the direction 1 to 3 (for 2D models) or 1 to 6 (for 3D models)
	int m_direction;
	// the scale factor for the computed displacement field
	double m_scale;
	// current mode
	int m_current_mode;
	// Ladruno ADR44 P1b: combination rule (-1 = none => stock per-mode behavior)
	// and per-mode damping ratios (for CQC).
	int m_combine = -1;
	std::vector<double> m_xi;
	// Ladruno WP-207: incremental mode (false => stock replace-and-commit).
	bool m_increment = false;
	double m_sign = 1.0;
	bool m_revert_at_end = false;
	// Ladruno WP-212: required cumulative mass ratio (0 => warn only).
	double m_require_mass = 0.0;
};

#endif
