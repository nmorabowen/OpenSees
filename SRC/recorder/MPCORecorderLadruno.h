/* ********************************************************************** **
**  MPCO_Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Phase-1 skeleton: registers `recorder mpcoLadruno`, writes a          **
**  schema-v1-valid .ladruno (INFO + MODEL_STAGE + MODEL/NODES). Node,    **
**  element, domain results + envelopes land in later phases.            **
**                                                                        **
**  All reusable machinery lives in namespace mpcol (MPCOL_*.h) to avoid  **
**  ODR clashes with the frozen MPCORecorder translation unit.            **
** ********************************************************************** */

#ifndef MPCORecorderLadruno_h
#define MPCORecorderLadruno_h

#include <Recorder.h>

class Domain;
class Channel;
class FEM_ObjectBroker;

class MPCORecorderLadruno : public Recorder
{
	friend void* OPS_MPCOLadrunoRecorder();

public:
	MPCORecorderLadruno();
	~MPCORecorderLadruno();

	int record(int commitTag, double timeStamp);
	virtual int restart(void);
	virtual int domainChanged(void);
	virtual int setDomain(Domain& theDomain);
	virtual int sendSelf(int commitTag, Channel& theChannel);
	virtual int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

private:
	int initialize();
	int writeModel();
	int writeModelNodes();
	int writeModelElements();
	int writeModelSets();
	int writeSections();

	int initNodeSources();
	int initElementSources();
	int initDomainSources();
	int recordResultsOnNodes();
	int recordResultsOnElements();
	int recordResultsOnDomain();

	// helpers (operate on the channel at the given index in m_data, so the
	// header never has to name private_data's nested channel types)
	void recordModeChannel(int node_channel_index);
	void writeElementColumnMap(int elem_channel_index);

	int clearSources();

	class private_data;
	private_data* m_data;
};

#endif // MPCORecorderLadruno_h
