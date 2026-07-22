#include "LadrunoCMSEigenSOE.h"
#include "LadrunoCMSEigenSolver.h"

#include <Graph.h>
#include <ID.h>
#include <Matrix.h>
#include <OPS_Stream.h>
#include <StandardStream.h>
#include <Vertex.h>

#include <mpi.h>

#include <algorithm>
#include <iostream>
#include <string>

StandardStream cmsTestErrorStream;
OPS_Stream *opserrPtr = &cmsTestErrorStream;

namespace {

int failures = 0;

void require(bool condition, const std::string &message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

Graph *makeChainGraph(int order)
{
    Graph *graph = new Graph(order);
    for (int equation = 0; equation < order; ++equation)
        graph->addVertex(new Vertex(equation, equation), false);
    for (int equation = 0; equation + 1 < order; ++equation)
        graph->addEdge(equation, equation + 1);
    return graph;
}

void addChainContributions(
    LadrunoCMSEigenSOE &soe,
    int order,
    int rankWithMismatch = -1,
    int worldRank = 0)
{
    Matrix stiffness(2, 2);
    stiffness(0, 0) = 2.0;
    stiffness(0, 1) = -1.0;
    stiffness(1, 0) = -1.0;
    stiffness(1, 1) = 2.0;
    Matrix mass(2, 2);
    mass(0, 0) = 1.0;
    mass(0, 1) = 0.1;
    mass(1, 0) = 0.1;
    mass(1, 1) = 1.0;
    for (int element = 0; element + 1 < order; ++element) {
        ID equations(2);
        equations(0) = element;
        equations(1) = element + 1;
        require(soe.addA(stiffness, equations, 1.0) == 0,
                "stiffness contribution was rejected");
        const double factor = rankWithMismatch == worldRank && element == order - 2
            ? 1.001
            : 1.0;
        require(soe.addM(mass, equations, factor) == 0,
                "mass contribution was rejected");
    }
}

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const int level1 = size >= 4 && size % 2 == 0 ? 2 : 1;
    const int level2 = size / level1;
    constexpr int order = 20;

    ladruno_cms::Options options;
    options.hierarchy = ladruno_cms::HierarchyMode::Logical;
    options.level1 = level1;
    options.level2 = level2;
    options.modesLevel2 = 6;
    options.modesLevel1 = 8;
    options.verifyAssembly = ladruno_cms::AssemblyVerification::Signature;
    LadrunoCMSEigenSolver *solver = new LadrunoCMSEigenSolver();
    LadrunoCMSEigenSOE *soe = new LadrunoCMSEigenSOE(*solver, options);
    Graph *graph = makeChainGraph(order);
    require(soe->setSize(*graph) == 0, "two-stage METIS topology failed");
    require(soe->getWorldSize() == size && soe->getWorldRank() == rank,
            "SOE stored an incorrect MPI identity");
    require(soe->getNumberOfCoarseSubdomains() == level1 &&
            soe->getFineSubdomainsPerCoarse() == level2,
            "logical hierarchy dimensions changed during setSize");
    require(soe->getFineLabels().size() == static_cast<std::size_t>(order),
            "fine owner map has the wrong length");
    require(std::all_of(
                soe->getFineLabels().begin(), soe->getFineLabels().end(),
                [size](int owner) { return owner >= 0 && owner < size; }),
            "fine owner map contains an invalid rank");

    addChainContributions(*soe, order);
    int localOwned = static_cast<int>(soe->getStiffnessContributions().size());
    int globalOwned = 0;
    MPI_Allreduce(&localOwned, &globalOwned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    require(globalOwned == order - 1,
            "a stiffness contribution was lost or assigned more than once");
    std::string message;
    require(soe->verifyReplicatedAssembly(message) == 0,
            "identical replicated assembly was rejected: " + message);
    require(soe->solve(3, true, true) == 0,
            "public SOE solve did not execute the distributed hierarchy");
    require(solver->getDiagnostics().appliedT2 &&
            solver->getDiagnostics().appliedS2 &&
            solver->getDiagnostics().appliedT1 &&
            solver->getDiagnostics().appliedS1,
            "public solver skipped a mandatory hierarchy phase");
    require(solver->getNormalizedResiduals().size() == 3u &&
            *std::max_element(solver->getNormalizedResiduals().begin(),
                              solver->getNormalizedResiduals().end()) <=
                options.tolerance,
            "public solver accepted an excessive residual");
    require(soe->getEigenvalue(1) > 0.0 &&
            soe->getEigenvector(1).Size() == order,
            "public OpenSees eigen result contract is incomplete");

    if (size > 1) {
        soe->zeroA();
        soe->zeroM();
        addChainContributions(*soe, order, 1, rank);
        message.clear();
        require(soe->verifyReplicatedAssembly(message) < 0,
                "a rank-dependent numerical assembly mismatch was accepted");
        require(message.find("rank") != std::string::npos,
                "assembly mismatch did not identify a failing rank");
    }

    delete soe;
    ladruno_cms::Options automaticOptions = options;
    automaticOptions.hierarchy = ladruno_cms::HierarchyMode::Auto;
    automaticOptions.level1 = 0;
    automaticOptions.level2 = 0;
    LadrunoCMSEigenSolver *automaticSolver = new LadrunoCMSEigenSolver();
    LadrunoCMSEigenSOE *automaticSOE =
        new LadrunoCMSEigenSOE(*automaticSolver, automaticOptions);
    require(automaticSOE->setSize(*graph) == 0,
            "shared-memory auto topology failed");
    require(automaticSOE->getNumberOfCoarseSubdomains() == 1 &&
            automaticSOE->getFineSubdomainsPerCoarse() == size,
            "auto topology did not recognize the single shared-memory node");
    delete automaticSOE;
    delete graph;
    int globalFailures = 0;
    MPI_Allreduce(&failures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Finalize();
    if (globalFailures != 0)
        return 1;
    if (rank == 0)
        std::cout << "CMS topology and ownership checks passed\n";
    return 0;
}
