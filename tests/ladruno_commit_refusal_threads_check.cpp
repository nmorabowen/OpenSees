// WP-224: the commit-refusal seam (SRC/material/LadrunoMaterialStatus.h) under threads.
// Pure C++ (no OpenSees build); compiled and run by tests/test_wp224_commit_refusal_threads.py.
// Fails on the pre-WP-224 header: a plain static int counter loses increments under concurrent
// ladrunoNoteCommitRefusal() calls (WP-192 threads the element commit loop), and a process-wide
// probe-scope depth lets a probe on one thread silence a real refusal on another.
#include "LadrunoMaterialStatus.h"
#include <cstdio>
#include <thread>
#include <vector>

int main()
{
    const int nThreads = 8, perThread = 200000;
    int worstLost = 0;
    for (int rep = 0; rep < 20; rep++) {
        ladrunoClearCommitRefusals();
        std::vector<std::thread> t;
        for (int k = 0; k < nThreads; k++)
            t.emplace_back([] { for (int i = 0; i < perThread; i++) ladrunoNoteCommitRefusal(); });
        for (auto &x : t) x.join();
        const int lost = nThreads * perThread - ladrunoPendingCommitRefusals();
        if (lost > worstLost) worstLost = lost;
    }

    // Thread A holds a probe scope; thread B declares a real refusal while it is open.
    bool suppressed = false;
    ladrunoClearCommitRefusals();
    std::thread a([&] {
        LadrunoProbeCommitScope probe;
        std::thread b([&] { ladrunoNoteCommitRefusal(); });
        b.join();
    });
    a.join();
    suppressed = (ladrunoPendingCommitRefusals() == 0);

    std::printf("lost increments (worst of 20 x %d x %d): %d; cross-thread probe suppression: %s\n",
                nThreads, perThread, worstLost, suppressed ? "YES" : "no");
    return (worstLost == 0 && !suppressed) ? 0 : 1;
}
