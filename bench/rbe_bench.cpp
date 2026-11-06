// Benchmark CLI. Example: rbe_bench --scene pile --bodies 2000 --steps 600 --broadphase tree --csv out.csv
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "rbe/world.hpp"
#include "scenes.hpp"

using namespace rbe;

int main(int argc, char** argv) {
    std::string scene = "pile", csv, broad = "tree", gyro = "implicit";
    int bodies = 2000, steps = 600, measureFrom = 300;
    bool sleep = false;
    for (int i = 1; i < argc; ++i) {
        auto next = [&] { return i + 1 < argc ? std::string(argv[++i]) : std::string(); };
        if (!std::strcmp(argv[i], "--scene")) scene = next();
        else if (!std::strcmp(argv[i], "--bodies")) bodies = std::atoi(next().c_str());
        else if (!std::strcmp(argv[i], "--steps")) steps = std::atoi(next().c_str());
        else if (!std::strcmp(argv[i], "--measure-from")) measureFrom = std::atoi(next().c_str());
        else if (!std::strcmp(argv[i], "--broadphase")) broad = next();
        else if (!std::strcmp(argv[i], "--gyro")) gyro = next();
        else if (!std::strcmp(argv[i], "--sleep")) sleep = true;
        else if (!std::strcmp(argv[i], "--csv")) csv = next();
        else {
            std::fprintf(stderr, "unknown option %s\n", argv[i]);
            return 2;
        }
    }
    if (measureFrom >= steps) measureFrom = steps / 2;

    WorldSettings s;
    s.allowSleep = sleep;
    s.broadphase = broad == "brute" ? BroadphaseMode::BruteForce : BroadphaseMode::Tree;
    s.gyro = gyro == "explicit" ? GyroMode::Explicit : gyro == "none" ? GyroMode::None : GyroMode::Implicit;
    World w(s);
    if (!scenes::build(w, scene, bodies)) {
        std::fprintf(stderr, "unknown scene %s\n", scene.c_str());
        return 2;
    }

    FILE* out = csv.empty() ? nullptr : std::fopen(csv.c_str(), "w");
    if (out) std::fprintf(out, "step,total_ms,broadphase_ms,narrowphase_ms,solve_ms,pairs,manifolds,contacts,awake,energy\n");

    double sum = 0, sumBroad = 0, sumNarrow = 0, sumSolve = 0, worst = 0;
    int measured = 0;
    for (int i = 0; i < steps; ++i) {
        w.step(1.0f / 60.0f);
        const StepStats& st = w.stats();
        if (out)
            std::fprintf(out, "%d,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%.6f\n", i, st.msTotal, st.msBroadphase,
                         st.msNarrowphase, st.msSolve, st.pairs, st.manifolds, st.contactPoints, st.awakeBodies,
                         w.kineticEnergy());
        if (i >= measureFrom) {
            sum += st.msTotal;
            sumBroad += st.msBroadphase;
            sumNarrow += st.msNarrowphase;
            sumSolve += st.msSolve;
            if (st.msTotal > worst) worst = st.msTotal;
            ++measured;
        }
    }
    if (out) std::fclose(out);
    const StepStats& st = w.stats();
    std::printf("scene=%s bodies=%d broadphase=%s steps %d-%d: mean %.3f ms (broad %.3f, narrow %.3f, solve %.3f), worst %.3f ms, "
                "contacts %d, awake %d, hash %016llx\n",
                scene.c_str(), w.bodyCount(), broad.c_str(), measureFrom, steps, sum / measured, sumBroad / measured,
                sumNarrow / measured, sumSolve / measured, worst, st.contactPoints, st.awakeBodies,
                static_cast<unsigned long long>(w.stateHash()));
    return 0;
}
