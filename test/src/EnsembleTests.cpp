#include "TestingUtility.h"

#include "pica/Ensemble.h"
#include "pica/utility/Utility.h"

#include <string>
#include <vector>

using namespace pica;
using std::string;
using std::vector;


class EnsembleTest : public BaseParticleFixture {
};


TEST_F(EnsembleTest, Constructor) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    Int3 expectedNumSystems = parameters.localGridSize + Int3(2, 2, 2);
    ASSERT_EQ_INT3(expectedNumSystems, ensemble.numSystems);
    ASSERT_FALSE(ensemble.inParticlesLoop);
    ASSERT_EQ(0, ensemble.getNumParticles());
}

bool eqEnsemble(Ensemble& a, Ensemble& b) {
    if (a.getNumParticles() != b.getNumParticles())
        return false;
    a.begin();
    b.begin();
    while (!a.end()) {
        if (!eqParticles(*a.currentParticle, *b.currentParticle))
            return false;
        a.next();
        b.next();
    }
    return true;
}

TEST_F(EnsembleTest, CopyConstructor) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int i = 0; i < 23; i++)
        ensemble.add(randomParticle());
    Ensemble copyEnsemble(ensemble);
    ASSERT_TRUE(eqEnsemble(ensemble, copyEnsemble));
}

TEST_F(EnsembleTest, Assignment) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    Ensemble copyEnsemble1(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    Ensemble copyEnsemble2(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int i = 0; i < 15; i++)
        ensemble.add(randomParticle());
    for (int i = 0; i < 22; i++)
        copyEnsemble1.add(randomParticle());
    copyEnsemble1 = ensemble;
    copyEnsemble2 = ensemble;
    ASSERT_TRUE(eqEnsemble(ensemble, copyEnsemble1));
    ASSERT_TRUE(eqEnsemble(ensemble, copyEnsemble2));
}

bool eqEnsembleParticleVector(Ensemble& ensemble, const vector<Particle>& particles) {
    if (particles.size() != ensemble.getNumParticles())
        return false;
    for (int j = 0; j < particles.size(); j++) {
        bool found = false;
        for (ensemble.begin(); !ensemble.end(); ensemble.next())
            if (eqParticles(particles[j], *ensemble.currentParticle)) {
                found = true;
                break;
            }
        if (!found)
            return false;
    }
    return true;
}

TEST_F(EnsembleTest, Add) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    vector<Particle> particles;
    for (int i = 0; i < 37; i++) {
        particles.push_back(randomParticle());
        ensemble.inParticlesLoop = false;
        ensemble.add(particles.back());
        ASSERT_TRUE(eqEnsembleParticleVector(ensemble, particles));
    }
}

TEST_F(EnsembleTest, AddPendingParticles) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    vector<Particle> particles;
    for (int i = 0; i < 37; i++) {
        ensemble.inParticlesLoop = true;
        particles.push_back(randomParticle());
        ensemble.add(particles.back());
        particles.push_back(randomParticle());
        ensemble.add(particles.back());
        ensemble.addPendingParticles();
        ASSERT_TRUE(eqEnsembleParticleVector(ensemble, particles));
    }
}

TEST_F(EnsembleTest, AddToCurrentSystem) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    Particle firstParticle = randomParticle();
    ensemble.add(firstParticle);
    for (int i = 0; i < 72; i++) {
        ensemble.begin();
        ensemble.inParticlesLoop = false;
        Particle newParticle = *ensemble.currentParticle;
        newParticle.setFactor(newParticle.getFactor() * (1.0 + (double)i / 50));
        ensemble.add(newParticle);
        // Check that the current particle stays the same even
        // in case of memory reallocation
        ASSERT_TRUE(eqParticles(firstParticle, *ensemble.currentParticle));
    }
}

TEST_F(EnsembleTest, AddParticle) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    vector<Particle> particles;
    for (int i = 0; i < 20; i++) {
        particles.push_back(randomParticle());
        ensemble.inParticlesLoop = false;
        ensemble.addParticle(particles.back());
        ASSERT_TRUE(eqEnsembleParticleVector(ensemble, particles));
    }
}

TEST_F(EnsembleTest, AddParticlePendingParticles) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    vector<Particle> particles;
    for (int i = 0; i < 37; i++) {
        ensemble.inParticlesLoop = true;
        particles.push_back(randomParticle());
        ensemble.addParticle(particles.back());
        particles.push_back(randomParticle());
        ensemble.addParticle(particles.back());
        ensemble.addPendingParticles();
        ASSERT_TRUE(eqEnsembleParticleVector(ensemble, particles));
    }
}

TEST_F(EnsembleTest, AddParticleToCurrentSystem) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    Particle firstParticle = randomParticle();
    ensemble.add(firstParticle);
    for (int i = 0; i < 72; i++) {
        ensemble.begin();
        ensemble.inParticlesLoop = false;
        Particle newParticle = *ensemble.currentParticle;
        newParticle.setFactor(newParticle.getFactor() * (1.0 + (double)i / 50));
        ensemble.addParticle(newParticle);
        // Check that the current particle stays the same even
        // in case of memory reallocation
        ASSERT_TRUE(eqParticles(firstParticle, *ensemble.currentParticle));
    }
}

TEST_F(EnsembleTest, IsMigratingParticle) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int i = 0; i < 51; i++) {
        Particle particle = randomParticle();
        Int3 cellIndex = ensemble.getCellIndex(particle);
        ASSERT_FALSE(ensemble.isMigratingParticle(particle, cellIndex));
        Int3 shifts[] = {Int3(1, 0, 0), Int3(0, 1, 0), Int3(0, 0, 1),
                         Int3(1, 1, 0), Int3(1, 0, 1), Int3(1, 1, 0),
                         parameters.localGridSize};
        for (int idx = 0; idx < sizeof(shifts) / sizeof(*shifts); idx++) {
            ASSERT_TRUE(ensemble.isMigratingParticle(particle, cellIndex + shifts[idx]));
            ASSERT_TRUE(ensemble.isMigratingParticle(particle, cellIndex - shifts[idx]));
        }
    }
}

TEST_F(EnsembleTest, Begin) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int idx = 0; idx < 41; idx++) {
        Particle particle = randomParticle();
        ensemble.add(particle);
    }
    ensemble.begin();
    Particle firstParticle = *ensemble.currentParticle;
    for (int i = 1; i < ensemble.getNumParticles(); i++) {
        ensemble.begin();
        for (int j = 0; j < i; j++)
            ensemble.next();
        ensemble.begin();
        ASSERT_TRUE(eqParticles(firstParticle, *ensemble.currentParticle));
    }
}

bool eqParticleVectorArray(const vector<Particle>& particleVector, const Particle* particleArray) {
    for (int i = 0; i < particleVector.size(); i++) {
        bool found = false;
        for (int j = 0; j < particleVector.size(); j++)
            if (eqParticles(particleVector[i], particleArray[j])) {
                found = true;
                break;
            }
        if (!found)
            return false;
    }
    return true;
}

TEST_F(EnsembleTest, Next) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    vector<Particle> particles;
    for (int idx = 0; idx < 9; idx++) {
        particles.push_back(randomParticle());
        ensemble.add(particles.back());
    }
    vector<Particle> traversedParticles;
    for (ensemble.begin(); !ensemble.end(); ensemble.next())
        traversedParticles.push_back(*ensemble.currentParticle);
    ASSERT_TRUE(eqParticleVectorArray(particles, ptr(traversedParticles)));
}

TEST_F(EnsembleTest, End) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    /// It might be a bug, probably begin() should be called in Ensemble constructor
    ensemble.begin();
    ASSERT_TRUE(ensemble.end());
    for (int idx = 0; idx < 31; idx++)
        ensemble.add(randomParticle());
    ensemble.begin();
    for (int i = 0; i < ensemble.getNumParticles(); i++) {
        ASSERT_FALSE(ensemble.end());
        ensemble.next();
    }
    ASSERT_TRUE(ensemble.end());
}

TEST_F(EnsembleTest, DeleteCurrentParticle) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int idx = 0; idx < 13; idx++)
        ensemble.add(randomParticle());
    int originalSize = ensemble.getNumParticles();
    for (int i = 0; i < originalSize; i++) {
        // copy ensemble into a vector
        vector<Particle> particles;
        for (ensemble.begin(); !ensemble.end(); ensemble.next())
            particles.push_back(*ensemble.currentParticle);
        // delete from the ensemble and the vector
        int deleteIdx = rand() % particles.size();
        ensemble.begin();
        for (int j = 0; j < deleteIdx; j++)
            ensemble.next();
        ensemble.deleteCurrentParticle();
        // Check that calling ensemble.next() will set ensemble.currentParticle
        // to the next particle
        if (deleteIdx < particles.size() - 1) {
            ensemble.next();
            ASSERT_TRUE(eqParticles(particles[deleteIdx + 1],
                                    *ensemble.currentParticle));
        }
        particles.erase(particles.begin() + deleteIdx);
        ASSERT_TRUE(eqEnsembleParticleVector(ensemble, particles));
    }
}

TEST_F(EnsembleTest, System) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int idx = 0; idx < 13; idx++)
        ensemble.add(randomParticle());
    for (int i = 0; i < ensemble.numSystems.x; i++)
    for (int j = 0; j < ensemble.numSystems.y; j++)
    for (int k = 0; k < ensemble.numSystems.z; k++) {
        ParticleSystem* systems[] = {ensemble.system(i, j, k), ensemble.system(Int3(i, j, k))};
        for (int systemIdx = 0; systemIdx < sizeof(systems) / sizeof(*systems); systemIdx++) {
            ParticleSystem* system = systems[systemIdx];
            for (system->begin(); !system->end(); system->next()) {
                Int3 cellIdx = ensemble.getCellIndex(*system->currentParticle);
                ASSERT_EQ_INT3(Int3(i, j, k), cellIdx);
            }
        }
    }
}

TEST_F(EnsembleTest, SystemConst) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int idx = 0; idx < 13; idx++)
        ensemble.add(randomParticle());
    const Ensemble& constEnsemble = ensemble;
    for (int i = 0; i < ensemble.numSystems.x; i++)
    for (int j = 0; j < ensemble.numSystems.y; j++)
    for (int k = 0; k < ensemble.numSystems.z; k++) {
        const ParticleSystem* constSystems[] = {constEnsemble.system(i, j, k),
                                                constEnsemble.system(Int3(i, j, k))};
        for (int systemIdx = 0; systemIdx < sizeof(constSystems) / 
                                            sizeof(*constSystems); systemIdx++) {
            const ParticleSystem* system = constSystems[systemIdx];
            for (int idx = 0; idx < system->size(); idx++) {
                Int3 cellIdx = constEnsemble.getCellIndex(system->raw()[idx]);
                ASSERT_EQ_INT3(Int3(i, j, k), cellIdx);
            }
        }
    }
}

TEST_F(EnsembleTest, Clear) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step,  3);
    for (int i = 0; i < 17; i++)
        ensemble.add(randomParticle());
    ensemble.clear();
    ASSERT_EQ(0, ensemble.getNumParticles());
}

TEST_F(EnsembleTest, GetNumAllocatedParticles) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int i = 0; i < 32; i++) {
        ensemble.add(randomParticle());
        ASSERT_LE(i + 1, ensemble.getNumAllocatedParticles());
    }
}

TEST_F(EnsembleTest, GetNumParticles) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int i = 0; i < 34; i++) {
        ensemble.add(randomParticle());
        ASSERT_EQ(i + 1, ensemble.getNumParticles());
    }
}

TEST_F(EnsembleTest, ShrinkToFit) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int i = 0; i < 71; i++) {
        ensemble.add(randomParticle());
        ensemble.shrinkToFit();
		// ShrinkToFit does not guarantee to decrease the amount of allocated memory
		// all the way to 0 extra memory, therefore, compare for <=
        ASSERT_LE(ensemble.getNumParticles(), ensemble.getNumAllocatedParticles());
    }
}

TEST_F(EnsembleTest, TryShrink) {
    // Here we can only check that memory overuse does not increase after tryShrink()
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    for (int i = 0; i < 66; i++) {
        ensemble.add(randomParticle());
        double beforeMemoryOveruse = (double)ensemble.getNumAllocatedParticles() /
                                     (double)ensemble.getNumParticles();
        ensemble.tryShrink();
        double afterMemoryOveruse = (double)ensemble.getNumAllocatedParticles() /
                                    (double)ensemble.getNumParticles();
        ASSERT_LE(afterMemoryOveruse, beforeMemoryOveruse);
    }
}

TEST_F(EnsembleTest, GetCellIndex) {
    // Here we can only check that memory overuse does not increase after tryShrink()
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    FP3 localMin = parameters.localMin;
    FP3 localMax = parameters.localMax;
    Int3 gridSize = parameters.localGridSize;
    FP3 step = parameters.step;
    FP3 coords[] = {localMin + step * 0.05, localMin + step * 0.9,
                    localMin + step * 1.1, localMin + step * FP3(1.2, 5.8, 4.3),
                    localMax - step * 1.2, localMax - step * 0.1};
    Int3 cellIdx[] = {Int3(1, 1, 1), Int3(1, 1, 1),
                      Int3(2, 2, 2), Int3(2, 6, 5),
                      gridSize - Int3(1, 1, 1), gridSize};
    for (int idx = 0; idx < sizeof(coords) / sizeof(*coords); idx++) {
        Particle particle = randomParticle();
        particle.setPosition(coords[idx]);
        ASSERT_EQ_INT3(cellIdx[idx], ensemble.getCellIndex(particle));
    }
}

TEST_F(EnsembleTest, DumpParticles) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    vector<Particle> particles;
    for (int i = 0; i < 73; i++) {
        particles.push_back(randomParticle());
        ensemble.add(particles[i]);
    }
    Int3 minCellIdxs[] = {Int3(0, 0, 0), Int3(1, 1, 1), parameters.localGridSize / 4};
    Int3 maxCellIdxs[] = {ensemble.numSystems, ensemble.numSystems - Int3(1, 1, 1),
        parameters.localGridSize / 2 + Int3(1, 1, 1)};
    for (int idx = 0; idx < sizeof(minCellIdxs) / sizeof(*minCellIdxs); idx++) {
        Int3 minCellIdx = minCellIdxs[idx];
        Int3 maxCellIdx = maxCellIdxs[idx];
        vector<Particle> expectedParticles;
        for (int i = 0; i < particles.size(); i++) {
            Int3 cellIdx = ensemble.getCellIndex(particles[i]);
            if ((cellIdx >= minCellIdx) && (cellIdx < maxCellIdx))
                expectedParticles.push_back(particles[i]);
        }
        vector<Particle> dumpedParticles(expectedParticles.size());
        ensemble.dumpParticles(ptr(dumpedParticles), &minCellIdx, &maxCellIdx);
        ASSERT_TRUE(eqParticleVectorArray(expectedParticles, ptr(dumpedParticles)));
    }
}

TEST_F(EnsembleTest, DumpParticleInCellCount) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    vector<Particle> particles;
    for (int i = 0; i < 81; i++) {
        particles.push_back(randomParticle());
        ensemble.add(particles[i]);
    }
    Int3 minCellIdxs[] = {Int3(0, 0, 0), Int3(1, 1, 1), parameters.localGridSize / 3};
    Int3 maxCellIdxs[] = {ensemble.numSystems, ensemble.numSystems - Int3(1, 1, 1),
        (parameters.localGridSize * 2) / 3 + Int3(1, 1, 1)};
    for (int idx = 0; idx < sizeof(minCellIdxs) / sizeof(*minCellIdxs); idx++) {
        Int3 minCellIdx = minCellIdxs[idx];
        Int3 maxCellIdx = maxCellIdxs[idx];
        Int3 size = maxCellIdx - minCellIdx;
        vector<int> expectedCounts(size.volume(), 0);
        for (int i = 0; i < particles.size(); i++) {
            Int3 cellIdx = ensemble.getCellIndex(particles[i]);
            if ((cellIdx >= minCellIdx) && (cellIdx < maxCellIdx)) {
                Int3 diff = cellIdx - minCellIdx;
                expectedCounts[diff.z + size.z * (diff.y + size.y * diff.x)]++;
            }
        }
        vector<int> dumpedCounts(expectedCounts.size(), -1);
        ensemble.dumpParticleInCellCount(ptr(dumpedCounts), &minCellIdx, &maxCellIdx);
        ASSERT_EQ(expectedCounts, dumpedCounts);
    }
}

TEST_F(EnsembleTest, DumpOutgoingParticles) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    vector<Particle> insideParticles;
    for (int i = 0; i < 23; i++) {
        insideParticles.push_back(randomParticle());
        ensemble.add(insideParticles.back());
    }
    vector<Particle> expectedOutsideParticles;
    for (int i = 0; i < 16; i++) {
        Particle particle = randomParticle();
        FP3 coords = particle.getPosition();
        for (int d = 0; d < 3; d++) {
            int r = rand() % 3;
            if (r == 0)
                coords[d] = parameters.localMin[d] - 0.1 * parameters.step[d];
            if (r == 2)
                coords[d] = parameters.localMax[d] + 0.8 * parameters.step[d];
        }
        // make sure new particles is outside
        int d = rand() % 3;
        coords[d] = parameters.localMax[d] + 0.2 * parameters.step[d];
        particle.setPosition(coords);
        expectedOutsideParticles.push_back(particle);
        ensemble.add(particle);
    }

    vector<Particle> outsideParticles;
    ensemble.dumpOutgoingParticles(outsideParticles);
    ASSERT_TRUE(eqParticleVectorArray(expectedOutsideParticles, ptr(outsideParticles)));
    ASSERT_TRUE(eqEnsembleParticleVector(ensemble, insideParticles));
}

TEST_F(EnsembleTest, GetTypeIdx) {
    Ensemble ensemble(parameters.localGridSize, parameters.localMin, parameters.step, 3);
    // hard coded names for now
    string typeName0 = "Electron";
    string typeName1 = "Proton";
    // check that the setup is right
    ASSERT_EQ(typeName0, Particle::types[0].name());
    ASSERT_EQ(typeName1, Particle::types[1].name());
    // actual test
    EXPECT_EQ(0, ensemble.getTypeIdx("Electron"));
    EXPECT_EQ(0, ensemble.getTypeIdx("electron"));
    EXPECT_EQ(1, ensemble.getTypeIdx("Proton"));
    EXPECT_EQ(1, ensemble.getTypeIdx("pROtON"));
}

