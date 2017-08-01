#include "TestingUtility.h"

#include "pica/ParticleSystem.h"

#include <vector>
using std::vector;

using namespace pica;


class ParticleSystemTest : public BaseParticleFixture {
};

TEST_F(ParticleSystemTest, DefaultConstructor) {
    ParticleSystem system;
    ASSERT_EQ(0, system.size());
}

TEST_F(ParticleSystemTest, Constructor) {
    ParticleSystem system(61);
    ASSERT_EQ(0, system.size());
}

TEST_F(ParticleSystemTest, CopyConstructor) {
    ParticleSystem system;
    for (int i = 0; i < 17; i++)
        system.add(randomParticle());
    ParticleSystem copySystem(system);
    ASSERT_TRUE(eqParticleSystems(system, copySystem));
}

TEST_F(ParticleSystemTest, Assignment) {
    ParticleSystem system, copySystem, copySystem2;
    for (int i = 0; i < 9; i++)
        system.add(randomParticle());
    for (int i = 0; i < 11; i++)
        copySystem.add(randomParticle());
    copySystem = system;
    copySystem2 = system;
    ASSERT_TRUE(eqParticleSystems(system, copySystem));
    ASSERT_TRUE(eqParticleSystems(system, copySystem2));
}

TEST_F(ParticleSystemTest, Add) {
    ParticleSystem system;
    for (int idx = 0; idx < 12; idx++) {
        Particle particle = randomParticle();
        system.add(particle);
        ASSERT_EQ(idx + 1, system.size());
        ASSERT_TRUE(eqParticles(particle, *(system.currentParticle + idx)));
    }
}

TEST_F(ParticleSystemTest, Raw) {
    ParticleSystem system;
    vector<Particle> particles;
    for (int idx = 0; idx < 23; idx++) {
        particles.push_back(randomParticle());
        system.add(particles.back());
    }
    ASSERT_EQ(particles.size(), system.size());
    Particle* raw = system.raw();
    const ParticleSystem* constSystemPtr = &system;
    const Particle* constRaw = constSystemPtr->raw();
    for (int idx = 0; idx < system.size(); idx++) {
        ASSERT_TRUE(eqParticles(particles[idx], raw[idx]));
        ASSERT_TRUE(eqParticles(particles[idx], constRaw[idx]));
    }
}

TEST_F(ParticleSystemTest, Begin) {
    ParticleSystem system;
    Particle firstParticle = randomParticle();
    system.add(firstParticle);
    for (int idx = 0; idx < 8; idx++) {
        Particle particle = randomParticle();
        system.add(particle);
    }
    for (int i = 1; i < system.size(); i++) {
        system.begin();
        for (int j = 0; j < i; j++)
            system.next();
        system.begin();
        ASSERT_TRUE(eqParticles(firstParticle, *system.currentParticle));
    }
}

TEST_F(ParticleSystemTest, Next) {
    ParticleSystem system;
    vector<Particle> particles;
    for (int idx = 0; idx < 9; idx++) {
        particles.push_back(randomParticle());
        system.add(particles.back());
    }
    for (int i = 0; i < system.size(); i++) {
        system.begin();
        for (int j = 0; j < i; j++)
            system.next();
        ASSERT_TRUE(eqParticles(particles[i], *system.currentParticle));
    }
}

TEST_F(ParticleSystemTest, End) {
    ParticleSystem emptySystem;
    ASSERT_TRUE(emptySystem.end());
    ParticleSystem system;
    for (int idx = 0; idx < 13; idx++)
        system.add(randomParticle());
    system.begin();
    for (int i = 0; i < system.size(); i++) {
        ASSERT_FALSE(system.end());
        system.next();
    }
    ASSERT_TRUE(system.end());
}

TEST_F(ParticleSystemTest, DeleteParticle) {
    ParticleSystem system;
    for (int idx = 0; idx < 14; idx++)
        system.add(randomParticle());
    int originalSize = system.size();
    for (int i = 0; i < originalSize; i++) {
        // copy system into a vector
        vector<Particle> particles;
        for (system.begin(); !system.end(); system.next())
            particles.push_back(*system.currentParticle);
        // delete from the system and the vector
        int deleteIdx = rand() % system.size();
        particles.erase(particles.begin() + deleteIdx);
        system.deleteParticle(deleteIdx);
        // check the system removed correct particle
        ASSERT_EQ(originalSize - i - 1, system.size());
        for (int j = 0; j < particles.size(); j++) {
            bool found = false;
            for (system.begin(); !system.end(); system.next())
                if (eqParticles(particles[j], *system.currentParticle)) {
                    found = true;
                    break;
                }
            ASSERT_TRUE(found);
        }
    }
}

TEST_F(ParticleSystemTest, DeleteCurrentParticle) {
    ParticleSystem system;
    for (int idx = 0; idx < 10; idx++)
        system.add(randomParticle());
    int originalSize = system.size();
    for (int i = 0; i < originalSize; i++) {
        // copy system into a vector
        vector<Particle> particles;
        for (system.begin(); !system.end(); system.next())
            particles.push_back(*system.currentParticle);
        // delete from the system and the vector
        int deleteIdx = rand() % system.size();
        particles.erase(particles.begin() + deleteIdx);
        system.begin();
        for (int j = 0; j < deleteIdx; j++)
            system.next();
        system.deleteCurrentParticle();
        // check the system removed correct particle
        ASSERT_EQ(originalSize - i - 1, system.size());
        for (int j = 0; j < particles.size(); j++) {
            bool found = false;
            for (system.begin(); !system.end(); system.next())
                if (eqParticles(particles[j], *system.currentParticle)) {
                    found = true;
                    break;
                }
            ASSERT_TRUE(found);
        }
    }
}

TEST_F(ParticleSystemTest, DeleteCurrentParticleTraversal) {
    ParticleSystem system;
    for (int idx = 0; idx < 10; idx++)
        system.add(randomParticle());
    int originalSize = system.size();
    int i = 0;
    for (system.begin(); !system.end(); system.next()) {
        // copy system into a vector
        vector<Particle> particles;
        for (int j = 0; j < system.size(); j++)
            particles.push_back(system.raw()[j]);
        particles.erase(particles.begin());
        system.deleteCurrentParticle();
        ASSERT_EQ(originalSize - i - 1, system.size());
        for (int j = 0; j < particles.size(); j++) {
            bool found = false;
            for (int k = 0; k < system.size(); k++)
                if (eqParticles(particles[j], system.raw()[k])) {
                    found = true;
                    break;
                }
            ASSERT_TRUE(found);
        }
        i++;
    }
}

TEST_F(ParticleSystemTest, DeleteCurrentParticleLast) {
    ParticleSystem system;
    for (int idx = 0; idx < 10; idx++)
        system.add(randomParticle());
    int originalSize = system.size();
    for (int i = 0; i < originalSize; i++) {
        // copy system into a vector
        vector<Particle> particles;
        for (int j = 0; j < system.size(); j++)
            particles.push_back(system.raw()[j]);
        particles.erase(particles.begin() + particles.size() - 1);
        system.begin();
        for (int j = 0; j < system.size() - 1; j++)
            system.next();
        system.deleteCurrentParticle();
        system.next();
        ASSERT_TRUE(system.end());
        ASSERT_EQ(originalSize - i - 1, system.size());
        for (int j = 0; j < particles.size(); j++) {
            bool found = false;
            for (int k = 0; k < system.size(); k++)
                if (eqParticles(particles[j], system.raw()[k])) {
                    found = true;
                    break;
                }
            ASSERT_TRUE(found);
        }
    }
}

TEST_F(ParticleSystemTest, Clear) {
    ParticleSystem emptySystem;
    emptySystem.clear();
    ASSERT_EQ(0, emptySystem.size());
    ParticleSystem system;
    for (int idx = 0; idx < 6; idx++)
        system.add(randomParticle());
    system.clear();
    ASSERT_EQ(0, system.size());
}

TEST_F(ParticleSystemTest, Size) {
    ParticleSystem emptySystem;
    ASSERT_EQ(0, emptySystem.size());
    ParticleSystem system;
    const int originalSize = 13;
    for (int idx = 0; idx < originalSize; idx++)
        system.add(randomParticle());
    ASSERT_EQ(originalSize, system.size());
    for (int idx = 0; idx < originalSize; idx++) {
        int deleteIdx = rand() % system.size();
        system.deleteParticle(deleteIdx);
        ASSERT_EQ(originalSize - idx - 1, system.size());
    }
}
