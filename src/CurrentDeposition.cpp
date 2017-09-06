#include "pica/CurrentDeposition.h"

#include "pica/fieldInterpolation/Formfactor.h"
#include "pica/grid/Grid.h"
#include "pica/math/Constants.h"
#include "pica/threading/OpenMPHelper.h"
#include "pica/Ensemble.h"

#include <algorithm>
#include <stdexcept>


namespace pica {

void CurrentDeposition::init(const Parameters& _parameters, Grid& _grid) {
    grid = &_grid;
    parameters = _parameters;
    currentDensityCoeff.resize(Particle::numTypes);
    esirkepovCurrentDensityCoeff.resize(Particle::numTypes);
    chargeCoeff.resize(Particle::numTypes);
    for (int i = 0; i < Particle::numTypes; ++i)
    {
        currentDensityCoeff[i] = constants::c * Particle::types[i].charge / grid->steps.volume();
        esirkepovCurrentDensityCoeff[i] = grid->steps * (Particle::types[i].charge /
            (grid->dt * grid->steps.volume()));
        chargeCoeff[i] = Particle::types[i].charge / grid->steps.volume() / grid->dt;
    }

    ///
    std::string type = "CIC";
    for (int i = 0; i < type.length(); i++)
        type[i] = ::toupper(type[i]);
    if (type == "CIC")
        grid->setDepositionType(Grid::Deposition_CIC);
    else if ((type == "VILLASENORBUNEMAN") || (type == "VB"))
        grid->setDepositionType(Grid::Deposition_VillasenorBuneman);
    else if ((type == "ZIGZAGFIRSTORDER") || (type == "ZZ1"))
        grid->setDepositionType(Grid::Deposition_ZigzagFirstOrder);
    else if (type == "TSC")
        grid->setDepositionType(Grid::Deposition_TSC);
    else if ((type == "ZIGZAGSECONDORDER") || (type == "ZZ2"))
        grid->setDepositionType(Grid::Deposition_ZigzagSecondOrder);
    else if (type == "PCS")
        grid->setDepositionType(Grid::Deposition_PCS);
    else if (type == "ESIRKEPOV")
        grid->setDepositionType(Grid::Deposition_Esirkepov);
    ///else
        ///pic->computationLog.writeError("unsupported current deposition type: ");
}


/* Base class for tile current depositors. It is responsible for allocating
memory for tile currents and writing them back to grid current values. */
template<int TileSize>
class TileDepositor
{
public:

    static const int tileSize = TileSize;

    TileDepositor(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
    {
        grid = currentDeposition->grid;
        currentDensityCoeff = ptr(currentDeposition->currentDensityCoeff);
        esirkepovCurrentDensityCoeff = ptr(currentDeposition->esirkepovCurrentDensityCoeff);
        dimensionality = currentDeposition->parameters.dimensionality;
        int tileOffset = (tileSize + 1) / 2;
        baseGridIdx = Int3(tileI, tileJ, tileK) - Int3(tileOffset, tileOffset, tileOffset)
            + grid->getNumExternalLeftCells();
        for (int d = dimensionality; d < 3; d++)
            baseGridIdx[d] = grid->getNumExternalLeftCells()[d];
        invSteps = FP3(1.0, 1.0, 1.0) / grid->steps;
        normalizedCellOrigin = grid->origin / grid->steps + FP3(baseGridIdx);
        normalizedMiddleCellOrigin = normalizedCellOrigin + FP3(0.5, 0.5, 0.5);
        for (int i = 0; i < tileSize; ++i)
        for (int j = 0; j < tileSize; ++j)
        for (int k = 0; k < tileSize; ++k)
        {
            jx[i][j][k] = 0;
            jy[i][j][k] = 0;
            jz[i][j][k] = 0;
        }
    }

    ~TileDepositor()
    {
        for (int i = 0; i < tileSize; ++i)
        for (int j = 0; j < tileSize; ++j)
        for (int k = 0; k < tileSize; ++k)
        {
            Int3 gridIdx = remainder(Int3(i, j, k) + baseGridIdx, grid->numCells);
            grid->Jx(gridIdx) += jx[i][j][k];
            grid->Jy(gridIdx) += jy[i][j][k];
            grid->Jz(gridIdx) += jz[i][j][k];
        }
    }

    // This is just to declare common interface for derived classes,
    // they are used as objects of the exact class anyway
    virtual void depositCurrent(const Particle& particle) = 0;

protected:
    
    FP jx[tileSize][tileSize][tileSize];
    FP jy[tileSize][tileSize][tileSize];
    FP jz[tileSize][tileSize][tileSize];
    Int3 baseGridIdx;
    Grid* grid;
    FP* currentDensityCoeff;
    FP3* esirkepovCurrentDensityCoeff;
    FP3 invSteps;
    FP3 normalizedCellOrigin, normalizedMiddleCellOrigin;
    int dimensionality;

};


class TileDepositorCIC: public TileDepositor<4>
{
public:

    TileDepositorCIC(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
        : TileDepositor<4>(tileI, tileJ, tileK, currentDeposition),
          halfDt((FP)0.5 * currentDeposition->parameters.timeStep) {}

    void depositCurrent(const Particle& particle)
    {
        FP3 current = particle.p * particle.invGamma
            * currentDensityCoeff[particle.getType()] * particle.getFactor();
        FP3 normalizedCoords =
            (particle.coords - particle.getVelocity() * halfDt) * invSteps;
        FP3 internalCoords = normalizedCoords - normalizedCellOrigin;
        Int3 idx = truncate(internalCoords);
        FP3 coeff = internalCoords - FP3(idx);
        internalCoords = normalizedCoords - normalizedMiddleCellOrigin;
        Int3 middleIdx = truncate(internalCoords);
        FP3 middleCoeff = internalCoords - FP3(middleIdx);
        depositComponent(current.x, coeff.x, middleCoeff.y, middleCoeff.z,
            idx.x, middleIdx.y, middleIdx.z, jx);
        depositComponent(current.y, middleCoeff.x, coeff.y, middleCoeff.z,
            middleIdx.x, idx.y, middleIdx.z, jy);
        depositComponent(current.z, middleCoeff.x, middleCoeff.y, coeff.z,
            middleIdx.x, middleIdx.y, idx.z, jz);
    }

private:

    // Deposit current component to base index (i, j, k) with coefficients
    // (x, y, z).
    void depositComponent(FP value, FP x, FP y, FP z, int i, int j,
        int k, FP component[tileSize][tileSize][tileSize])
    {
        const FP xv = x * value;
        const FP ixv = value - xv;
        const FP iy = (FP)1 - y;
        const FP iz = (FP)1 - z;
        const FP iyiz = iy * iz;
        const FP iyz = iy * z;
        const FP yiz = y * iz;
        const FP yz = y * z;
        component[i][j][k] += ixv * iyiz;
        component[i][j][k + 1] += ixv * iyz;
        component[i][j + 1][k] += ixv * yiz;
        component[i][j + 1][k + 1] += ixv * yz;
        component[i + 1][j][k] += xv * iyiz;
        component[i + 1][j][k + 1] += xv * iyz;
        component[i + 1][j + 1][k] += xv * yiz;
        component[i + 1][j + 1][k + 1] += xv * yz;
    }

    const FP halfDt;

};


class TileDepositorVillasenorBuneman: public TileDepositor<3>
{
public:

    TileDepositorVillasenorBuneman(int tileI, int tileJ, int tileK,
        CurrentDeposition* currentDeposition)
        : TileDepositor<3>(tileI, tileJ, tileK, currentDeposition)
    {
        chargeCoeff = ptr(currentDeposition->chargeCoeff);
        steps = grid->steps;
        minCoords = grid->origin;
        middleCellOrigin = grid->origin + 0.5 * grid->steps;
        dt = grid->dt;
    }

    void depositCurrent(const Particle& particle)
    {
        FP charge = chargeCoeff[particle.getType()] * particle.getFactor();
        FP3 coords = particle.coords;
        FP3 velocity = particle.getVelocity();
        FP3 oldCoords = particle.coords - velocity * dt;
        Int3 oldLocalOrigin = truncate((oldCoords - middleCellOrigin) * invSteps);
        Int3 localOrigin = truncate((coords - middleCellOrigin) * invSteps);
        for (int d = dimensionality; d < 3; d++)
            oldLocalOrigin[d] = localOrigin[d] = 0;
        depositCurrentRecursive(coords, oldCoords, localOrigin, oldLocalOrigin, velocity, charge);
    }

private:

    void depositCurrentRecursive(FP3 coords, FP3 oldCoords, Int3 localOrigin,
                                 Int3 oldLocalOrigin, FP3 velocity, FP charge)
    {
        /* Find the axis that is first intersected.
        Take max time as it is measured for backward movement. */
        FP maxTimeToIntersection = 0;
        int earliestIntersectionAxis = -1;
        for (int d = 0; d < dimensionality; d++)
        {
            if (localOrigin[d] != oldLocalOrigin[d])
            {
                FP newOldCoords = minCoords[d] + 0.5 * steps[d] +
                    std::max(oldLocalOrigin[d], localOrigin[d]) * steps[d];
                FP timeToIntersection = (coords[d] - newOldCoords) / velocity[d];
                if (timeToIntersection > maxTimeToIntersection)
                {
                    maxTimeToIntersection = timeToIntersection;
                    earliestIntersectionAxis = d;
                }
            }
        }
        if (earliestIntersectionAxis >= 0)
        {
            /* Split particle along axis dim0; dim1 and dim2 are two other axes
            (0 = x, 1 = y, 2 = z) */
            const int dim0 = earliestIntersectionAxis;
            const int dim1 = (dim0 + 1) % 3;
            const int dim2 = (dim0 + 2) % 3;
            FP3 newOldCoords;
            newOldCoords[dim0] = middleCellOrigin[dim0] +
                std::max(oldLocalOrigin[dim0], localOrigin[dim0]) * steps[dim0];
            newOldCoords[dim1] = coords[dim1] - velocity[dim1] * maxTimeToIntersection;
            newOldCoords[dim2] = coords[dim2] - velocity[dim2] * maxTimeToIntersection;
            Int3 newOldLocalOrigin = oldLocalOrigin;
            newOldLocalOrigin[dim0] = localOrigin[dim0];
            depositCurrentRecursive(coords, newOldCoords, localOrigin, newOldLocalOrigin, velocity, charge);
            coords = newOldCoords;
            localOrigin = oldLocalOrigin;
        }

        FP3 delta = (coords - oldCoords) * invSteps;
        FP3 midway = (0.5 * (coords + oldCoords) - minCoords) * invSteps -
            FP3(localOrigin) - FP3(0.5, 0.5, 0.5);

        FP3 iMidway = FP3(1.0, 1.0, 1.0) - midway;
        FP3 csteps = charge * steps;
        FP3 delcs = delta * csteps;
        FP3 csdelmid, csdelimid;
        csdelmid[0] = delcs.x * midway.z;
        csdelmid[1] = delcs.y * midway.z;
        csdelmid[2] = delcs.z * midway.x;
        csdelimid[0] = delcs.x * iMidway.z;
        csdelimid[1] = delcs.y * iMidway.z;
        csdelimid[2] = delcs.z * iMidway.x;
        FP3 csdel = csteps * delta.x * delta.y * delta.z / 12;
        Int3 idx = localOrigin - baseGridIdx;
        Int3 iIdx = idx + Int3(1, 1, 1);

        jx[iIdx.x][iIdx.y][iIdx.z] += csdelmid[0] * midway.y + csdel.x;
        jx[iIdx.x][idx.y][iIdx.z] += csdelmid[0] * iMidway.y - csdel.x;
        jx[iIdx.x][iIdx.y][idx.z] += csdelimid[0] * midway.y - csdel.x;
        jx[iIdx.x][idx.y][idx.z] += csdelimid[0] * iMidway.y + csdel.x;

        jy[iIdx.x][iIdx.y][iIdx.z] += csdelmid[1] * midway.x + csdel.y;
        jy[idx.x][iIdx.y][iIdx.z] += csdelmid[1] * iMidway.x - csdel.y;
        jy[iIdx.x][iIdx.y][idx.z] += csdelimid[1] * midway.x - csdel.y;
        jy[idx.x][iIdx.y][idx.z] += csdelimid[1] * iMidway.x + csdel.y;

        jz[iIdx.x][iIdx.y][iIdx.z] += csdelmid[2] * midway.y + csdel.z;
        jz[iIdx.x][idx.y][iIdx.z] += csdelmid[2] * iMidway.y - csdel.z;
        jz[idx.x][iIdx.y][iIdx.z] += csdelimid[2] * midway.y - csdel.z;
        jz[idx.x][idx.y][iIdx.z] += csdelimid[2] * iMidway.y + csdel.z;
    }

    FP* chargeCoeff;
    FP3 steps, minCoords, middleCellOrigin;
    FP dt;
};


class TileDepositorZigzagFirstOrder: public TileDepositor<3>
{
public:

    TileDepositorZigzagFirstOrder(int tileI, int tileJ, int tileK,
        CurrentDeposition* currentDeposition)
        : TileDepositor<3>(tileI, tileJ, tileK, currentDeposition) {}

    void depositCurrent(const Particle& particle)
    {
        FP dt = grid->dt;
        FP3 steps = grid->steps;
        FP3 coords = particle.coords;
        FP3 minCoords = grid->origin;
        FP3 velocity = particle.getVelocity();

        FP3 oldCoords = coords - velocity * dt;
        FP charge = particle.charge() * particle.getFactor() / steps.volume() / dt;
        Int3 localOrigin = truncate((coords - minCoords) * invSteps + FP3(0.5, 0.5, 0.5));
        Int3 oldLocalOrigin = truncate((oldCoords - minCoords) * invSteps + FP3(0.5, 0.5, 0.5));

        FP3 r;
        for(int i = 0; i < 3; i++)
        {
            if(oldLocalOrigin[i] == localOrigin[i]) r[i] = (coords[i] + oldCoords[i]) * (FP)0.5;
            else r[i] = (oldLocalOrigin[i] + localOrigin[i]) * steps[i] * (FP)0.5 + minCoords[i];
        }

        FP3 F[2], W1[2], W2[2];
        W1[0] = ((oldCoords + r) * (FP)0.5 - minCoords) * invSteps - (FP3)oldLocalOrigin + FP3(0.5, 0.5, 0.5);
        W1[1] = ((coords + r) * (FP)0.5 - minCoords) * invSteps - (FP3)localOrigin + FP3(0.5, 0.5, 0.5);
        W2[0] = FP3(1.0, 1.0, 1.0) - W1[0];
        W2[1] = FP3(1.0, 1.0, 1.0) - W1[1];
        F[0] = charge * (r - oldCoords);
        F[1] = charge * (coords - r);

        Int3 idx[2];
        idx[0] = oldLocalOrigin - baseGridIdx;
        idx[1] = localOrigin - baseGridIdx;
        for (int d = dimensionality; d < 3; d++)
            idx[0][d] = idx[1][d] = 1;
        for(int i = 0; i < 2; i++)
        {
            jx[idx[i].x][idx[i].y]    [idx[i].z]     += F[i].x * W1[i].y * W1[i].z;
            jx[idx[i].x][idx[i].y - 1][idx[i].z]     += F[i].x * W2[i].y * W1[i].z;
            jx[idx[i].x][idx[i].y]    [idx[i].z - 1] += F[i].x * W1[i].y * W2[i].z;
            jx[idx[i].x][idx[i].y - 1][idx[i].z - 1] += F[i].x * W2[i].y * W2[i].z;

            jy[idx[i].x]    [idx[i].y][idx[i].z]     += F[i].y * W1[i].x * W1[i].z;
            jy[idx[i].x - 1][idx[i].y][idx[i].z]     += F[i].y * W2[i].x * W1[i].z;
            jy[idx[i].x]    [idx[i].y][idx[i].z - 1] += F[i].y * W1[i].x * W2[i].z;
            jy[idx[i].x - 1][idx[i].y][idx[i].z - 1] += F[i].y * W2[i].x * W2[i].z;

            jz[idx[i].x]    [idx[i].y]    [idx[i].z] += F[i].z * W1[i].x * W1[i].y;
            jz[idx[i].x - 1][idx[i].y]    [idx[i].z] += F[i].z * W2[i].x * W1[i].y;
            jz[idx[i].x]    [idx[i].y - 1][idx[i].z] += F[i].z * W1[i].x * W2[i].y;
            jz[idx[i].x - 1][idx[i].y - 1][idx[i].z] += F[i].z * W2[i].x * W2[i].y;
        }
    }
};


class TileDepositorZigzagSecondOrder: public TileDepositor<4>
{
public:

    TileDepositorZigzagSecondOrder(int tileI, int tileJ, int tileK,
        CurrentDeposition* currentDeposition)
        : TileDepositor<4>(tileI, tileJ, tileK, currentDeposition) {}

    void depositCurrent(const Particle& particle)
    {
        FP dt = grid->dt;
        FP3 steps = grid->steps;
        FP3 coords = particle.coords;
        FP3 minCoords = grid->origin;
        FP3 velocity = particle.getVelocity();

        FP3 oldCoords = coords - velocity * dt;
        FP charge = particle.charge() * particle.getFactor() / steps.volume();
        Int3 localOrigin = truncate((coords - minCoords) * invSteps);
        Int3 oldLocalOrigin = truncate((oldCoords - minCoords) * invSteps);

        FP3 r;
        for (int i = 0; i < 3; i++)
        {
            if (oldLocalOrigin[i] == localOrigin[i]) r[i] = (coords[i] + oldCoords[i]) * (FP)0.5;
            else r[i] = std::max(oldLocalOrigin[i], localOrigin[i]) * steps[i] + minCoords[i];
        }

        FP3 W[2];
        W[0] = ((oldCoords + r) * (FP)0.5 - minCoords) * invSteps - (FP3)oldLocalOrigin -
            FP3(0.5, 0.5, 0.5);
        W[1] = ((coords + r) * (FP)0.5 - minCoords) * invSteps - (FP3)localOrigin -
            FP3(0.5, 0.5, 0.5);

        FP3 F[2][2];
        FP3 cv = charge * velocity * (FP)0.5;
        F[0][0] = cv * (FP3(0.5, 0.5, 0.5) - W[0]);
        F[1][0] = cv * (FP3(0.5, 0.5, 0.5) - W[1]);
        F[0][1] = cv * (FP3(0.5, 0.5, 0.5) + W[0]);
        F[1][1] = cv * (FP3(0.5, 0.5, 0.5) + W[1]);

        FP3 W1[2], W2[2], W3[2];
        W1[0] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) - W[0]) * (FP3(0.5, 0.5, 0.5) - W[0]);
        W1[1] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) - W[1]) * (FP3(0.5, 0.5, 0.5) - W[1]);

        W2[0] = FP3(0.75, 0.75, 0.75) - W[0] * W[0];
        W2[1] = FP3(0.75, 0.75, 0.75) - W[1] * W[1];

        W3[0] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) + W[0]) * (FP3(0.5, 0.5, 0.5) + W[0]);
        W3[1] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) + W[1]) * (FP3(0.5, 0.5, 0.5) + W[1]);

        Int3 idx[2];
        idx[0] = oldLocalOrigin - baseGridIdx;
        idx[1] = localOrigin - baseGridIdx;
        for (int d = dimensionality; d < 3; d++)
            idx[0][d] = idx[1][d] = 1;

        for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
        {
            jx[idx[i].x + j][idx[i].y - 1][idx[i].z - 1] += F[i][j].x * W1[i].y * W1[i].z;
            jx[idx[i].x + j][idx[i].y - 1][idx[i].z]     += F[i][j].x * W1[i].y * W2[i].z;
            jx[idx[i].x + j][idx[i].y - 1][idx[i].z + 1] += F[i][j].x * W1[i].y * W3[i].z;
            jx[idx[i].x + j][idx[i].y]    [idx[i].z - 1] += F[i][j].x * W2[i].y * W1[i].z;
            jx[idx[i].x + j][idx[i].y]    [idx[i].z]     += F[i][j].x * W2[i].y * W2[i].z;
            jx[idx[i].x + j][idx[i].y]    [idx[i].z + 1] += F[i][j].x * W2[i].y * W3[i].z;
            jx[idx[i].x + j][idx[i].y + 1][idx[i].z - 1] += F[i][j].x * W3[i].y * W1[i].z;
            jx[idx[i].x + j][idx[i].y + 1][idx[i].z]     += F[i][j].x * W3[i].y * W2[i].z;
            jx[idx[i].x + j][idx[i].y + 1][idx[i].z + 1] += F[i][j].x * W3[i].y * W3[i].z;
            
            jy[idx[i].x - 1][idx[i].y + j][idx[i].z - 1] += F[i][j].y * W1[i].x * W1[i].z;
            jy[idx[i].x - 1][idx[i].y + j][idx[i].z]     += F[i][j].y * W1[i].x * W2[i].z;
            jy[idx[i].x - 1][idx[i].y + j][idx[i].z + 1] += F[i][j].y * W1[i].x * W3[i].z;
            jy[idx[i].x]    [idx[i].y + j][idx[i].z - 1] += F[i][j].y * W2[i].x * W1[i].z;
            jy[idx[i].x]    [idx[i].y + j][idx[i].z]     += F[i][j].y * W2[i].x * W2[i].z;
            jy[idx[i].x]    [idx[i].y + j][idx[i].z + 1] += F[i][j].y * W2[i].x * W3[i].z;
            jy[idx[i].x + 1][idx[i].y + j][idx[i].z - 1] += F[i][j].y * W3[i].x * W1[i].z;
            jy[idx[i].x + 1][idx[i].y + j][idx[i].z]     += F[i][j].y * W3[i].x * W2[i].z;
            jy[idx[i].x + 1][idx[i].y + j][idx[i].z + 1] += F[i][j].y * W3[i].x * W3[i].z;

            jz[idx[i].x - 1][idx[i].y - 1][idx[i].z + j] += F[i][j].z * W1[i].y * W1[i].x;
            jz[idx[i].x]    [idx[i].y - 1][idx[i].z + j] += F[i][j].z * W1[i].y * W2[i].x;
            jz[idx[i].x + 1][idx[i].y - 1][idx[i].z + j] += F[i][j].z * W1[i].y * W3[i].x;
            jz[idx[i].x - 1][idx[i].y]    [idx[i].z + j] += F[i][j].z * W2[i].y * W1[i].x;
            jz[idx[i].x]    [idx[i].y]    [idx[i].z + j] += F[i][j].z * W2[i].y * W2[i].x;
            jz[idx[i].x + 1][idx[i].y]    [idx[i].z + j] += F[i][j].z * W2[i].y * W3[i].x;
            jz[idx[i].x - 1][idx[i].y + 1][idx[i].z + j] += F[i][j].z * W3[i].y * W1[i].x;
            jz[idx[i].x]    [idx[i].y + 1][idx[i].z + j] += F[i][j].z * W3[i].y * W2[i].x;
            jz[idx[i].x + 1][idx[i].y + 1][idx[i].z + j] += F[i][j].z * W3[i].y * W3[i].x;
        }
    }
};


class TileDepositorTSC: public TileDepositor<5>
{
public:

    TileDepositorTSC(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
        : TileDepositor<5>(tileI, tileJ, tileK, currentDeposition),
          halfDt((FP)0.5 * currentDeposition->parameters.timeStep) {}

    void depositCurrent(const Particle& particle)
    {
        FP3 current = particle.p * particle.invGamma
            * currentDensityCoeff[particle.getType()] * particle.getFactor();
        FP3 normalizedCoords =
            (particle.coords - particle.getVelocity() * halfDt) * invSteps;
        FP3 internalCoords = normalizedCoords - normalizedCellOrigin;
        Int3 idx = truncate(internalCoords + FP3(0.5, 0.5, 0.5));
        FP3 coeff = internalCoords - FP3(idx);
        internalCoords = normalizedCoords - normalizedMiddleCellOrigin;
        Int3 middleIdx = truncate(internalCoords + FP3(0.5, 0.5, 0.5));
        FP3 middleCoeff = internalCoords - FP3(middleIdx);
        depositComponent(current.x, coeff.x, middleCoeff.y, middleCoeff.z,
            idx.x, middleIdx.y, middleIdx.z, jx);
        depositComponent(current.y, middleCoeff.x, coeff.y, middleCoeff.z,
            middleIdx.x, idx.y, middleIdx.z, jy);
        depositComponent(current.z, middleCoeff.x, middleCoeff.y, coeff.z,
            middleIdx.x, middleIdx.y, idx.z, jz);
    }

private:

    void depositComponent(FP value, FP x, FP y, FP z, int i, int j,
        int k, FP component[tileSize][tileSize][tileSize])
    {
        FP cx[3], cy[3], cz[3];
        for (int ii = 0; ii < 3; ii++)
            cx[ii] = formfactorTSC(FP(ii - 1) - x);
        for (int jj = 0; jj < 3; jj++)
            cy[jj] = formfactorTSC(FP(jj - 1) - y);
        for (int kk = 0; kk < 3; kk++)
            cz[kk] = formfactorTSC(FP(kk - 1) - z);
        for (int ii = -1; ii <= 1; ii++)
        for (int jj = -1; jj <= 1; jj++)
        for (int kk = -1; kk <= 1; kk++)
            component[i + ii][j + jj][k + kk] += cx[ii + 1] * cy[jj + 1] *
                cz[kk + 1] * value;
    }

    const FP halfDt;

};



class TileDepositorEsirkepov : public TileDepositor<5>
{
public:

    TileDepositorEsirkepov(int tileI, int tileJ, int tileK,
        CurrentDeposition* currentDeposition)
        : TileDepositor<5>(tileI, tileJ, tileK, currentDeposition)
    {
        centralNodeIdx = Int3(tileI + 1, tileJ + 1, tileK + 1);
        centralNodePosition = grid->RhoPosition(centralNodeIdx.x, centralNodeIdx.y, centralNodeIdx.z);
    }

    void depositCurrent(const Particle& particle)
    {
        FP3 position = particle.getPosition();
        FP3 internalPosition = (position - centralNodePosition) * invSteps;
        FP3 S1[5];
        formFactorEsirkepov(internalPosition, S1, Int3(0, 0, 0));
        FP3 oldPosition = position - particle.getVelocity() * grid->dt;
        FP3 oldInternalPosition = (oldPosition - centralNodePosition) * invSteps;
        Int3 shifts;
        for (int d = 0; d < 3; ++d)
            if (oldInternalPosition[d] < (FP)-0.5)
            {
                oldInternalPosition[d] += (FP)1.0;
                --shifts[d];
            }
            else
                if (oldInternalPosition[d] > (FP)0.5)
                {
                    oldInternalPosition[d] -= (FP)1.0;
                    ++shifts[d];
                }
        FP3 S0[5];
        formFactorEsirkepov(oldInternalPosition, S0, shifts);
        FP3 DS[5];
        for (int i = 0; i < 5; i++)
            DS[i] = S1[i] - S0[i];

        const FP oneThird = 1.0 / 3.0;
        FP3 normalizedQ = esirkepovCurrentDensityCoeff[particle.getType()] * particle.getFactor();
        for (int i = 0; i < 5; i++)
        {
            FP DSix = DS[i].x * oneThird;
            for (int j = 0; j < 5; j++)
            {
                FP DSixDSjy = DSix * DS[j].y;
#pragma novector
                for (int k = 0; k < 5; k++)
                {
                    FP DSixDSjyDSkz = DSixDSjy * DS[k].z;
                    jx[i][j][k] += normalizedQ.x * ((FP)0.5 * DS[i].x *
                        (S0[j].y * S1[k].z + S1[j].y * S0[k].z) +
                        DSixDSjyDSkz);
                    jy[i][j][k] += normalizedQ.y * ((FP)0.5 * DS[j].y *
                        (S0[i].x * S1[k].z + S1[i].x * S0[k].z) +
                        DSixDSjyDSkz);
                    jz[i][j][k] += normalizedQ.z * ((FP)0.5 * DS[k].z *
                        (S0[i].x * S1[j].y + S1[i].x * S0[j].y) +
                        DSixDSjyDSkz);
                }
            }
        }
    }

    ~TileDepositorEsirkepov()
    {
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 5; k++)
            {
                FP element = 0;
                for (int i = 4; i >= 0; i--)
                {
                    element += jx[i][j][k];
                    jx[i][j][k] = element;
                }
            }
        for (int i = 0; i < 5; i++)
            for (int k = 0; k < 5; k++)
            {
                FP element = 0;
                for (int j = 4; j >= 0; j--)
                {
                    element += jy[i][j][k];
                    jy[i][j][k] = element;
                }
            }
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
            {
                FP element = 0;
                for (int k = 4; k >= 0; k--)
                {
                    element += jz[i][j][k];
                    jz[i][j][k] = element;
                }
            }
    }

private:

    Int3 centralNodeIdx;
    FP3 centralNodePosition;

};


class TileDepositorPCS: public TileDepositor<6>
{
public:

    TileDepositorPCS(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
        : TileDepositor<6>(tileI, tileJ, tileK, currentDeposition) {}

    void depositCurrent(const Particle& particle)
    {
        FP3 current = particle.p * particle.invGamma
            * currentDensityCoeff[particle.getType()] * particle.getFactor();
        grid->depositJPCS(current, particle.coords);
    }

};


/* Base class for tile current depositors in XY 2D case. It is responsible for allocating
memory for tile currents and writing them back to grid current values. */
template<int TileSize>
class XYTileDepositor
{
public:
    
    static const int tileSize = TileSize;

    XYTileDepositor(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
    {
        grid = currentDeposition->grid;
        currentDensityCoeff = ptr(currentDeposition->currentDensityCoeff);
        esirkepovCurrentDensityCoeff = ptr(currentDeposition->esirkepovCurrentDensityCoeff);
        int tileOffset = (tileSize + 1) / 2;
        baseGridIdx = Int3(tileI, tileJ, 0) - Int3(tileOffset, tileOffset, 0)
            + grid->getNumExternalLeftCells();
        invSteps = FP3(1.0, 1.0, 1.0) / grid->steps;
        normalizedCellOrigin = grid->origin / grid->steps + FP3(baseGridIdx);
        normalizedMiddleCellOrigin = normalizedCellOrigin + FP3(0.5, 0.5, 0.5);
        for (int i = 0; i < tileSize; ++i)
        for (int j = 0; j < tileSize; ++j)
        {
            jx[i][j] = 0;
            jy[i][j] = 0;
            jz[i][j] = 0;
        }
    }

    ~XYTileDepositor()
    {
        for (int i = 0; i < tileSize; ++i)
        for (int j = 0; j < tileSize; ++j)
        {
            Int3 gridIdx = Int3(i, j, 0) + baseGridIdx;
            grid->Jx(gridIdx) += jx[i][j];
            grid->Jy(gridIdx) += jy[i][j];
            grid->Jz(gridIdx) += jz[i][j];
        }
    }

    // Same as in TileDepositor
    virtual void depositCurrent(const Particle& particle) = 0;

protected:
    
    FP jx[tileSize][tileSize];
    FP jy[tileSize][tileSize];
    FP jz[tileSize][tileSize];
    Int3 baseGridIdx;
    Grid* grid;
    FP* currentDensityCoeff;
    FP3* esirkepovCurrentDensityCoeff;
    FP3 invSteps;
    FP3 normalizedCellOrigin, normalizedMiddleCellOrigin;

};


class XYTileDepositorTSC: public XYTileDepositor<5>
{
public:

    XYTileDepositorTSC(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
        : XYTileDepositor<5>(tileI, tileJ, tileK, currentDeposition),
          halfDt((FP)0.5 * currentDeposition->parameters.timeStep) {}

    void depositCurrent(const Particle& particle)
    {
        FP3 current = particle.p * particle.invGamma
            * currentDensityCoeff[particle.getType()] * particle.getFactor();
        FP3 halfV = particle.getVelocity() * halfDt;
        FP3 normalizedCoords((particle.coords.x - halfV.x) * invSteps.x,
                             (particle.coords.y - halfV.y) * invSteps.y, 0);
        FP3 internalCoords(normalizedCoords.x - normalizedCellOrigin.x,
                           normalizedCoords.y - normalizedCellOrigin.y, 0);
        Int3 idx((int)(internalCoords.x + (FP)0.5),
                 (int)(internalCoords.y + (FP)0.5), 0);
        FP3 coeff(internalCoords.x - (FP)idx.x, internalCoords.y - (FP)idx.y, 0);
        internalCoords.x = normalizedCoords.x - normalizedMiddleCellOrigin.x;
        internalCoords.y = normalizedCoords.y - normalizedMiddleCellOrigin.y;
        Int3 middleIdx((int)(internalCoords.x + (FP)0.5),
                       (int)(internalCoords.y + (FP)0.5), 0);
        FP3 middleCoeff(internalCoords.x - (FP)middleIdx.x,
                        internalCoords.y - (FP)middleIdx.y, 0);
        depositComponent(current.x, coeff.x, middleCoeff.y, idx.x, middleIdx.y,
            jx);
        depositComponent(current.y, middleCoeff.x, coeff.y, middleIdx.x, idx.y,
            jy);
        depositComponent(current.z, middleCoeff.x, middleCoeff.y, middleIdx.x,
            middleIdx.y, jz);
    }

private:

    void depositComponent(FP value, FP x, FP y, int i, int j,
        FP component[tileSize][tileSize])
    {
        FP cx[3], cy[3];
        for (int ii = 0; ii < 3; ii++)
            cx[ii] = formfactorTSC(FP(ii - 1) - x);
        for (int jj = 0; jj < 3; jj++)
            cy[jj] = formfactorTSC(FP(jj - 1) - y);
        for (int ii = -1; ii <= 1; ii++)
        for (int jj = -1; jj <= 1; jj++)
            component[i + ii][j + jj] += cx[ii + 1] * cy[jj + 1] * value;
    }

    const FP halfDt;

};


class XYTileDepositorEsirkepov : public XYTileDepositor<5>
{
public:
    XYTileDepositorEsirkepov(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
        : XYTileDepositor<5>(tileI, tileJ, tileK, currentDeposition)
    {
        centralNodeIdx = Int3(tileI + 1, tileJ + 1, 0);
        centralNodePosition = grid->RhoPosition(centralNodeIdx.x, centralNodeIdx.y, centralNodeIdx.z);
    }

    void depositCurrent(const Particle& particle)
    {
        FP3 position = particle.getPosition();
        FP3 internalPosition = (position - centralNodePosition) * invSteps;
        FP3 S1[5];
        formFactorEsirkepovXY(internalPosition, S1, Int3(0, 0, 0));
        FP3 oldPosition = position - particle.getVelocity() * grid->dt;
        FP3 oldInternalPosition = (oldPosition - centralNodePosition) * invSteps;
        Int3 shifts;
        for (int d = 0; d < 2; ++d)
            if (oldInternalPosition[d] < (FP)-0.5)
            {
                oldInternalPosition[d] += (FP)1.0;
                --shifts[d];
            }
            else
                if (oldInternalPosition[d] > (FP)0.5)
                {
                    oldInternalPosition[d] -= (FP)1.0;
                    ++shifts[d];
                }
        FP3 S0[5];
        formFactorEsirkepovXY(oldInternalPosition, S0, shifts);
        FP3 DS[5];
        for (int i = 0; i < 5; i++)
        {
            DS[i].x = S1[i].x - S0[i].x;
            DS[i].y = S1[i].y - S0[i].y;
        }

        // The following computations are done according to (35-36) in the paper
        const FP oneThird = 1.0 / 3.0;
        FP3 normalizedQ = esirkepovCurrentDensityCoeff[particle.getType()] * particle.getFactor();
        normalizedQ.z = particle.getVelocity().z * particle.charge()
            * particle.getFactor() / grid->steps.volume();
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
            {
                jx[i][j] += normalizedQ.x * DS[i].x *
                    (S0[j].y + (FP)0.5 * DS[j].y);
                jy[i][j] += normalizedQ.y * DS[j].y *
                    (S0[i].x + (FP)0.5 * DS[i].x);
                jz[i][j] += normalizedQ.z * (S0[i].x * S0[j].y +
                    (FP)0.5 * (DS[i].x * S0[j].y + DS[j].y * S0[i].x) +
                    oneThird * DS[i].x * DS[j].y);
            }
    }

    ~XYTileDepositorEsirkepov()
    {
        for (int j = 0; j < 5; j++)
        {
            FP element = 0;
            for (int i = 4; i >= 0; i--)
            {
                element += jx[i][j];
                jx[i][j] = element;
            }
        }
        for (int i = 0; i < 5; i++)
        {
            FP element = 0;
            for (int j = 4; j >= 0; j--)
            {
                element += jy[i][j];
                jy[i][j] = element;
            }
        }
    }

private:

    Int3 centralNodeIdx;
    FP3 centralNodePosition;

};


/* Base class for tile current depositors in X 1D case. It is responsible for allocating
memory for tile currents and writing them back to grid current values. */
template<int TileSize>
class XTileDepositor
{
public:
    
    static const int tileSize = TileSize;

    XTileDepositor(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
    {
        grid = currentDeposition->grid;
        currentDensityCoeff = ptr(currentDeposition->currentDensityCoeff);
        int tileOffset = (tileSize + 1) / 2;
        baseGridIdx = Int3(tileI, 0, 0) - Int3(tileOffset, 0, 0)
            + grid->getNumExternalLeftCells();
        invStep = (FP)1.0  / grid->steps.x;
        normalizedCellOrigin = grid->origin.x / grid->steps.x + (FP)baseGridIdx.x;
        normalizedMiddleCellOrigin = normalizedCellOrigin + (FP)0.5;
        for (int i = 0; i < tileSize; ++i)
        {
            jx[i] = 0;
            jy[i] = 0;
            jz[i] = 0;
        }
    }

    ~XTileDepositor()
    {
        for (int i = 0; i < tileSize; ++i)
         {
            Int3 gridIdx = Int3(i, 0, 0) + baseGridIdx;
            grid->Jx(gridIdx) += jx[i];
            grid->Jy(gridIdx) += jy[i];
            grid->Jz(gridIdx) += jz[i];
        }
    }

    // Same as in TileDepositor
    virtual void depositCurrent(const Particle& particle) = 0;

protected:
    
    FP jx[tileSize];
    FP jy[tileSize];
    FP jz[tileSize];
    Int3 baseGridIdx;
    Grid* grid;
    FP* currentDensityCoeff;
    FP invStep;
    FP normalizedCellOrigin, normalizedMiddleCellOrigin;

};


class XTileDepositorTSC: public XTileDepositor<5>
{
public:

    XTileDepositorTSC(int tileI, int tileJ, int tileK, CurrentDeposition* currentDeposition)
        : XTileDepositor<5>(tileI, tileJ, tileK, currentDeposition),
          halfDt((FP)0.5 * currentDeposition->parameters.timeStep) {}

    void depositCurrent(const Particle& particle)
    {
        FP3 current = particle.p * particle.invGamma
            * currentDensityCoeff[particle.getType()] * particle.getFactor();
        FP normalizedCoords =
            (particle.coords.x - particle.getVelocity().x * halfDt) * invStep;
        FP internalCoords = normalizedCoords - normalizedCellOrigin;
        int idx = (int)(internalCoords + (FP)0.5);
        FP coeff = internalCoords - (FP)idx;
        internalCoords = normalizedCoords - normalizedMiddleCellOrigin;
        int middleIdx = (int)(internalCoords + (FP)0.5);
        FP middleCoeff = internalCoords - (FP)middleIdx;
        depositComponent(current.x, coeff, idx, jx);
        depositComponent(current.y, middleCoeff, middleIdx, jy);
        depositComponent(current.z, middleCoeff, middleIdx, jz);
    }

private:

    void depositComponent(FP value, FP x, int i, FP component[tileSize])
    {
        component[i - 1] += formfactorTSC((FP)1.0 + x) * value;
        component[i] += formfactorTSC(x) * value;
        component[i + 1] += formfactorTSC((FP)1.0 - x) * value;
    }

    const FP halfDt;

};


void CurrentDeposition::run(const Ensemble& ensemble, Grid& grid)
{
    grid.zeroizeJ();
    Grid::DepositionType depositionType = grid.getDepositionType();
    switch (depositionType)
    {
        case Grid::Deposition_CIC:
            addCurrents<TileDepositorCIC>(ensemble, grid);
            break;
        case Grid::Deposition_VillasenorBuneman:
            addCurrents<TileDepositorVillasenorBuneman>(ensemble, grid);
            break;
        case Grid::Deposition_ZigzagFirstOrder:
            addCurrents<TileDepositorZigzagFirstOrder>(ensemble, grid);
            break;
        case Grid::Deposition_TSC:
            if (parameters.globalGridSize.z > 1)
                addCurrents<TileDepositorTSC>(ensemble, grid);
            else
                if (parameters.globalGridSize.y > 1)
                    addCurrents<XYTileDepositorTSC>(ensemble, grid);
                else
                    addCurrents<XTileDepositorTSC>(ensemble, grid);
            break;
        case Grid::Deposition_ZigzagSecondOrder:
            addCurrents<TileDepositorZigzagSecondOrder>(ensemble, grid);
            break;
        case Grid::Deposition_PCS:
            addCurrents<TileDepositorPCS>(ensemble, grid);
            break;
        case Grid::Deposition_Esirkepov:
            if (parameters.globalGridSize.z > 1)
                addCurrents<TileDepositorEsirkepov>(ensemble, grid);
            else
                addCurrents<XYTileDepositorEsirkepov>(ensemble, grid);
            break;
        default:
            throw std::domain_error(
                "Invalid Grid::depositionType not supported in current deposition");
    }
}


template <typename TileDepositorType>
void CurrentDeposition::addCurrents(const Ensemble& ensemble, Grid& _grid)
{
    grid = &_grid; ///
    const int tileSize = TileDepositorType::tileSize;
    for (int ii = 1; ii <= tileSize; ++ii)
    for (int jj = 1; jj <= tileSize; ++jj)
    for (int kk = 1; kk <= tileSize; ++kk)
    {
        Int3 numSystems = ensemble.numSystems;
        #pragma omp parallel for collapse(3) schedule(static, 1)
        for (int i = ii; i < numSystems.x; i += tileSize)
        for (int j = jj; j < numSystems.y; j += tileSize)
        for (int k = kk; k < numSystems.z; k += tileSize)
        {
            const ParticleSystem* system = ensemble.system(i, j, k);
            if (system->size() == 0)
                continue;
            addCurrentsFromCell<TileDepositorType>(i, j, k, ensemble);
        }
    }
}


/* Add currents induced by particles in cell (i, j, k). */
template <typename TileDepositorType>
void CurrentDeposition::addCurrentsFromCell(int i, int j, int k, const Ensemble& ensemble)
{
    TileDepositorType depositor(i, j, k, this);
    const ParticleSystem* system = ensemble.system(i, j, k);
    const int numParticles = system->size();
    const Particle *particles = system->raw();
    for (int idx = 0; idx < numParticles; ++idx)
        depositor.depositCurrent(particles[idx]);
}

} // namespace pica
