/*
 * PyClaw.cpp
 *
 *  Created on: Feb 18, 2012
 *      Author: kristof
 */

#include "peanoclaw/native/SWEKernel.h"
#include "peanoclaw/Patch.h"
#include "tarch/timing/Watch.h"
#include "tarch/parallel/Node.h"

#include "peanoclaw/native/SWE_WavePropagationBlock_patch.hh"

tarch::logging::Log peanoclaw::native::SWEKernel::_log("peanoclaw::native::SWEKernel");

peanoclaw::native::SWEKernel::SWEKernel(
  SWEKernelScenario& scenario,
  peanoclaw::interSubgridCommunication::Interpolation*  interpolation,
  peanoclaw::interSubgridCommunication::Restriction*    restriction,
  peanoclaw::interSubgridCommunication::FluxCorrection* fluxCorrection
) : Numerics(interpolation, restriction, fluxCorrection),
_scenario(scenario),
_totalSolverCallbackTime(0.0)
{
  //import_array();
}

peanoclaw::native::SWEKernel::~SWEKernel()
{
}


void peanoclaw::native::SWEKernel::initializePatch(
  Patch& patch
) {
  logTraceIn( "initializePatch(...)");

  _scenario.initializePatch(patch);

  logTraceOutWith1Argument( "initializePatch(...)", demandedMeshWidth);
}

void peanoclaw::native::SWEKernel::solveTimestep(Patch& patch, double maximumTimestepSize, bool useDimensionalSplitting) {
  logTraceInWith2Arguments( "solveTimestep(...)", maximumTimestepSize, useDimensionalSplitting);

  assertion2(tarch::la::greater(maximumTimestepSize, 0.0), "Timestepsize == 0 should be checked outside.", patch.getTimeIntervals().getMinimalNeighborTimeConstraint());

  tarch::timing::Watch pyclawWatch("", "", false);
  pyclawWatch.startTimer();
  double dtAndEstimatedNextDt[2];

  SWE_WavePropagationBlock_patch swe(patch);

  swe.computeNumericalFluxes();

  double dt = fmin(swe.getMaxTimestep(), maximumTimestepSize);
  double estimatedNextTimestepSize = swe.getMaxTimestep();

  swe.updateUnknowns(dt);

  pyclawWatch.stopTimer();
  _totalSolverCallbackTime += pyclawWatch.getCalendarTime();

  assertion4(
      tarch::la::greater(patch.getTimeIntervals().getTimestepSize(), 0.0)
      || tarch::la::greater(estimatedNextTimestepSize, 0.0)
      || tarch::la::equals(maximumTimestepSize, 0.0)
      || tarch::la::equals(patch.getTimeIntervals().getEstimatedNextTimestepSize(), 0.0),
      patch, maximumTimestepSize, estimatedNextTimestepSize, patch.toStringUNew());
  assertion(patch.getTimeIntervals().getTimestepSize() < std::numeric_limits<double>::infinity());

  if (tarch::la::greater(dt, 0.0)) {
    patch.getTimeIntervals().advanceInTime();
    patch.getTimeIntervals().setTimestepSize(dt);
  }
  patch.getTimeIntervals().setEstimatedNextTimestepSize(estimatedNextTimestepSize);

  logTraceOut( "solveTimestep(...)");
}

tarch::la::Vector<DIMENSIONS, double> peanoclaw::native::SWEKernel::getDemandedMeshWidth(Patch& patch, bool isInitializing) {
  return _scenario.computeDemandedMeshWidth(patch, isInitializing);
}

void peanoclaw::native::SWEKernel::addPatchToSolution(Patch& patch) {
}

void peanoclaw::native::SWEKernel::fillBoundaryLayer(Patch& patch, int dimension, bool setUpper) {
  logTraceInWith3Arguments("fillBoundaryLayerInPyClaw", patch, dimension, setUpper);

  logDebug("fillBoundaryLayerInPyClaw", "Setting left boundary for " << patch.getPosition() << ", dim=" << dimension << ", setUpper=" << setUpper);

   //std::cout << "------ setUpper " << setUpper << " dimension " << dimension << std::endl;
   //std::cout << patch << std::endl;
   //std::cout << "++++++" << std::endl;
   //std::cout << patch.toStringUOldWithGhostLayer() << std::endl;
   //std::cout << "||||||" << std::endl;

   // implement a wall boundary
    tarch::la::Vector<DIMENSIONS, int> src_subcellIndex;
    tarch::la::Vector<DIMENSIONS, int> dest_subcellIndex;

    if (dimension == 0) {
        for (int yi = -1; yi < patch.getSubdivisionFactor()(1)+1; yi++) {
            int xi = setUpper ? patch.getSubdivisionFactor()(0) : -1;
            src_subcellIndex(0) = xi;
            src_subcellIndex(1) = yi;
            src_subcellIndex(dimension) += setUpper ? -1 : +1; 

            dest_subcellIndex(0) = xi;
            dest_subcellIndex(1) = yi;
     
            for (int unknown=0; unknown < patch.getUnknownsPerSubcell(); unknown++) {
                double q = patch.getValueUOld(src_subcellIndex, unknown);

                if (unknown == dimension + 1) {
                    patch.setValueUOld(dest_subcellIndex, unknown, -q);
                } else {
                    patch.setValueUOld(dest_subcellIndex, unknown, q);
                }
            }
        }

    } else {
        for (int xi = -1; xi < patch.getSubdivisionFactor()(0)+1; xi++) {
            int yi = setUpper ? patch.getSubdivisionFactor()(1) : -1;
            src_subcellIndex(0) = xi;
            src_subcellIndex(1) = yi;
            src_subcellIndex(dimension) += setUpper ? -1 : +1; 

            dest_subcellIndex(0) = xi;
            dest_subcellIndex(1) = yi;
     
            for (int unknown=0; unknown < patch.getUnknownsPerSubcell(); unknown++) {
                double q = patch.getValueUOld(src_subcellIndex, unknown);

                if (unknown == dimension + 1) {
                    patch.setValueUOld(dest_subcellIndex, unknown, -q);
                } else {
                    patch.setValueUOld(dest_subcellIndex, unknown, q);
                }
            }
        }
    }

   //std::cout << "++++++" << std::endl;
   //std::cout << patch.toStringUOldWithGhostLayer() << std::endl;
   //std::cout << "||||||" << std::endl;

      logTraceOut("fillBoundaryLayerInPyClaw");
    }

