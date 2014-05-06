/*
 * NativeKernel.cpp
 *
 *  Created on: Jun 24, 2013
 *      Author: kristof
 */
#include "peanoclaw/native/NativeKernel.h"

#include "peanoclaw/Patch.h"
#include "peanoclaw/interSubgridCommunication/DefaultInterpolation.h"
#include "peanoclaw/interSubgridCommunication/DefaultRestriction.h"
#include "peanoclaw/interSubgridCommunication/DefaultFluxCorrection.h"

peanoclaw::native::NativeKernel::NativeKernel()
: Numerics(
    new peanoclaw::interSubgridCommunication::DefaultInterpolation,
    new peanoclaw::interSubgridCommunication::DefaultRestriction,
    new peanoclaw::interSubgridCommunication::DefaultFluxCorrection
  )
{

}


void peanoclaw::native::NativeKernel::addPatchToSolution(Patch& patch) {
}

void peanoclaw::native::NativeKernel::initializePatch(Patch& patch) {
}

void peanoclaw::native::NativeKernel::fillBoundaryLayer(
  Patch& patch,
  int dimension,
  bool setUpper
) {
}

void peanoclaw::native::NativeKernel::solveTimestep(
  Patch& patch,
  double maximumTimestepSize,
  bool useDimensionalSplitting
) {

}

tarch::la::Vector<DIMENSIONS, double> peanoclaw::native::NativeKernel::getDemandedMeshWidth(Patch& patch, bool isInitializing) {
  return patch.getSubcellSize()(0);
}
