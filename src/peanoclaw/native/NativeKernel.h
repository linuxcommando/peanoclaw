/*
 * NativeKernel.h
 *
 *  Created on: Jun 24, 2013
 *      Author: kristof
 */

#ifndef PEANOCLAW_NATIVE_NATIVEKERNEL_H_
#define PEANOCLAW_NATIVE_NATIVEKERNEL_H_

#include "peanoclaw/Numerics.h"

namespace peanoclaw {

  class Patch;

  namespace native {
    class NativeKernel;
  }
}

class peanoclaw::native::NativeKernel : public peanoclaw::Numerics {

  public:

    NativeKernel();

    /**
     * @see peanoclaw::Numerics
     */
    void addPatchToSolution(Patch& patch);

    /**
     * @see peanoclaw::Numerics
     */
    void initializePatch(Patch& patch);

    /**
     * @see peanoclaw::Numerics
     */
    void fillBoundaryLayer(
      Patch& patch,
      int dimension,
      bool setUpper
    );

    /**
     * @see peanoclaw::Numerics
     */
    void solveTimestep(
      Patch& patch,
      double maximumTimestepSize,
      bool useDimensionalSplitting
    );

    tarch::la::Vector<DIMENSIONS, double> getDemandedMeshWidth(
        Patch& patch,
        bool isInitializing
    );

};

#endif /* PEANOCLAW_NATIVE_NATIVEKERNEL_H_ */
