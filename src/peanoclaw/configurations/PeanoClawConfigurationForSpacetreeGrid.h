// Copyright (C) 2009 Technische Universitaet Muenchen 
// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef PEANO_APPLICATIONS_PEANOCLAW_CONFIGURATIONSPeanoClawConfiguration_FOR_SPACETREE_GRID_H_
#define PEANO_APPLICATIONS_PEANOCLAW_CONFIGURATIONSPeanoClawConfiguration_FOR_SPACETREE_GRID_H_

#include "tarch/logging/Log.h"

namespace peanoclaw {
  namespace configurations {
      class PeanoClawConfigurationForSpacetreeGrid;
  } 
}
 
 
class peanoclaw::configurations::PeanoClawConfigurationForSpacetreeGrid {

  private:
    /**
     * Logging device.
     */
    static tarch::logging::Log _log;

    bool _isValid;

    bool _plotAtOutputTimes;

    bool _plotSubsteps;

    int _plotSubstepsAfterOutputTime;

    int _additionalLevelsForPredefinedRefinement;

    bool _disableDimensionalSplittingOptimization;

  public: 
    PeanoClawConfigurationForSpacetreeGrid(); 
    virtual ~PeanoClawConfigurationForSpacetreeGrid();
     
    virtual bool isValid() const;
    
    /**
     * Indicates wether the PyClaw output times should be plotted to VTK.
     */
    bool plotAtOutputTimes() const;

    /**
     * Indicates wether the substeps of the adaptive timestepping should be plotted.
     */
    bool plotSubsteps() const;

    /**
     * Turns the substeps plotting on after the given output time has been reached.
     * -1 disables this functionality
     */
    int plotSubstepsAfterOutputTime() const;

    int getAdditionalLevelsForPredefinedRefinement() const;

    bool disableDimensionalSplittingOptimization() const;
};


#endif
