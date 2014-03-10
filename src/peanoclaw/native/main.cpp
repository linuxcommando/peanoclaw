/*
 * pyclawBindings.c
 *
 *  Created on: Feb 7, 2012
 *      Author: Kristof Unterweger
 */

#include "peano/utils/Globals.h"
#include "tarch/logging/Log.h"
#include "tarch/logging/CommandLineLogger.h"

#include "peano/peano.h"

#include <list>

#include "peanoclaw/Patch.h"
#include "peanoclaw/Numerics.h"
#include "peanoclaw/NumericsFactory.h"
#include "peanoclaw/configurations/PeanoClawConfigurationForSpacetreeGrid.h"
#include "peanoclaw/runners/PeanoClawLibraryRunner.h"
#include "tarch/tests/TestCaseRegistry.h"

#include "peanoclaw/native/SWEKernel.h"

#if USE_VALGRIND
#include <callgrind.h>
#endif

#if defined(SWE)
#include "peanoclaw/native/BreakingDam.h"
#endif

#if defined(PEANOCLAW_FULLSWOF2D)
#include "peanoclaw/native/MekkaFlood.h"
#include "peanoclaw/native/dem.h"
#endif

static peanoclaw::configurations::PeanoClawConfigurationForSpacetreeGrid* _configuration;

/*void importArrays() {
  import_array();
}*/

int main(int argc, char **argv) {
  peano::fillLookupTables();


#if defined(Parallel)
  int parallelSetup = peano::initParallelEnvironment(&argc,(char ***)&argv);
  int sharedMemorySetup = peano::initSharedMemoryEnvironment();
#endif

  //importArrays();

  //Initialize Logger
  static tarch::logging::Log _log("peanoclaw");
  logInfo("main(...)", "Initializing Peano");

  // Configure the output
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", true ) );
  //tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", true ) );
//  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "trace", true ) );

  //Validation
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", -1, "peanoclaw::statistics::ParallelGridValidator", true ) );

  //Selective Tracing
//  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "peanoclaw::mappings::Remesh", false ) );
//  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "peanoclaw::mappings::Remesh::destroyVertex", false ) );
//  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "peanoclaw::mappings::Remesh::endIteration", false ) );
//  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "peanoclaw::mappings::Remesh::touchVertex", false ) );

  //tarch::logging::CommandLineLogger::getInstance().setLogFormat( ... please consult source code documentation );

  std::ostringstream logFileName;
  #ifdef Parallel
  logFileName << "rank-" << tarch::parallel::Node::getInstance().getRank() << "-trace.txt";
  #endif
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, logFileName.str() );

  //Tests
  if(false) {
    tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().run();
  }

  //PyClaw - this object is copied to the runner and is stored there.
  peanoclaw::NumericsFactory numericsFactory;

#if defined(SWE) || defined(PEANOCLAW_FULLSWOF2D)
  _configuration = new peanoclaw::configurations::PeanoClawConfigurationForSpacetreeGrid;
  // assertion1(_configuration->isValid(), _configuration);

  //Construct parameters
 
#if defined(PEANOCLAW_FULLSWOF2D)
    DEM dem;

    dem.load("DEM_050cm.bin");

    MekkaFlood_SWEKernelScenario scenario(dem);
    peanoclaw::Numerics* numerics = numericsFactory.createFullSWOF2DNumerics(scenario);
 
    tarch::la::Vector<DIMENSIONS, double> domainOffset;

    // TODO: aaarg Y U NO PLOT CORRECTLY! -> work around established
    domainOffset(0) = 0.0; //dem.lower_left(0);
    domainOffset(1) = 0.0; //dem.lower_left(1);

    tarch::la::Vector<DIMENSIONS, double> domainSize;
    double upper_right_0 = dem.upper_right(0);
    double upper_right_1 = dem.upper_right(1);
 
    double lower_left_0 = dem.lower_left(0);
    double lower_left_1 = dem.lower_left(1);
 
    double x_size = upper_right_0 - lower_left_0;
    double y_size = upper_right_1 - lower_left_1;
 
    domainSize(0) = x_size;
    domainSize(1) = y_size;

    // keep aspect ratio of map: 4000 3000: ratio 4:3
    tarch::la::Vector<DIMENSIONS, int> subdivisionFactor;
    subdivisionFactor(0) = static_cast<int>(24); //  6 * 4
    subdivisionFactor(1) = static_cast<int>(18); //  6 * 3

    double min_domainSize = std::min(domainSize(0),domainSize(1));
    int min_subdivisionFactor = std::min(subdivisionFactor(0),subdivisionFactor(1));

    tarch::la::Vector<DIMENSIONS, double> initialMinimalMeshWidth(min_domainSize/min_subdivisionFactor);

    int ghostlayerWidth = 2;
    int unknownsPerSubcell = 6;

    int initialTimestepSize = 10.0;


#else
    BreakingDam_SWEKernelScenario scenario;
    peanoclaw::Numerics* numerics = numericsFactory.createSWENumerics(scenario);
 
    tarch::la::Vector<DIMENSIONS, double> domainOffset(0);
    tarch::la::Vector<DIMENSIONS, double> domainSize(1.0);
    tarch::la::Vector<DIMENSIONS, double> initialMinimalMeshWidth(domainSize(0)/6/27);

    int ghostlayerWidth = 1;
    int unknownsPerSubcell = 3;

    int initialTimestepSize = 0.1;
 
    tarch::la::Vector<DIMENSIONS, int> subdivisionFactor(6);
#endif
  
  //tarch::la::Vector<DIMENSIONS, double> initialMinimalMeshWidth(10.0/130/27);
  //tarch::la::Vector<DIMENSIONS, int> subdivisionFactor(130);
  int auxiliarFieldsPerSubcell = 1;
  bool useDimensionalSplittingOptimization = true;
 

  //Check parameters
  assertion1(tarch::la::greater(domainSize(0), 0.0) && tarch::la::greater(domainSize(1), 0.0), domainSize);
  if(initialMinimalMeshWidth(0) > domainSize(0) || initialMinimalMeshWidth(1) > domainSize(1)) {
    logError("main(...)", "Domainsize or initialMinimalMeshWidth not set properly.");
  }
  if(tarch::la::oneGreater(tarch::la::Vector<DIMENSIONS, int>(1), subdivisionFactor(0)) ) {
    logError("main(...)", "subdivisionFactor not set properly.");
  }
 
  //Create runner
  peanoclaw::runners::PeanoClawLibraryRunner* runner
    = new peanoclaw::runners::PeanoClawLibraryRunner(
    *_configuration,
    *numerics,
    domainOffset,
    domainSize,
    initialMinimalMeshWidth,
    subdivisionFactor,
    ghostlayerWidth,
    unknownsPerSubcell,
    auxiliarFieldsPerSubcell,
    initialTimestepSize,
    useDimensionalSplittingOptimization,
    1
  );

#if defined(Parallel) 
  std::cout << tarch::parallel::Node::getInstance().getRank() << ": peano instance created" << std::endl;
#endif

  assertion(runner != 0);
 
  // run experiment
  double timestep = 10.0;
  double endtime = 10000.0; //1.0; //2.0;
#if defined(Parallel)
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
#endif
      for (double time=0.0; time < endtime; time+=timestep) {
        runner->evolveToTime(time);
        //runner->gatherCurrentSolution();
        std::cout << "time " << time << " done " << std::endl;
      }
 
      /*for (int i=0; i < 20; ++i) {
          runner->runNextPossibleTimestep();
      }*/

#if defined(Parallel)
  } else {
    runner->runWorker();
  }
#endif

  // experiment done -> cleanup

  delete runner;
 
  if(_configuration != 0) {
    delete _configuration;
  }
#endif
  return 0;
}
