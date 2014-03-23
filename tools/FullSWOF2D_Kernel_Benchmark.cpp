#include <mpi.h>
#include <iostream>
#include <fstream>

#include "tarch/logging/Log.h"

#include "peanoclaw/native/MekkaFlood_solver.h"

#include "peanoclaw/native/FullSWOF2D.h"
#include "peanoclaw/native/FullSWOF_2D/Headers/libschemes/choice_scheme.hpp"


#define DEVELOPMENT

void setupScenario(Scheme* scheme, Parameters& par) {
        /** Water height for breaking dam*/
        TAB& h = scheme->getH();
        for (int x = 0; x < par.get_Nxcell()+2; x++) {
            for (int y = 0; y < par.get_Nycell()+2; y++) {
                h[x][y] = 0.0;
            }
        }
        for (int x = par.get_Nxcell()/2+1; x < par.get_Nxcell()+2; x++) {
            for (int y = par.get_Nycell()/2+1; y < par.get_Nycell()+2; y++) {
                h[x][y] = 1.0;
            }
        }
         
        /** X Velocity.*/
        TAB& u = scheme->getU();
        for (int x = 0; x < par.get_Nxcell()+2; x++) {
            for (int y = 0; y < par.get_Nycell()+2; y++) {
                u[x][y] = 0.0;
            }
        }

        /** Y Velocity.*/
        TAB& v = scheme->getV();
        for (int x = 0; x < par.get_Nxcell()+2; x++) {
            for (int y = 0; y < par.get_Nycell()+2; y++) {
                v[x][y] = 0.0;
            }
        }
         
        /** Topography.*/
        TAB& z = scheme->getZ();
        for (int x = 0; x < par.get_Nxcell()+2; x++) {
            for (int y = 0; y < par.get_Nycell()+2; y++) {
                z[x][y] = 0.0;
            }
        }

        for (int x = par.get_Nxcell()/2+1; x < par.get_Nxcell()+2; x++) {
            for (int y = par.get_Nycell()/2+1; y < par.get_Nycell()+2; y++) {
                z[x][y] = 1.0;
            }
        }
 
}


void setupScenario(const int patchid, int dim, unsigned int *strideinfo, const MekkaFlood_solver::Constants& constants, MekkaFlood_solver::InputArrays& input) {
        /** Water height for breaking dam*/
        for (int x = 0; x < constants.NXCELL+2; x++) {
            for (int y = 0; y < constants.NYCELL+2; y++) {
                    unsigned int index[3];
                    index[0] = y;
                    index[1] = x;
                    index[2] = patchid;
                    unsigned int centerIndex = MekkaFlood_solver::linearizeIndex(3, index, strideinfo);
                
                    input.h[centerIndex] = 0.0;
                }
            }
            for (int x = constants.NXCELL/2+1; x < constants.NXCELL+2; x++) {
                for (int y = constants.NYCELL/2+1; y < constants.NYCELL+2; y++) {
                    unsigned int index[3];
                    index[0] = y;
                    index[1] = x;
                    index[2] = patchid;
                    unsigned int centerIndex = MekkaFlood_solver::linearizeIndex(3, index, strideinfo);
 
                    input.h[centerIndex] = 1.0;
                }
            }
             
            /** X Velocity.*/
            for (int x = 0; x < constants.NXCELL+2; x++) {
                for (int y = 0; y < constants.NYCELL+2; y++) {
                    unsigned int index[3];
                    index[0] = y;
                    index[1] = x;
                    index[2] = patchid;
                    unsigned int centerIndex = MekkaFlood_solver::linearizeIndex(3, index, strideinfo);
 
                    input.u[centerIndex] = 0.0;
                }
            }

            /** Y Velocity.*/
            for (int x = 0; x < constants.NXCELL+2; x++) {
                for (int y = 0; y < constants.NYCELL+2; y++) {
                    unsigned int index[3];
                    index[0] = y;
                    index[1] = x;
                    index[2] = patchid;
                    unsigned int centerIndex = MekkaFlood_solver::linearizeIndex(3, index, strideinfo);
 
                    input.v[centerIndex] = 0.0;
                }
            }
             
            /** Topography.*/
            for (int x = 0; x < constants.NXCELL+2; x++) {
                for (int y = 0; y < constants.NYCELL+2; y++) {
                    unsigned int index[3];
                    index[0] = y;
                    index[1] = x;
                    index[2] = patchid;
                    unsigned int centerIndex = MekkaFlood_solver::linearizeIndex(3, index, strideinfo);

                    input.z[centerIndex] = 0.0;
            }

            for (int x = constants.NXCELL/2+1; x < constants.NXCELL+2; x++) {
                for (int y = constants.NYCELL/2+1; y < constants.NYCELL+2; y++) {
                    unsigned int index[3];
                    index[0] = y;
                    index[1] = x;
                    index[2] = patchid;
                    unsigned int centerIndex = MekkaFlood_solver::linearizeIndex(3, index, strideinfo);
 
                    input.z[centerIndex] = 1.0;
                }
            }
 
        }
}

void compareResults(int nx, int ny, unsigned int *strideinfo, TAB& scheme_unknown, double *kernel_unknown) {
    const int patchid = 0;
        for (int y = 0; y < ny+2; y++) {
            for (int x = 0; x < nx+2; x++) {
               std::cout << " " << scheme_unknown[x][y];
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;

        for (int y = 0; y < ny+2; y++) {
            for (int x = 0; x < nx+2; x++) {
              unsigned int index[3];

               index[0] = y;
               index[1] = x;
               index[2] = patchid;

               unsigned int centerIndex = MekkaFlood_solver::linearizeIndex(3, index, strideinfo);
                
               std::cout << " " << kernel_unknown[centerIndex];
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;
  
        for (int y = 0; y < ny+2; y++) {
            for (int x = 0; x < nx+2; x++) {
              unsigned int index[3];

               index[0] = y;
               index[1] = x;
               index[2] = patchid;

               unsigned int centerIndex = MekkaFlood_solver::linearizeIndex(3, index, strideinfo);
                
               std::cout << " " << std::abs(kernel_unknown[centerIndex]-scheme_unknown[x][y]);
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;

}

void runScenario(int iterations, Choice_scheme* wrapper_scheme, double dtMax) {
        Scheme *scheme = wrapper_scheme->getInternalScheme();

        scheme->resetTimings();
    
        for (int i=0; i < iterations; i++) {
            scheme->setMaxTimestep(dtMax);
            scheme->resetN();
  
            wrapper_scheme->calcul();

            //std::cout << "reference: possible timestep: " << scheme->getTimestep() << std::endl;

#if defined(DEVELOPMENT)
            std::cout << "i=" << i << " reference: current dt " << scheme->getTimestep() << " current dt_max " << scheme->getMaxTimestep() << std::endl;
#endif
        }

        std::cout << std::endl;
}


void runScenario(int iterations, const int patchid, int dim, unsigned int *strideinfo,  MekkaFlood_solver::InputArrays& input, MekkaFlood_solver::TempArrays& temp, const MekkaFlood_solver::Constants& constants, double dtMax) {
    double dt_used = 0.0;
    MekkaFlood_solver solver;
    for (int i=0; i < iterations; i++) {
        dt_used = solver.calcul(patchid, dim, strideinfo, input, temp, constants, dtMax);
        //std::cout << "optimized: possible timestep: " << dt_used << std::endl;
  
#if defined(DEVELOPMENT)
        std::cout << "i= " << i << " optimized: current dt " << dt_used << " current dt_max " << dtMax << std::endl;
#endif
    }
}

void writeTimings(Scheme* scheme, Parameters& par, ofstream& resultfile, const std::string& name) {
        Scheme::timings_t& time_meas = scheme->getTimings();
 
        resultfile << name
                   << " nx " << par.get_Nxcell()
                   << " ny " << par.get_Nycell()
                   << " time " << (time_meas.total_time / time_meas.total_samples)

                   // maincalcflux times
                   << " horizontal_flux " << (time_meas.horizontal_flux_time / time_meas.maincalcflux_samples)
                   << " vertical_flux " << (time_meas.vertical_flux_time / time_meas.maincalcflux_samples)

                   // maincalcscheme times
                   << " rain " << (time_meas.rain_time / time_meas.maincalcscheme_samples)
                   << " mass_conservation " << (time_meas.mass_conservation_time / time_meas.maincalcscheme_samples)
                   << " infiltration " << (time_meas.infiltration_time / time_meas.maincalcscheme_samples)
                   << " momentum " << (time_meas.momentum_time / time_meas.maincalcscheme_samples)
                   << " accumulated_rain " << (time_meas.accumulated_rain_time / time_meas.maincalcscheme_samples)

                   // times used in order 2
                   << " delz " << (time_meas.delz_time / time_meas.delz_samples)
                   << " boundary " << (time_meas.boundary_time / time_meas.boundary_samples)
                   << " reconstruction " << (time_meas.reconstruction_time / time_meas.reconstruction_samples) 
                   << " vincopy " << (time_meas.vincopy_time / time_meas.vincopy_samples) 
                   << " filter " << (time_meas.filter_time / time_meas.filter_samples) 
                   << " heun " <<  (time_meas.heun_time / time_meas.heun_samples) 
                   << std::endl;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    static tarch::logging::Log _log("FullSWOF2D_Kernel_Benchmark");
 
#if defined(DEVELOPMENT)
    const int nr_settings = 1;
    const int cells_x = 8; 
    const int cells_y = 8;
    const int iterations = 20;
#else
    const int nr_settings = 35;
    const int cells_x = 16; 
    const int cells_y = 9;
    const int iterations = 100;
#endif

    double domainSize_x = 4.0;
    double domainSize_y = 3.0;

    // kick of computation
    const int ghostlayerWidth = 1;
    const double meshwidth_x = domainSize_x / cells_x;
    const double meshwidth_y = domainSize_y / cells_y;
    double maximumTimestepSize = 10;
 
    std::ofstream resultfile("result.dat");

    // TODO: run development version here as well:

    // TODO: currently the solver is initialized once and then used
    // in peanoclaw we reconstruct it all the time

    for (size_t i=1; i <= nr_settings; i++) {
        int nx = i * cells_x;
        int ny = i * cells_y;

        // initialize reference solver
        peanoclaw::native::FullSWOF2D_Parameters par(ghostlayerWidth, nx, ny, meshwidth_x, meshwidth_y, 2, 1); // order2 + MUSCL
        Choice_scheme *wrapper_scheme = new Choice_scheme(par);
        Scheme *scheme = wrapper_scheme->getInternalScheme();
  
        // setup and run scenario for solver
        setupScenario(scheme, par);
        runScenario(iterations, wrapper_scheme, maximumTimestepSize);
 
        // write results
        //writeTimings(scheme, par, resultfile, "FullSWOF2D_ref");
  
        // ----------------------------------------------------------------------------------------------------------

#if 0
        // initialize optimized solver
        peanoclaw::native::FullSWOF2D_Parameters par_optimized(ghostlayerWidth, nx, ny, meshwidth_x, meshwidth_y, 3, 4); // order2_opt + MUSCL_OPT
        Choice_scheme *wrapper_scheme_optimized = new Choice_scheme(par_optimized);
        Scheme *scheme_optimized = wrapper_scheme_optimized->getInternalScheme();
           
        // setup and run scenario for optimized solver
        setupScenario(scheme_optimized, par_optimized);
        runScenario(wrapper_scheme_optimized, maximumTimestepSize);

        // write results
        writeTimings(scheme_optimized, par_optimized, resultfile, "FullSWOF2D_opt");
#endif

        // initialize optimized MekkaFlood solver:
        unsigned int strideinfo[3];
        const int nr_patches = 1;
        const int patchid = 0;
        MekkaFlood_solver::InputArrays input;
        MekkaFlood_solver::TempArrays temp;
        MekkaFlood_solver::Constants constants(nx, ny, meshwidth_x, meshwidth_y);

        MekkaFlood_solver::initializeStrideinfo(constants, 3, strideinfo);
        MekkaFlood_solver::allocateInput(nr_patches, 3, strideinfo, input);
        MekkaFlood_solver::allocateTemp(nr_patches, 3, strideinfo, temp);
 
        setupScenario(0, 3, strideinfo, constants, input);
  
        struct timeval start_tv;
        struct timeval stop_tv;

        gettimeofday(&start_tv, NULL);

        runScenario(iterations, 0, 3, strideinfo, input, temp, constants, maximumTimestepSize);
 
        gettimeofday(&stop_tv, NULL);
        double optimized_time = (stop_tv.tv_sec - start_tv.tv_sec) + (stop_tv.tv_usec - start_tv.tv_usec) / 1000000.0;
        optimized_time = optimized_time / iterations;

        std::cout << "nx " << nx << " ny " << ny << " reference " << (scheme->getTimings().total_time / scheme->getTimings().total_samples) << " optimized " << optimized_time << std::endl;

#if defined(DEVELOPMENT)

        std::cout << "comparison: nx " << nx << " ny " << ny << " | " << par.get_Nxcell() << " " << par.get_Nycell() <<  std::endl;

        std::cout << "final h: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getH(), input.h);

        std::cout << "final u: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getU(), input.u);
 
        std::cout << "final v: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getV(), input.v);

        std::cout << "first temp h: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getHs(), temp.hs);
  
        std::cout << "first temp u: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getUs(), temp.us);
  
        std::cout << "first temp v: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getVs(), temp.vs);
 
        std::cout << "second temp h: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getHsa(), temp.hsa);
  
        std::cout << "second temp u: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getUsa(), temp.usa);
  
        std::cout << "second temp v: " << std::endl;
        compareResults(nx, ny, strideinfo, scheme->getVsa(), temp.vsa);
 
 #endif
 
        delete wrapper_scheme;

        MekkaFlood_solver::freeInput(input);
        MekkaFlood_solver::freeTemp(temp);
    }

    MPI_Finalize();

    return 0;
}
