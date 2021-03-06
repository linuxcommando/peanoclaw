# Peano/SConstruct

import os
from os.path import join
from os.path import dirname
import sys
import shutil

#########################################################################
##### FUNCTION DEFINITIONS
#########################################################################

def addPeanoClawFlags(libpath, libs, cpppath, cppdefines):
   ccflags.append('-g3')
   ccflags.append('-g')
   if(environment['PLATFORM'] != 'darwin'):
     ccflags.append('-march=native')
   
   if sys.version_info[0] == 2 and sys.version_info[1] < 7:
       pythonVersion = str(sys.version_info[0]) + '.' + str(sys.version_info[1]) #For Python 2.6
   else: 
       pythonVersion = str(sys.version_info.major) + '.' + str(sys.version_info.minor) #For Python 2.7
   # Determine python version from environment variable:
   peanoClawPythonVersion = os.getenv ('PEANOCLAW_PYTHONVERSION')
   if (peanoClawPythonVersion != None):
      pythonVersion = peanoClawPythonVersion

   # Determine python root path from environment variable:
   pythonHome = os.getenv ('PYTHONHOME')
   peanoClawPythonHome = os.getenv ('PEANOCLAW_PYTHONHOME')
   if (peanoClawPythonHome != None):
      print 'Using environment variable PEANOCLAW_PYTHONHOME =', peanoClawPythonHome
      pythonHome = peanoClawPythonHome
   elif (pythonHome != None):
      print 'Using environment variable PYTHONHOME =', pythonHome
   else:
      print('Environment variables PYTHONHOME and PEANOCLAW_PYTHONHOME not defined. Using path depending on the interpreter\'s path ' + sys.executable)
      pythonHome = join(dirname(sys.executable), '..')
      print pythonHome
   cppdefines.append('NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION')
      
   # Add paths and lib
   libpath.append(pythonHome + '/lib')
   libpath.append(pythonHome + '/lib/python' + pythonVersion)
   cpppath.append(pythonHome + '/include/python' + pythonVersion)
   cpppath.append(pythonHome + '/lib/python' + pythonVersion + '/site-packages/numpy/core/include')
   cpppath.append(os.getenv("HOME") + '/.local/lib/python' + pythonVersion + '/site-packages/numpy/core/include')
   libs.append('python' + pythonVersion)
   if(environment['PLATFORM'] == 'darwin'):
     ccflags.append('-flat_namespace')
     linkerflags.append('-flat_namespace')
   elif build == 'release':
     #cppdefines.append('_GLIBCXX_DEBUG')
     cppdefines.append('NDEBUG')
     
     
   if '-Werror' in ccflags:
     ccflags.remove('-Werror')
     
#########################################################################
##### MAIN CODE
#########################################################################

##### Initialize build variables
#
environment = Environment()
cxx = ''
cppdefines = []
cpppath = ['./src']
ccflags = []
linkerflags = []
libpath = []
libs = []

filenameSuffix = ''

#Configure Peano 3
p3Path = '../peano3'
try:
  import peanoConfiguration
  p3Path = peanoConfiguration.getPeano3Path()
except ImportError:
  pass
p3SourcePath = join(p3Path, 'src')
toolboxSourcePath = join(p3Path, 'toolboxes')
cpppath.append(p3SourcePath)
cpppath.append(toolboxSourcePath)
cpppath.append(join(toolboxSourcePath, 'ControlLoopLoadBalancer'))

if not os.path.isdir(join(toolboxSourcePath, 'ControlLoopLoadBalancer')):
  shutil.copytree('tools/ControlLoopLoadBalancer/ControlLoopLoadBalancer', join(toolboxSourcePath, 'ControlLoopLoadBalancer'))

# Platform specific settings
# Only include library rt if not compiling on Mac OS.
if(environment['PLATFORM'] != 'darwin'):
    libs.append('rt')
    
##### Determine dimension for which to build
#
dim = ARGUMENTS.get('dim', 2)  # Read command line parameter
if int(dim) == 2:
   cppdefines.append('Dim2')
elif int(dim) == 3:
   cppdefines.append('Dim3')
elif int(dim) == 4:
   cppdefines.append('Dim4')
elif int(dim) == 5:
   cppdefines.append('Dim5')
else:
   print "ERROR: dim must be either '2', '3', '4', or '5'!"
   sys.exit(1)

##### Add build parameter specific build variable settings:
# This section only defines Peano-specific flags. It does not
# set compiler specific stuff.
#
build = ARGUMENTS.get('build', 'debug')  # Read command line parameter
if build == 'debug':
   cppdefines.append('Debug')
   cppdefines.append('Asserts')
   cppdefines.append('LogTrace')
   cppdefines.append('LogSeparator')
elif build == 'release':
   pass
elif build == 'asserts':
   cppdefines.append('Asserts')
   pass
else:
   print "ERROR: build must be 'debug', 'asserts', or 'release'!"
   sys.exit(1)
   
##### Determine MPI-Parallelization
#
mpiConfigurationFile = ARGUMENTS.get('mpiconfig', 'openMPIConfiguration')
mpiConfiguration = __import__(mpiConfigurationFile)

parallel = ARGUMENTS.get('parallel', 'parallel_no')  # Read command line parameter
if parallel == 'yes' or parallel == 'parallel_yes':
   cppdefines.append('Parallel')
   cppdefines.append('MPICH_IGNORE_CXX_SEEK')
   cppdefines.append('MPICH_SKIP_MPICXX')
   cpppath.extend(mpiConfiguration.getMPIIncludes())
   libpath.extend(mpiConfiguration.getMPILibrarypaths())
   libs.extend(mpiConfiguration.getMPILibraries())
elif parallel == 'no' or parallel == 'parallel_no':
   pass
else:
   print "ERROR: parallel must be = 'yes', 'parallel_yes', 'no' or 'parallel_no'!"
   sys.exit(1)

##### Determine Multicore usage
#   
multicore = ARGUMENTS.get('multicore', 'multicore_no')  # Read command line parameter

if multicore == 'no' or multicore == 'multicore_no':
   pass
elif multicore == 'openmp':
   ompDir = os.getenv ('OMP_DIR', '')
   cppdefines.append('SharedOMP')
   cpppath.append(ompDir + '/include')   
   pass
elif multicore == 'tbb':
   libs.append('pthread')
   libs.append('dl')
   # Determine tbb directory and architecture from environment variables:
   tbbDir = os.getenv ('TBB_DIR')
          
   libs.append ('tbb')

   cppdefines.append('SharedTBB')
elif multicore == 'opencl':
   libs.append('OpenCL')
   libs.append ('pthread')
   cppdefines.append('SIMD_OpenCL')
else:
   print "ERROR: multicore must be = 'tbb',  'openmp', 'no' or 'multicore_no'!"
   sys.exit(1)

##### Determine Valgrind usage
# 
valgrind = ARGUMENTS.get('valgrind', 'no')
if valgrind == 'no':
   pass
elif valgrind == 'yes':
   ccflags.append('-g')
   cppdefines.append('USE_VALGRIND')
   valgrindRoot = os.getenv ('VALGRIND_ROOT')
   if(valgrindRoot == None):
     valgrindRoot = "/usr"
   cpppath.append(join(valgrindRoot, "include"))
   cpppath.append(join(valgrindRoot, "callgrind"))
else:
   print "ERROR: valgrind must be = 'yes' or 'no'!"
   sys.exit(1)

##### Determine gprof usage
# 
gprof = ARGUMENTS.get('gprof', 'no')
if gprof == 'no':
   pass
elif gprof == 'yes':
   ccflags.append('-pg')
   linkerflags.append('-pg')
else:
   print "ERROR: gprof must be = 'yes' or 'no'!"
   sys.exit(1)
   
##### Switch Compiler
#
compiler = ARGUMENTS.get('compiler', 'gcc')  # Read command line parameter
if compiler == 'gcc':
   if(parallel == 'parallel_no' or parallel == 'no'):
     cxx = 'g++'
   else:
     cxx = 'mpicxx'
     cppdefines.append('MPICH_SKIP_MPICXX')
   ccflags.append('-Wall')
   ccflags.append('-Wstrict-aliasing')
   ccflags.append('-fstrict-aliasing')
   # ccflags.append('-fno-exceptions')
   # ccflags.append('-fno-rtti')
   ccflags.append('-Wno-long-long')
   ccflags.append('-Wno-unknown-pragmas')
   # if multicore == 'no' or multicore == 'multicore_no':
      # ccflags.append('-Wconversion')
   ccflags.append('-Wno-non-virtual-dtor')

   gccversion = environment['CCVERSION'].split('.')
   if int(gccversion[0]) > 4 or int(gccversion[1]) > 6:
     ccflags.append('-std=c++11')
   else:
     ccflags.append('-std=c++0x')
   if build == 'debug':
      ccflags.append('-g3')
      ccflags.append('-O0')
   elif build == 'asserts"':
      ccflags.append('-O2')
      ccflags.append('-g3') 
      ccflags.append('-ggdb')
   elif build == 'release':
      ccflags.append('-O3') 
   if multicore == 'openmp':
      ccflags.append('-fopenmp')
      linkerflags.append('-fopenmp')
elif compiler == 'xlc':
   if(parallel == 'parallel_no' or parallel == 'no'):
     cxx = 'xlc++'
   else:
     cxx = 'mpixlcxx'
   if build == 'debug':
      ccflags.append('-g3')
      ccflags.append('-O0')
   elif build == 'asserts':
      ccflags.append('-qstrict')
      ccflags.append('-O2')
   elif build == 'release':
      ccflags.append('-qstrict')
      ccflags.append('-O3')
   if multicore == 'openmp':
      ccflags.append('-qsmp=omp')
      linkerflags.append('-qsmp=omp')
      cxx = cxx + '_r'
elif compiler == 'icc':
   if(parallel == 'parallel_no' or parallel == 'no'):
     cxx = 'icpc'
   else:
     cxx = 'mpiCC'
   ccflags.append('-fstrict-aliasing')
   ccflags.append('-qpack_semantic=gnu')
   ccflags.append('-std=c++11')
   if build == 'debug':
      ccflags.append('-O0')
   elif build == 'asserts':
      ccflags.append('-w')
      ccflags.append('-align')
      ccflags.append('-ansi-alias')
      ccflags.append('-O2')
   elif build == 'release':
      ccflags.append('-w')
      ccflags.append('-align')
      ccflags.append('-ansi-alias')
      ccflags.append('-O3')
   if multicore == 'openmp':
      ccflags.append('-openmp')
      linkerflags.append('-openmp')
else:
   print "ERROR: compiler must be = 'gcc', 'xlc' or 'icc'!"
   sys.exit(1)
   
##### Determine Scalasca Usage
#
scalasca = ARGUMENTS.get('scalasca', 'scalasca_no')  # Read command line parameter
if scalasca == 'yes' or scalasca == 'scalasca_yes':
   cxx = 'scalasca -instrument ' + cxx
elif scalasca == 'no' or scalasca == 'scalasca_no':
   pass
else:
   print "ERROR: scalasca must be = 'scalasca_yes', 'yes', 'scalasca_no' or 'no'!"
   sys.exit(1)
   
##### Determine Solver
#
solver = ARGUMENTS.get('solver', 'pyclaw')
if solver == 'pyclaw':
  cppdefines.append('PYCLAW')
  cppdefines.append('PEANOCLAW_PYCLAW')
  cppdefines.append('AssertForPositiveValues')
elif solver == 'swe':
  #Configure SWE-Sources
  swePath = '../SWE/src'
  try:
    import sweConfiguration
    swePath = sweConfiguration.getSWEPath()
  except ImportError:
    pass
  cpppath.append(swePath)
  cppdefines.append('SWE')
  cppdefines.append('NDEBUG')
  
  cppdefines.append('WAVE_PROPAGATION_SOLVER=1')
  cppdefines.append('VECTORIZE')
  
  cppdefines.append('AssertForPositiveValues')
elif solver == 'fullswof2d':
  #Configure FullSWOF-Sources
  fullSWOF2DPath = '../FullSWOF_2D'
  try:
    import fullSWOF2DConfiguration
    fullSWOF2DPath = fullSWOF2DConfiguration.getFullSWOF2DPath()
  except ImportError:
    pass
  cpppath.append(fullSWOF2DPath)
  cppdefines.append('NDEBUG')

  cppdefines.append('WAVE_PROPAGATION_SOLVER=1')
  cppdefines.append('VECTORIZE')
  cppdefines.append('DoNotAssertForPositiveValues')

  #Configure fullswof2d
  # add FULLSWOF2D mode
  cppdefines.append('PEANOCLAW_FULLSWOF2D')
  cpppath.append( join(fullSWOF2DPath, 'Headers/liblimitations') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libfrictions') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libparser') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libflux') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libsave') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libschemes') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libreconstructions') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libinitializations') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/librain_infiltration') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libboundaryconditions') )
  cpppath.append( join(fullSWOF2DPath, 'Headers/libparameters') )
  libs.append('png') # for texture file
else:
  raise Exception("ERROR: solver must be 'pyclaw', 'swe', or 'fullswof2d'")

##### Determine Heap Compression
#
heapCompression = ARGUMENTS.get('heapCompression', 'yes')
if heapCompression == 'no':
  cppdefines.append('noPackedEmptyHeapMessages')
  filenameSuffix += '_noHeapCompression'
elif heapCompression == 'yes':
  pass
else:
  raise Exception("ERROR: heapCompression must be 'yes' or 'no'")
   
##### Determine build path
#
build_offset = ARGUMENTS.get('buildoffset', 'build')
buildpath = build_offset + '/' + str(build) + '/dim' + str(dim) + '/' 
if multicore == 'tbb':
   buildpath = join(buildpath, 'tbb')
elif multicore == 'openmp':
   buildpath = join(buildpath, 'openmp')
elif multicore == 'opencl':
   buildpath = join(buildpath, 'openCL')
else:
   buildpath = join(buildpath, 'multicore_no')
if parallel == 'yes' or parallel == 'parallel_yes':
   buildpath = join(buildpath, 'parallel_yes')
else:
   buildpath = join(buildpath, 'parallel_no')
buildpath = join(buildpath, compiler)
buildpath = join(buildpath, solver)
if scalasca == 'yes' or scalasca == 'scalasca_yes':
   buildpath = join(buildpath, 'scalasca')
if heapCompression == 'no':
   buildpath = join(buildpath, 'noHeapCompression')
if gprof == 'yes':
  buildpath = join(buildpath, 'gprof') 

buildpath = buildpath + '/'
   
##### Specify build settings
#
addPeanoClawFlags(libpath, libs, cpppath, cppdefines)

##### Print options used to build
#
print
print "Building PeanoClaw"
print "Options: build = " + str(build) + ", dim = " + str(dim) + ", build-offset = " + str(build_offset) + ", parallel = " + str(parallel) + ", multicore = " + str(multicore) + ", compiler = " + str(compiler)
print "Buildpath: " + buildpath
print

VariantDir (buildpath, './src', duplicate=0)  # Set build directory for PeanoClaw sources
VariantDir (join(buildpath, 'kernel'), p3SourcePath, duplicate=0)  # Set build directory for Peano sources
if solver == 'swe':
  VariantDir (join(buildpath, 'swe'), swePath, duplicate=0)  # Set build directory for SWE sources
if solver == 'fullswof2d':
  
  print "VariantDir", join(buildpath, 'fullswof2d'), fullSWOF2DPath
  
  VariantDir (join(buildpath, 'fullswof2d'), fullSWOF2DPath, duplicate=0)  # Set build directory for SWE sources
  
  
##### Setup construction environment:
#
env = Environment (
   CPPDEFINES=cppdefines,
   LIBPATH=libpath,
   LIBS=libs,
   CPPPATH=cpppath,
   CCFLAGS=ccflags,
   LINKFLAGS=linkerflags,
   CXX=cxx,
   ENV=os.environ  # Makes environment variables visible to scons
   # tools      = compiler_tools
   )

################################################################################
#
# Define sources
#

##### Sub T-components

sourcesTLa = [
   Glob(join(buildpath, 'kernel/tarch/la/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/la/tests/*.cpp'))
   ]

sourcesTCompiler = [
   Glob(join(buildpath, 'kernel/tarch/compiler/*.cpp'))
   ]

sourcesTConfiguration = [
   Glob(join(buildpath, 'kernel/tarch/configuration/*.cpp'))
   ]

sourcesTIrr = [
   Glob(join(buildpath, 'kernel/tarch/irr/*.cpp'))
 ]

sourcesTLogging = [
  Glob(join(buildpath, 'kernel/tarch/logging/*.cpp')),
  Glob(join(buildpath, 'kernel/tarch/logging/configurations/*.cpp'))
]

sourcesTServices = [
  Glob(join(buildpath, 'kernel/tarch/services/*.cpp'))
]

sourcesTTests = [
  Glob(join(buildpath, 'kernel/tarch/tests/*.cpp')),
  Glob(join(buildpath, 'kernel/tarch/tests/configurations/*.cpp'))
  ]

sourcesTUtils = [
  Glob(join(buildpath, 'kernel/tarch/utils/*.cpp'))
]

sourcesTTiming = [
  Glob(join(buildpath, 'kernel/tarch/timing/*.cpp'))
]

sourcesTPlotter = [ 
   Glob(join(buildpath, 'kernel/tarch/plotter/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/plotter/globaldata/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/plotter/globaldata/tests/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/plotter/griddata/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/plotter/griddata/multiscale/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/plotter/griddata/unstructured/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/plotter/griddata/unstructured/configurations/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/plotter/griddata/unstructured/vtk/*.cpp')),
   Glob(join(buildpath, 'kernel/tarch/plotter/griddata/unstructured/vtk/tests/*.cpp'))
   ]   
   
##### Define sources T-components
#            
sourcesTComponents = [
   sourcesTCompiler,
   sourcesTConfiguration,
   sourcesTIrr,
   sourcesTLa,
   sourcesTLogging,
   sourcesTPlotter,
   sourcesTServices,
   sourcesTTests,
   sourcesTTiming,
   sourcesTUtils
   ]


##### Define sources for multicore support
#    
sourcesDatatraversal = [
    Glob(join(buildpath, 'kernel/peano/datatraversal/*.cpp')),
    Glob(join(buildpath, 'kernel/peano/datatraversal/configurations/*.cpp')),
    Glob(join(buildpath, 'kernel/peano/datatraversal/tests/*.cpp')),
    Glob(join(buildpath, 'kernel/peano/datatraversal/autotuning/*.cpp')),
    Glob(join(buildpath, 'kernel/tarch/multicore/configurations/*.cpp')),
    Glob(join(buildpath, 'kernel/tarch/multicore/*.cpp'))
  ]       
      
if multicore == 'no' or multicore == 'multicore_no':
   pass
elif multicore == 'openmp':
   sourcesDatatraversal = sourcesDatatraversal + [
     Glob(join(buildpath, 'kernel/tarch/multicore/openMP/*.cpp'))
   ]
elif multicore == 'tbb':
   sourcesDatatraversal = sourcesDatatraversal + [
     Glob(join(buildpath, 'kernel/tarch/multicore/tbb/*.cpp'))
   ]

sourcesParallel = [
 Glob(join(buildpath, 'kernel/tarch/parallel/configuration/*.cpp')),
 Glob(join(buildpath, 'kernel/tarch/parallel/*.cpp')),
 Glob(join(buildpath, 'kernel/tarch/parallel/strategy/*.cpp')),
 Glob(join(buildpath, 'kernel/tarch/parallel/messages/*.cpp')),
 Glob(join(buildpath, 'kernel/tarch/parallel/dastgen/*.cpp')),
 Glob(join(buildpath, 'kernel/tarch/parallel/configurations/*.cpp')),
 Glob(join(buildpath, 'kernel/peano/parallel/*.cpp')),
 Glob(join(buildpath, 'kernel/peano/parallel/*.cc')),
 Glob(join(buildpath, 'kernel/peano/parallel/configurations/*.cpp')),
 Glob(join(buildpath, 'kernel/peano/parallel/loadbalancing/*.cpp')),
 Glob(join(buildpath, 'kernel/peano/parallel/messages/*.cpp')),
 Glob(join(buildpath, 'kernel/peano/parallel/tests/*.cpp')),
 Glob(join(buildpath, 'kernel/tarch/mpianalysis/*.cpp')),
 Glob(join(buildpath, 'kernel/tarch/analysis/*.cpp'))
]


#### Peano Utils
sourcesPeanoUtils = [
  Glob(join(buildpath, 'kernel/peano/utils/*.cpp'))
]


# ## Peano partition coupling
sourcesPartitionCoupling = [
  Glob(join(buildpath, 'kernel/peano/integration/partitioncoupling/*.cpp')),
  Glob(join(buildpath, 'kernel/peano/integration/partitioncoupling/builtin/*.cpp')),
  Glob(join(buildpath, 'kernel/peano/integration/partitioncoupling/builtin/configurations/*.cpp')),
  Glob(join(buildpath, 'kernel/peano/integration/partitioncoupling/builtin/tests/*.cpp')),
  Glob(join(buildpath, 'kernel/peano/integration/partitioncoupling/builtin/records/*.cpp')),
  Glob(join(buildpath, 'kernel/peano/integration/partitioncoupling/services/*.cpp'))
]

# ## Kernel
sourcesKernelConfiguration = [
   Glob(join(buildpath, 'kernel/peano/configurations/*.cpp'))
   ]

sourcesGridInterface = [
   Glob(join(buildpath, 'kernel/peano/gridinterface/*.cpp'))
   ]
   
sourcesGrid = [
   Glob(join(buildpath, 'kernel/peano/grid/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/grid/aspects/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/grid/nodes/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/grid/nodes/loops/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/grid/nodes/tasks/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/grid/tests/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/grid/tests/records/*.cpp'))
   ]
sourcesStacks = [
    Glob(join(buildpath, 'kernel/peano/stacks/*.cpp')),
    Glob(join(buildpath, 'kernel/peano/stacks/implementation/*.cpp'))
    ]
sourcesHeap = [
    Glob(join(buildpath, 'kernel/peano/heap/records/*.cpp'))
    ]
sourcesPeanoKernel = [
   sourcesKernelConfiguration,
   sourcesGridInterface,
   sourcesGrid,
   sourcesStacks,
   sourcesHeap
   ]

#### Geometry
##### Builtin Geometry
sourcesBuiltinGeometry = [
   Glob(join(buildpath, 'kernel/peano/geometry/builtin/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/geometry/builtin/services/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/geometry/extensions/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/geometry/builtin/configurations/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/geometry/builtin/tests/*.cpp'))
   ]


sourcesPeanoGeometry = [
   Glob(join(buildpath, 'kernel/peano/geometry/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/geometry/tests/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/geometry/configurations/*.cpp')),
   Glob(join(buildpath, 'kernel/peano/geometry/services/*.cpp')),
   sourcesBuiltinGeometry
   ]

sourcesPeanoBase = [
  sourcesPeanoKernel,
  sourcesPeanoGeometry,
  sourcesPeanoUtils,
  sourcesDatatraversal,
  # sourcesQueries,
  Glob(join(buildpath, 'kernel/peano/*.cpp')),
  Glob(join(buildpath, 'kernel/*.cpp'))
]

sourcesToolBox = [
  Glob(join(buildpath, 'toolboxes/ControlLoopLoadBalancer/*.cpp')),
  Glob(join(buildpath, 'toolboxes/ControlLoopLoadBalancer/strategies/*.cpp'))
]

sourcesToolBoxVHH = [
  Glob(join(buildpath, 'kernel/peano/toolbox/solver/vhh/*.cpp')),
  Glob(join(buildpath, 'kernel/peano/toolbox/solver/vhh/tests/*.cpp')),
  Glob(join(buildpath, 'kernel/peano/toolbox/solver/vhh/opencl/*.cpp')),
  Glob(join(buildpath, 'kernel/peano/toolbox/solver/vhh/opencl/tests/*.cpp'))
]

# ## Applications
sourcesPeanoClaw = [
  Glob(join(buildpath, 'peanoclaw/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/adapters/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/configurations/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/grid/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/interSubgridCommunication/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/interSubgridCommunication/aspects/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/mappings/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/parallel/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/records/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/repositories/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/runners/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/statistics/*.cpp')),
  Glob(join(buildpath, 'peanoclaw/tests/*.cpp')),
	]

##### Define sources of application peanoclaw
if solver == 'swe':
  sourcesSolver = [
    Glob(join(buildpath, 'peanoclaw/native/main.cpp')),
    Glob(join(buildpath, 'peanoclaw/native/SWEKernel.cpp')),
    Glob(join(buildpath, 'peanoclaw/native/SWE_WavePropagationBlock_patch.cpp')),
    Glob(join(buildpath, 'peanoclaw/native/BreakingDam.cpp')),
    Glob(join(buildpath, 'swe/blocks/SWE_Block.cpp')),
    #Glob(join(buildpath, 'swe/blocks/SWE_WaveAccumulationBlock.cpp'))
    Glob(join(buildpath, 'swe/blocks/SWE_WavePropagationBlock.cpp'))
    ]
elif solver == 'pyclaw':
  sourcesSolver = [
     Glob(join(buildpath, 'peanoclaw/pyclaw/*.cpp'))
     ]
elif solver == 'fullswof2d':
  sourcesSolver = [
     Glob(join(buildpath, 'peanoclaw/native/FullSWOF2D.cpp')),
     Glob(join(buildpath, 'peanoclaw/native/MekkaFlood_solver.cpp')),
     Glob(join(buildpath, 'peanoclaw/native/MekkaFlood.cpp')),
     Glob(join(buildpath, 'peanoclaw/native/BreakingDam.cpp')),
     Glob(join(buildpath, 'peanoclaw/native/dem.cpp')),
     Glob(join(buildpath, 'peanoclaw/native/main.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/liblimitations/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libfrictions/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libparser/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libflux/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libsave/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libschemes/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libreconstructions/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libinitializations/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/librain_infiltration/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libboundaryconditions/*.cpp')),
     Glob(join(buildpath, 'fullswof2d/Sources/libparameters/*.cpp'))
     ]
sourcesPeanoClaw.extend(sourcesSolver)

################################################################################

##### Configure
configure = Configure(env)
#if configure.CheckCXXHeader('peano/parallel/MeshCommunication.h'):
if os.path.isfile(join(p3SourcePath, 'peano/parallel/MeshCommunication.h')):
  print 'Using RMK'
  env['CPPDEFINES'].append('UseBlockedMeshCommunication')
else:
  print 'Using Peano classic communication'
  env['CPPDEFINES'].append('DoNotUseBlockedMeshCommunication')
#env = configure.Finish()

##### Build selected target
#
source = [
   sourcesTComponents,
   sourcesPeanoBase,
   sourcesPeanoClaw,
   sourcesParallel,
   sourcesToolBox
   ]

if solver == 'pyclaw':
  targetfilename = 'libpeano-claw-' + str(dim) + 'd' + filenameSuffix
  target = buildpath + targetfilename
  library = env.SharedLibrary (
    target=target,
    source=source
    )
    
  ##### Copy library to Clawpack
  #
  installation = env.Alias('install', env.Install('src/python/peanoclaw', library))
elif solver == 'swe':
  targetfilename = 'peano-claw-' + str(dim) + 'd'
  target = buildpath + targetfilename
  executable = env.Program ( 
    target=target,
    source=source
    )
  ##### Copy executable to bin directory
  #
  installation = env.Alias('install', env.Install('bin', executable))    
elif solver == 'fullswof2d':
  targetfilename = 'peano-claw-' + str(dim) + 'd'
  target = buildpath + targetfilename
  executable = env.Program ( 
    target=target,
    source=source
    )
  ##### Copy executable to bin directory
  #
  installation = env.Alias('install', env.Install('bin', executable))

Default(installation)

