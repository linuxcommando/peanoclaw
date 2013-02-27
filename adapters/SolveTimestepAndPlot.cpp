#include "peanoclaw/adapters/SolveTimestepAndPlot.h"



peano::MappingSpecification   peanoclaw::adapters::SolveTimestepAndPlot::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMostGeneralSpecification()
   & peanoclaw::mappings::Remesh::touchVertexLastTimeSpecification()
   & peanoclaw::mappings::SolveTimestep::touchVertexLastTimeSpecification()
   & peanoclaw::mappings::Plot::touchVertexLastTimeSpecification()


;
}


peano::MappingSpecification   peanoclaw::adapters::SolveTimestepAndPlot::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification::getMostGeneralSpecification()
   & peanoclaw::mappings::Remesh::touchVertexFirstTimeSpecification()
   & peanoclaw::mappings::SolveTimestep::touchVertexFirstTimeSpecification()
   & peanoclaw::mappings::Plot::touchVertexFirstTimeSpecification()


;
}


peano::MappingSpecification   peanoclaw::adapters::SolveTimestepAndPlot::enterCellSpecification() {
  return peano::MappingSpecification::getMostGeneralSpecification()
   & peanoclaw::mappings::Remesh::enterCellSpecification()
   & peanoclaw::mappings::SolveTimestep::enterCellSpecification()
   & peanoclaw::mappings::Plot::enterCellSpecification()


;
}


peano::MappingSpecification   peanoclaw::adapters::SolveTimestepAndPlot::leaveCellSpecification() {
  return peano::MappingSpecification::getMostGeneralSpecification()
   & peanoclaw::mappings::Remesh::leaveCellSpecification()
   & peanoclaw::mappings::SolveTimestep::leaveCellSpecification()
   & peanoclaw::mappings::Plot::leaveCellSpecification()


;
}


peano::MappingSpecification   peanoclaw::adapters::SolveTimestepAndPlot::ascendSpecification() {
  return peano::MappingSpecification::getMostGeneralSpecification()
   & peanoclaw::mappings::Remesh::ascendSpecification()
   & peanoclaw::mappings::SolveTimestep::ascendSpecification()
   & peanoclaw::mappings::Plot::ascendSpecification()


;
}


peano::MappingSpecification   peanoclaw::adapters::SolveTimestepAndPlot::descendSpecification() {
  return peano::MappingSpecification::getMostGeneralSpecification()
   & peanoclaw::mappings::Remesh::descendSpecification()
   & peanoclaw::mappings::SolveTimestep::descendSpecification()
   & peanoclaw::mappings::Plot::descendSpecification()


;
}


peanoclaw::adapters::SolveTimestepAndPlot::SolveTimestepAndPlot() {
}


peanoclaw::adapters::SolveTimestepAndPlot::~SolveTimestepAndPlot() {
}


#if defined(SharedMemoryParallelisation)
peanoclaw::adapters::SolveTimestepAndPlot::SolveTimestepAndPlot(const SolveTimestepAndPlot&  masterThread):
  _map2Remesh(masterThread._map2Remesh) , 
  _map2SolveTimestep(masterThread._map2SolveTimestep) , 
  _map2Plot(masterThread._map2Plot) 

 

{
}


void peanoclaw::adapters::SolveTimestepAndPlot::mergeWithWorkerThread(const SolveTimestepAndPlot& workerThread) {

  _map2Remesh.mergeWithWorkerThread(workerThread._map2Remesh);
  _map2SolveTimestep.mergeWithWorkerThread(workerThread._map2SolveTimestep);
  _map2Plot.mergeWithWorkerThread(workerThread._map2Plot);


}
#endif


void peanoclaw::adapters::SolveTimestepAndPlot::createHangingVertex(
      peanoclaw::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      peanoclaw::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      peanoclaw::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {

  _map2Remesh.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2SolveTimestep.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Plot.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void peanoclaw::adapters::SolveTimestepAndPlot::destroyHangingVertex(
      const peanoclaw::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2Remesh.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2SolveTimestep.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Plot.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void peanoclaw::adapters::SolveTimestepAndPlot::createInnerVertex(
      peanoclaw::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {

  _map2Remesh.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2SolveTimestep.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Plot.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void peanoclaw::adapters::SolveTimestepAndPlot::createBoundaryVertex(
      peanoclaw::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {

  _map2Remesh.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2SolveTimestep.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Plot.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void peanoclaw::adapters::SolveTimestepAndPlot::destroyVertex(
      const peanoclaw::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2Remesh.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2SolveTimestep.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Plot.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void peanoclaw::adapters::SolveTimestepAndPlot::createCell(
      peanoclaw::Cell&                 fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {

  _map2Remesh.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2SolveTimestep.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Plot.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


}


void peanoclaw::adapters::SolveTimestepAndPlot::destroyCell(
      const peanoclaw::Cell&           fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {

  _map2Remesh.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2SolveTimestep.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Plot.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


}

#ifdef Parallel
void peanoclaw::adapters::SolveTimestepAndPlot::mergeWithNeighbour(
  peanoclaw::Vertex&  vertex,
  const peanoclaw::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {

   _map2Remesh.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2SolveTimestep.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2Plot.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );


}

void peanoclaw::adapters::SolveTimestepAndPlot::prepareSendToNeighbour(
  peanoclaw::Vertex&  vertex,
  int  toRank,
  int  level
) {

   _map2Remesh.prepareSendToNeighbour( vertex, toRank, level );
   _map2SolveTimestep.prepareSendToNeighbour( vertex, toRank, level );
   _map2Plot.prepareSendToNeighbour( vertex, toRank, level );


}

void peanoclaw::adapters::SolveTimestepAndPlot::prepareCopyToRemoteNode(
  peanoclaw::Vertex&  localVertex,
  int  toRank
) {

   _map2Remesh.prepareCopyToRemoteNode( localVertex, toRank );
   _map2SolveTimestep.prepareCopyToRemoteNode( localVertex, toRank );
   _map2Plot.prepareCopyToRemoteNode( localVertex, toRank );


}

void peanoclaw::adapters::SolveTimestepAndPlot::prepareCopyToRemoteNode(
  peanoclaw::Cell&  localCell,
  int  toRank
) {

   _map2Remesh.prepareCopyToRemoteNode( localCell, toRank );
   _map2SolveTimestep.prepareCopyToRemoteNode( localCell, toRank );
   _map2Plot.prepareCopyToRemoteNode( localCell, toRank );


}

void peanoclaw::adapters::SolveTimestepAndPlot::mergeWithRemoteDataDueToForkOrJoin(
  peanoclaw::Vertex&  localVertex,
  const peanoclaw::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {

   _map2Remesh.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2SolveTimestep.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2Plot.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );


}

void peanoclaw::adapters::SolveTimestepAndPlot::mergeWithRemoteDataDueToForkOrJoin(
  peanoclaw::Cell&  localCell,
  const peanoclaw::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {

   _map2Remesh.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2SolveTimestep.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2Plot.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );


}

void peanoclaw::adapters::SolveTimestepAndPlot::prepareSendToWorker(
  peanoclaw::Cell&                 fineGridCell,
  peanoclaw::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  peanoclaw::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  peanoclaw::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {

   _map2Remesh.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   _map2SolveTimestep.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   _map2Plot.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );


}

void peanoclaw::adapters::SolveTimestepAndPlot::prepareSendToMaster(
  peanoclaw::Cell&     localCell,
  peanoclaw::Vertex *  vertices,
  const peano::grid::VertexEnumerator&  verticesEnumerator
) {

   _map2Remesh.prepareSendToMaster( localCell, vertices, verticesEnumerator );
   _map2SolveTimestep.prepareSendToMaster( localCell, vertices, verticesEnumerator );
   _map2Plot.prepareSendToMaster( localCell, vertices, verticesEnumerator );


}

void peanoclaw::adapters::SolveTimestepAndPlot::mergeWithMaster(
  const peanoclaw::Cell&           workerGridCell,
  peanoclaw::Vertex * const        workerGridVertices,
  const peano::grid::VertexEnumerator& workerEnumerator,
  peanoclaw::Cell&                 fineGridCell,
  peanoclaw::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  peanoclaw::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  peanoclaw::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker,
    const peanoclaw::State&          workerState,
  peanoclaw::State&                masterState
) {

   _map2Remesh.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2SolveTimestep.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2Plot.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );


}

void peanoclaw::adapters::SolveTimestepAndPlot::receiveDataFromMaster(
  peanoclaw::Cell&                    receivedCell, 
  peanoclaw::Vertex *                 receivedVertices,
  const peano::grid::VertexEnumerator&    verticesEnumerator
) {

   _map2Remesh.receiveDataFromMaster( receivedCell, receivedVertices, verticesEnumerator );
   _map2SolveTimestep.receiveDataFromMaster( receivedCell, receivedVertices, verticesEnumerator );
   _map2Plot.receiveDataFromMaster( receivedCell, receivedVertices, verticesEnumerator );


}


void peanoclaw::adapters::SolveTimestepAndPlot::mergeWithWorker(
  peanoclaw::Cell&           localCell, 
  const peanoclaw::Cell&     receivedMasterCell
) {

   _map2Remesh.mergeWithWorker( localCell, receivedMasterCell );
   _map2SolveTimestep.mergeWithWorker( localCell, receivedMasterCell );
   _map2Plot.mergeWithWorker( localCell, receivedMasterCell );


}

void peanoclaw::adapters::SolveTimestepAndPlot::mergeWithWorker(
  peanoclaw::Vertex&        localVertex,
  const peanoclaw::Vertex&  receivedMasterVertex
) {

   _map2Remesh.mergeWithWorker( localVertex, receivedMasterVertex );
   _map2SolveTimestep.mergeWithWorker( localVertex, receivedMasterVertex );
   _map2Plot.mergeWithWorker( localVertex, receivedMasterVertex );


}
#endif

void peanoclaw::adapters::SolveTimestepAndPlot::touchVertexFirstTime(
      peanoclaw::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {

  _map2Remesh.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2SolveTimestep.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Plot.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void peanoclaw::adapters::SolveTimestepAndPlot::touchVertexLastTime(
      peanoclaw::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2Remesh.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2SolveTimestep.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Plot.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void peanoclaw::adapters::SolveTimestepAndPlot::enterCell(
      peanoclaw::Cell&                 fineGridCell,
      peanoclaw::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      peanoclaw::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      peanoclaw::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {

  _map2Remesh.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2SolveTimestep.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Plot.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


}


void peanoclaw::adapters::SolveTimestepAndPlot::leaveCell(
      peanoclaw::Cell&           fineGridCell,
      peanoclaw::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      peanoclaw::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      peanoclaw::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {

  _map2Remesh.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2SolveTimestep.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Plot.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


}


void peanoclaw::adapters::SolveTimestepAndPlot::beginIteration(
  peanoclaw::State&  solverState
) {

  _map2Remesh.beginIteration( solverState );
  _map2SolveTimestep.beginIteration( solverState );
  _map2Plot.beginIteration( solverState );


}


void peanoclaw::adapters::SolveTimestepAndPlot::endIteration(
  peanoclaw::State&  solverState
) {

  _map2Remesh.endIteration( solverState );
  _map2SolveTimestep.endIteration( solverState );
  _map2Plot.endIteration( solverState );


}




void peanoclaw::adapters::SolveTimestepAndPlot::descend(
  peanoclaw::Cell * const          fineGridCells,
  peanoclaw::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  peanoclaw::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  peanoclaw::Cell&                 coarseGridCell
) {

  _map2Remesh.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2SolveTimestep.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2Plot.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );


}


void peanoclaw::adapters::SolveTimestepAndPlot::ascend(
  peanoclaw::Cell * const    fineGridCells,
  peanoclaw::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  peanoclaw::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  peanoclaw::Cell&           coarseGridCell
) {

  _map2Remesh.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2SolveTimestep.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2Plot.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );


}
