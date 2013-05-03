/*
 * NeighbourCommunicator.cpp
 *
 *  Created on: Mar 15, 2013
 *      Author: unterweg
 */
#include "NeighbourCommunicator.h"

#include "peanoclaw/records/CellDescription.h"
#include "peanoclaw/records/Data.h"

tarch::logging::Log peanoclaw::parallel::NeighbourCommunicator::_log("peanoclaw::parallel::NeighbourCommunicator");

void peanoclaw::parallel::NeighbourCommunicator::sendCellDescription(int cellDescriptionIndex) {
  logTraceInWith1Argument("sendCellDescription", cellDescriptionIndex);
  _cellDescriptionHeap.sendData(cellDescriptionIndex, _remoteRank, _position, _level, peano::heap::NeighbourCommunication);
  logTraceOut("sendCellDescription");
}

void peanoclaw::parallel::NeighbourCommunicator::sendPaddingCellDescription() {
  logTraceIn("sendPaddingCellDescription");
  int cellDescriptionIndex = _cellDescriptionHeap.createData();
  sendCellDescription(cellDescriptionIndex);
  _cellDescriptionHeap.deleteData(cellDescriptionIndex);
  logTraceOut("sendPaddingCellDescription");
}

void peanoclaw::parallel::NeighbourCommunicator::sendDataArray(int index) {
  logTraceInWith3Arguments("sendDataArray", index, _position, _level);
  _dataHeap.sendData(index, _remoteRank, _position, _level, peano::heap::NeighbourCommunication);
  logTraceOut("sendDataArray");
}

void peanoclaw::parallel::NeighbourCommunicator::sendPaddingDataArray() {
  logTraceInWith2Arguments("sendPaddingDataArray", _position, _level);
  int index = _dataHeap.createData();
  sendDataArray(index);
  _dataHeap.deleteData(index);
  logTraceOut("sendPaddingDataArray");
}

void peanoclaw::parallel::NeighbourCommunicator::deleteArraysFromPatch(int cellDescriptionIndex) {
  logTraceInWith1Argument("deleteArraysFromPatch", cellDescriptionIndex);
  if(cellDescriptionIndex != -1) {
    assertion2(_cellDescriptionHeap.isValidIndex(cellDescriptionIndex), _position, _level);
    CellDescription cellDescription = _cellDescriptionHeap.getData(cellDescriptionIndex).at(0);

    if(cellDescription.getUNewIndex() != -1) {
      _dataHeap.deleteData(cellDescription.getUNewIndex());
    }
    if(cellDescription.getUOldIndex() != -1) {
      _dataHeap.deleteData(cellDescription.getUOldIndex());
    }
    if(cellDescription.getAuxIndex() != -1) {
      _dataHeap.deleteData(cellDescription.getAuxIndex());
    }
  }
  logTraceOut("deleteArraysFromPatch");
}

int peanoclaw::parallel::NeighbourCommunicator::receiveDataArray() {
  logTraceIn("receiveDataArray");
  std::vector<Data> remoteArray = _dataHeap.receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);

  int localIndex = _dataHeap.createData();
  std::vector<Data>& localArray = _dataHeap.getData(localIndex);

  // Copy array
  std::vector<Data>::iterator it = remoteArray.begin();
  localArray.assign(it, remoteArray.end());

  logTraceOut("receiveDataArray");
  return localIndex;
}

peanoclaw::parallel::NeighbourCommunicator::NeighbourCommunicator(
  int remoteRank,
  const tarch::la::Vector<DIMENSIONS,double> position,
  int level
) : _remoteRank(remoteRank),
    _position(position),
    _level(level),
    _cellDescriptionHeap(peano::heap::Heap<CellDescription>::getInstance()),
    _dataHeap(peano::heap::Heap<Data>::getInstance()) {
  logTraceInWith3Arguments("NeighbourCommunicator", remoteRank, position, level);

  logTraceOut("NeighbourCommunicator");
}

void peanoclaw::parallel::NeighbourCommunicator::sendPatch(
  int cellDescriptionIndex
) {
  logTraceInWith3Arguments("sendPatch", cellDescriptionIndex, _position, _level);
  if(cellDescriptionIndex != -1) {
    sendCellDescription(cellDescriptionIndex);

    CellDescription cellDescription = _cellDescriptionHeap.getData(cellDescriptionIndex).at(0);

    if(cellDescription.getUNewIndex() != -1) {
      sendDataArray(cellDescription.getUNewIndex());
    } else {
      sendPaddingDataArray();
    }

    if(cellDescription.getUOldIndex() != -1) {
      sendDataArray(cellDescription.getUOldIndex());
    } else {
      sendPaddingDataArray();
    }

    if(cellDescription.getAuxIndex() != -1) {
      sendDataArray(cellDescription.getAuxIndex());
    } else {
      sendPaddingDataArray();
    }
  } else {
    sendPaddingPatch();
  }
  logTraceOut("sendPatch");
}

void peanoclaw::parallel::NeighbourCommunicator::sendPaddingPatch() {
  logTraceInWith2Arguments("sendPaddingPatch", _position, _level);
  sendPaddingCellDescription();
  sendPaddingDataArray(); //UNew
  sendPaddingDataArray(); //UOld
  sendPaddingDataArray(); //Aux
  logTraceOut("sendPaddingPatch");
}

void peanoclaw::parallel::NeighbourCommunicator::receivePatch(int localCellDescriptionIndex) {
  logTraceInWith3Arguments("receivePatch", localCellDescriptionIndex, _position, _level);

  std::vector<CellDescription> remoteCellDescriptionVector = _cellDescriptionHeap.receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);
  assertionEquals2(remoteCellDescriptionVector.size(), 1, _position, _level);
  CellDescription remoteCellDescription = remoteCellDescriptionVector[0];

  //Load arrays and stores according indices in cell description
  if(remoteCellDescription.getAuxIndex() != -1) {
    remoteCellDescription.setAuxIndex(receiveDataArray());
  } else {
    receiveDataArray();
  }
  if(remoteCellDescription.getUOldIndex() != -1) {
    remoteCellDescription.setUOldIndex(receiveDataArray());
  } else {
    receiveDataArray();
  }
  if(remoteCellDescription.getUNewIndex() != -1) {
    remoteCellDescription.setUNewIndex(receiveDataArray());
  } else {
    receiveDataArray();
  }

  //Copy remote cell description to local cell description
  assertion2(localCellDescriptionIndex > 0, _position, _level);
  deleteArraysFromPatch(localCellDescriptionIndex);
  remoteCellDescription.setCellDescriptionIndex(localCellDescriptionIndex);
  _cellDescriptionHeap.getData(localCellDescriptionIndex).at(0) = remoteCellDescription;
  assertionEquals(_cellDescriptionHeap.getData(localCellDescriptionIndex).size(), 1);

  assertionEquals(_cellDescriptionHeap.getData(localCellDescriptionIndex).at(0).getCellDescriptionIndex(), localCellDescriptionIndex);
  logTraceOut("receivePatch");
}

void peanoclaw::parallel::NeighbourCommunicator::receivePaddingPatch() {
  logTraceInWith2Arguments("receivePaddingPatch", _position, _level);
  //Receive padding patch
  _cellDescriptionHeap.receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);

  //Aux
  _dataHeap.receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);

  //UOld
  _dataHeap.receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);

  //UNew
  _dataHeap.receiveData(_remoteRank, _position, _level, peano::heap::NeighbourCommunication);
  logTraceOut("receivePaddingPatch");
}