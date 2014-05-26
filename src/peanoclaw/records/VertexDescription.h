#ifndef _PEANOCLAW_RECORDS_VERTEXDESCRIPTION_H
#define _PEANOCLAW_RECORDS_VERTEXDESCRIPTION_H

#include "peano/utils/Globals.h"
#include "tarch/compiler/CompilerSpecificSettings.h"
#include "peano/utils/PeanoOptimisations.h"
#ifdef Parallel
	#include "tarch/parallel/Node.h"
#endif
#ifdef Parallel
	#include <mpi.h>
#endif
#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"
#include <bitset>
#include <complex>
#include <string>
#include <iostream>

namespace peanoclaw {
   namespace records {
      class VertexDescription;
      class VertexDescriptionPacked;
   }
}

/**
 * @author This class is generated by DaStGen
 * 		   DataStructureGenerator (DaStGen)
 * 		   2007-2009 Wolfgang Eckhardt
 * 		   2012      Tobias Weinzierl
 *
 * 		   build date: 09-02-2014 14:40
 *
 * @date   26/05/2014 09:58
 */
class peanoclaw::records::VertexDescription { 
   
   public:
      
      typedef peanoclaw::records::VertexDescriptionPacked Packed;
      
      enum IterationParity {
         EVEN = 0, ODD = 1
      };
      
      struct PersistentRecords {
         tarch::la::Vector<TWO_POWER_D,int> _indicesOfAdjacentCellDescriptions;
         bool _touched;
         /**
          * Generated
          */
         PersistentRecords();
         
         /**
          * Generated
          */
         PersistentRecords(const tarch::la::Vector<TWO_POWER_D,int>& indicesOfAdjacentCellDescriptions, const bool& touched);
         
         
         inline tarch::la::Vector<TWO_POWER_D,int> getIndicesOfAdjacentCellDescriptions() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _indicesOfAdjacentCellDescriptions;
         }
         
         
         
         inline void setIndicesOfAdjacentCellDescriptions(const tarch::la::Vector<TWO_POWER_D,int>& indicesOfAdjacentCellDescriptions) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _indicesOfAdjacentCellDescriptions = (indicesOfAdjacentCellDescriptions);
         }
         
         
         
         inline bool getTouched() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _touched;
         }
         
         
         
         inline void setTouched(const bool& touched) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _touched = touched;
         }
         
         
         
      };
      
   private: 
      PersistentRecords _persistentRecords;
      
   public:
      /**
       * Generated
       */
      VertexDescription();
      
      /**
       * Generated
       */
      VertexDescription(const PersistentRecords& persistentRecords);
      
      /**
       * Generated
       */
      VertexDescription(const tarch::la::Vector<TWO_POWER_D,int>& indicesOfAdjacentCellDescriptions, const bool& touched);
      
      /**
       * Generated
       */
      ~VertexDescription();
      
      
      inline tarch::la::Vector<TWO_POWER_D,int> getIndicesOfAdjacentCellDescriptions() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         return _persistentRecords._indicesOfAdjacentCellDescriptions;
      }
      
      
      
      inline void setIndicesOfAdjacentCellDescriptions(const tarch::la::Vector<TWO_POWER_D,int>& indicesOfAdjacentCellDescriptions) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         _persistentRecords._indicesOfAdjacentCellDescriptions = (indicesOfAdjacentCellDescriptions);
      }
      
      
      
      inline int getIndicesOfAdjacentCellDescriptions(int elementIndex) const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         assertion(elementIndex>=0);
         assertion(elementIndex<TWO_POWER_D);
         return _persistentRecords._indicesOfAdjacentCellDescriptions[elementIndex];
         
      }
      
      
      
      inline void setIndicesOfAdjacentCellDescriptions(int elementIndex, const int& indicesOfAdjacentCellDescriptions) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         assertion(elementIndex>=0);
         assertion(elementIndex<TWO_POWER_D);
         _persistentRecords._indicesOfAdjacentCellDescriptions[elementIndex]= indicesOfAdjacentCellDescriptions;
         
      }
      
      
      
      inline bool getTouched() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         return _persistentRecords._touched;
      }
      
      
      
      inline void setTouched(const bool& touched) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         _persistentRecords._touched = touched;
      }
      
      
      /**
       * Generated
       */
      static std::string toString(const IterationParity& param);
      
      /**
       * Generated
       */
      static std::string getIterationParityMapping();
      
      /**
       * Generated
       */
      std::string toString() const;
      
      /**
       * Generated
       */
      void toString(std::ostream& out) const;
      
      
      PersistentRecords getPersistentRecords() const;
      /**
       * Generated
       */
      VertexDescriptionPacked convert() const;
      
      
   #ifdef Parallel
      protected:
         static tarch::logging::Log _log;
         
      public:
         
         /**
          * Global that represents the mpi datatype.
          * There are two variants: Datatype identifies only those attributes marked with
          * parallelise. FullDatatype instead identifies the whole record with all fields.
          */
         static MPI_Datatype Datatype;
         static MPI_Datatype FullDatatype;
         
         /**
          * Initializes the data type for the mpi operations. Has to be called
          * before the very first send or receive operation is called.
          */
         static void initDatatype();
         
         static void shutdownDatatype();
         
         void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
         
         void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
         
         static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
         
         #endif
            
         };
         
         /**
          * @author This class is generated by DaStGen
          * 		   DataStructureGenerator (DaStGen)
          * 		   2007-2009 Wolfgang Eckhardt
          * 		   2012      Tobias Weinzierl
          *
          * 		   build date: 09-02-2014 14:40
          *
          * @date   26/05/2014 09:58
          */
         class peanoclaw::records::VertexDescriptionPacked { 
            
            public:
               
               typedef peanoclaw::records::VertexDescription::IterationParity IterationParity;
               
               struct PersistentRecords {
                  tarch::la::Vector<TWO_POWER_D,int> _indicesOfAdjacentCellDescriptions;
                  bool _touched;
                  /**
                   * Generated
                   */
                  PersistentRecords();
                  
                  /**
                   * Generated
                   */
                  PersistentRecords(const tarch::la::Vector<TWO_POWER_D,int>& indicesOfAdjacentCellDescriptions, const bool& touched);
                  
                  
                  inline tarch::la::Vector<TWO_POWER_D,int> getIndicesOfAdjacentCellDescriptions() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     return _indicesOfAdjacentCellDescriptions;
                  }
                  
                  
                  
                  inline void setIndicesOfAdjacentCellDescriptions(const tarch::la::Vector<TWO_POWER_D,int>& indicesOfAdjacentCellDescriptions) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     _indicesOfAdjacentCellDescriptions = (indicesOfAdjacentCellDescriptions);
                  }
                  
                  
                  
                  inline bool getTouched() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     return _touched;
                  }
                  
                  
                  
                  inline void setTouched(const bool& touched) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     _touched = touched;
                  }
                  
                  
                  
               };
               
            private: 
               PersistentRecords _persistentRecords;
               
            public:
               /**
                * Generated
                */
               VertexDescriptionPacked();
               
               /**
                * Generated
                */
               VertexDescriptionPacked(const PersistentRecords& persistentRecords);
               
               /**
                * Generated
                */
               VertexDescriptionPacked(const tarch::la::Vector<TWO_POWER_D,int>& indicesOfAdjacentCellDescriptions, const bool& touched);
               
               /**
                * Generated
                */
               ~VertexDescriptionPacked();
               
               
               inline tarch::la::Vector<TWO_POWER_D,int> getIndicesOfAdjacentCellDescriptions() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  return _persistentRecords._indicesOfAdjacentCellDescriptions;
               }
               
               
               
               inline void setIndicesOfAdjacentCellDescriptions(const tarch::la::Vector<TWO_POWER_D,int>& indicesOfAdjacentCellDescriptions) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  _persistentRecords._indicesOfAdjacentCellDescriptions = (indicesOfAdjacentCellDescriptions);
               }
               
               
               
               inline int getIndicesOfAdjacentCellDescriptions(int elementIndex) const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  assertion(elementIndex>=0);
                  assertion(elementIndex<TWO_POWER_D);
                  return _persistentRecords._indicesOfAdjacentCellDescriptions[elementIndex];
                  
               }
               
               
               
               inline void setIndicesOfAdjacentCellDescriptions(int elementIndex, const int& indicesOfAdjacentCellDescriptions) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  assertion(elementIndex>=0);
                  assertion(elementIndex<TWO_POWER_D);
                  _persistentRecords._indicesOfAdjacentCellDescriptions[elementIndex]= indicesOfAdjacentCellDescriptions;
                  
               }
               
               
               
               inline bool getTouched() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  return _persistentRecords._touched;
               }
               
               
               
               inline void setTouched(const bool& touched) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  _persistentRecords._touched = touched;
               }
               
               
               /**
                * Generated
                */
               static std::string toString(const IterationParity& param);
               
               /**
                * Generated
                */
               static std::string getIterationParityMapping();
               
               /**
                * Generated
                */
               std::string toString() const;
               
               /**
                * Generated
                */
               void toString(std::ostream& out) const;
               
               
               PersistentRecords getPersistentRecords() const;
               /**
                * Generated
                */
               VertexDescription convert() const;
               
               
            #ifdef Parallel
               protected:
                  static tarch::logging::Log _log;
                  
               public:
                  
                  /**
                   * Global that represents the mpi datatype.
                   * There are two variants: Datatype identifies only those attributes marked with
                   * parallelise. FullDatatype instead identifies the whole record with all fields.
                   */
                  static MPI_Datatype Datatype;
                  static MPI_Datatype FullDatatype;
                  
                  /**
                   * Initializes the data type for the mpi operations. Has to be called
                   * before the very first send or receive operation is called.
                   */
                  static void initDatatype();
                  
                  static void shutdownDatatype();
                  
                  void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
                  
                  void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
                  
                  static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
                  
                  #endif
                     
                  };
                  
                  #endif
                  
