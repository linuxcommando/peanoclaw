#ifndef _PEANOCLAW_RECORDS_REPOSITORYSTATE_H
#define _PEANOCLAW_RECORDS_REPOSITORYSTATE_H

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
      class RepositoryState;
      class RepositoryStatePacked;
   }
}

/**
 * @author This class is generated by DaStGen
 * 		   DataStructureGenerator (DaStGen)
 * 		   2007-2009 Wolfgang Eckhardt
 * 		   2012      Tobias Weinzierl
 *
 * 		   build date: 12-04-2013 09:18
 *
 * @date   29/04/2013 12:57
 */
class peanoclaw::records::RepositoryState { 
   
   public:
      
      typedef peanoclaw::records::RepositoryStatePacked Packed;
      
      enum Action {
         WriteCheckpoint = 0, ReadCheckpoint = 1, Terminate = 2, UseAdapterInitialiseGrid = 3, UseAdapterPlot = 4, UseAdapterRemesh = 5, UseAdapterSolveTimestep = 6, UseAdapterSolveTimestepAndPlot = 7, UseAdapterGatherCurrentSolution = 8, UseAdapterCleanup = 9
      };
      
      struct PersistentRecords {
         Action _action;
         bool _reduceState;
         /**
          * Generated
          */
         PersistentRecords();
         
         /**
          * Generated
          */
         PersistentRecords(const Action& action, const bool& reduceState);
         
         
         inline Action getAction() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _action;
         }
         
         
         
         inline void setAction(const Action& action) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _action = action;
         }
         
         
         
         inline bool getReduceState() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _reduceState;
         }
         
         
         
         inline void setReduceState(const bool& reduceState) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _reduceState = reduceState;
         }
         
         
         
      };
      
   private: 
      PersistentRecords _persistentRecords;
      
   public:
      /**
       * Generated
       */
      RepositoryState();
      
      /**
       * Generated
       */
      RepositoryState(const PersistentRecords& persistentRecords);
      
      /**
       * Generated
       */
      RepositoryState(const Action& action, const bool& reduceState);
      
      /**
       * Generated
       */
      virtual ~RepositoryState();
      
      
      inline Action getAction() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         return _persistentRecords._action;
      }
      
      
      
      inline void setAction(const Action& action) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         _persistentRecords._action = action;
      }
      
      
      
      inline bool getReduceState() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         return _persistentRecords._reduceState;
      }
      
      
      
      inline void setReduceState(const bool& reduceState) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         _persistentRecords._reduceState = reduceState;
      }
      
      
      /**
       * Generated
       */
      static std::string toString(const Action& param);
      
      /**
       * Generated
       */
      static std::string getActionMapping();
      
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
      RepositoryStatePacked convert() const;
      
      
   #ifdef Parallel
      protected:
         static tarch::logging::Log _log;
         
         int _senderDestinationRank;
         
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
         
         void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
         
         void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
         
         static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
         
         int getSenderRank() const;
         
   #endif
      
   };
   
   #ifndef DaStGenPackedPadding
     #define DaStGenPackedPadding 1      // 32 bit version
     // #define DaStGenPackedPadding 2   // 64 bit version
   #endif
   
   
   #ifdef PackedRecords
      #pragma pack (push, DaStGenPackedPadding)
   #endif
   
   /**
    * @author This class is generated by DaStGen
    * 		   DataStructureGenerator (DaStGen)
    * 		   2007-2009 Wolfgang Eckhardt
    * 		   2012      Tobias Weinzierl
    *
    * 		   build date: 12-04-2013 09:18
    *
    * @date   29/04/2013 12:57
    */
   class peanoclaw::records::RepositoryStatePacked { 
      
      public:
         
         typedef peanoclaw::records::RepositoryState::Action Action;
         
         struct PersistentRecords {
            Action _action;
            bool _reduceState;
            /**
             * Generated
             */
            PersistentRecords();
            
            /**
             * Generated
             */
            PersistentRecords(const Action& action, const bool& reduceState);
            
            
            inline Action getAction() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
               return _action;
            }
            
            
            
            inline void setAction(const Action& action) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
               _action = action;
            }
            
            
            
            inline bool getReduceState() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
               return _reduceState;
            }
            
            
            
            inline void setReduceState(const bool& reduceState) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
               _reduceState = reduceState;
            }
            
            
            
         };
         
      private: 
         PersistentRecords _persistentRecords;
         
      public:
         /**
          * Generated
          */
         RepositoryStatePacked();
         
         /**
          * Generated
          */
         RepositoryStatePacked(const PersistentRecords& persistentRecords);
         
         /**
          * Generated
          */
         RepositoryStatePacked(const Action& action, const bool& reduceState);
         
         /**
          * Generated
          */
         virtual ~RepositoryStatePacked();
         
         
         inline Action getAction() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _persistentRecords._action;
         }
         
         
         
         inline void setAction(const Action& action) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _persistentRecords._action = action;
         }
         
         
         
         inline bool getReduceState() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _persistentRecords._reduceState;
         }
         
         
         
         inline void setReduceState(const bool& reduceState) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _persistentRecords._reduceState = reduceState;
         }
         
         
         /**
          * Generated
          */
         static std::string toString(const Action& param);
         
         /**
          * Generated
          */
         static std::string getActionMapping();
         
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
         RepositoryState convert() const;
         
         
      #ifdef Parallel
         protected:
            static tarch::logging::Log _log;
            
            int _senderDestinationRank;
            
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
            
            void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
            
            void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
            
            static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
            
            int getSenderRank() const;
            
      #endif
         
      };
      
      #ifdef PackedRecords
      #pragma pack (pop)
      #endif
      
      
      #endif
      
