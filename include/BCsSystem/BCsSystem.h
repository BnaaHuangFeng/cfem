# pragma once
#include "MeshSystem/MeshSystem.h"
#include <vector>
#include <map>
using namespace std;
class BCsSystem
{
protected:
    bool m_ifSetMeshSysPtr; 
    bool m_ifHasReadBCDes;                            /**< if has read the boundary condition description*/
    bool m_ifHasInit;                                 /**< if has completed inition*/
    MeshSystem *m_meshSysPtr;                         /**< ptr to mesh system it relies on*/
    BCDescription *m_bcDesPtr;                        /**< ptr to boundary condition description*/
    PetscInt m_mConstrainedDof;                       /**< num of constrained dofs*/
    map<PetscInt,PetscScalar> m_bcValsMap;            /**< (constrained dof id in this rank, prescribed dof val)*/
    PetscInt *m_arrayConstrainedRows;                 /**< ptr to dynamic array storing constrained dof id in this rank (for acceleration)*/
    PetscScalar *m_arrayPresetVals;                   /**< ptr to dynamic array storing prescribed dof val in this rank (for acceleration)*/
    PetscScalar *m_arrayZero;                         /**< ptr to dynamic array of the same size arrayConstrainedRows as storing zero val in this rank (for acceleration)*/
protected:
    inline void setMeshSysPtr(MeshSystem *t_meshSysPtr){
        m_meshSysPtr=t_meshSysPtr;
        m_ifSetMeshSysPtr=true;
    }    
public:
    BCsSystem():
        m_ifSetMeshSysPtr(false),m_ifHasReadBCDes(false),m_ifHasInit(false),
        m_meshSysPtr(nullptr),m_bcDesPtr(nullptr),
        m_arrayConstrainedRows(nullptr),m_arrayPresetVals(nullptr),m_arrayZero(nullptr){};

    BCsSystem(BCDescription *t_bcDesPtr):
        m_ifSetMeshSysPtr(false),m_ifHasReadBCDes(true),m_ifHasInit(false),
        m_meshSysPtr(nullptr),m_bcDesPtr(t_bcDesPtr),
        m_arrayConstrainedRows(nullptr),m_arrayPresetVals(nullptr),m_arrayZero(nullptr){};

    BCsSystem(BCDescription *t_bcDesPtr, MeshSystem *t_meshSysPtr):
        m_ifHasReadBCDes(true),m_ifHasInit(false),m_meshSysPtr(nullptr),
        m_bcDesPtr(t_bcDesPtr),
        m_arrayConstrainedRows(nullptr),m_arrayPresetVals(nullptr),m_arrayZero(nullptr){
        setMeshSysPtr(t_meshSysPtr);
    };
    virtual PetscErrorCode init()=0;
    virtual PetscErrorCode init(BCDescription *t_bcDesPtr)=0;
    virtual PetscErrorCode init(BCDescription *t_bcDesPtr, MeshSystem *t_meshSysPtr)=0;
    /**
     * check solution system's inition
    */
    virtual PetscErrorCode checkInit()=0;
    virtual ~BCsSystem(){};
    /**
     * set arrayPresetVals
     * @param facInc > Incremental boundary condition scaling factor
    */
    virtual PetscErrorCode setArrayPresetVals(PetscScalar facInc)=0;
    /**
     * set solution initial guess, it's a must do not just for efficiency
     * @param uInc1Ptr > ptr to Vec to store the solution guess
    */
    virtual PetscErrorCode setInitialSolution()=0;
    /**
     * apply boundary condition to global residual Vec
     * @param residualPtr > ptr to global residual Vec
    */
    virtual PetscErrorCode applyResidualBoundaryCondition(Vec *residualPtr)=0;
    /**
     * apply boundary condition to global jacobian Mat
     * @param AMatrixPtr > ptr to global Jacobian Mat
    */
    virtual PetscErrorCode applyJacobianBoundaryCondition(Mat *AMatrixPtr)=0;
public:
    Vec m_uIncInitial;                                /**< Vec for initial incremental dof*/
};
