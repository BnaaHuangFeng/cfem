# pragma once
#include "MeshSystem/MeshSystem.h"
#include <vector>
#include <map>
using namespace std;
enum class DirichletMethod{
    SETUNIT,
    SETLARGE
};
class BCsSystem
{
protected:
    bool m_ifSetMeshSysPtr; 
    bool m_ifHasReadBCDes;                            /**< if has read the boundary condition description*/
    bool m_ifHasInit;                                 /**< if has completed inition*/
    Vec m_max_entry_vec;                              /**< Vec storing the maximum of every row of jacobian*/
    MeshSystem *m_meshSysPtr;                         /**< ptr to mesh system it relies on*/
    BCDescription *m_bcDesPtr;                        /**< ptr to boundary condition description*/
    PetscInt m_mConstrainedDof;                       /**< num of constrained dofs*/
    map<PetscInt,PetscScalar> m_bcValsMap;            /**< (constrained dof id in this rank, prescribed dof val)*/
    PetscInt *m_arrayConstrainedRows;                 /**< ptr to dynamic array storing constrained dof id in this rank (for acceleration)*/
    const PetscScalar m_penalty_coef=1.0e10;
    PetscScalar *m_arrayPresetVals;                   /**< ptr to dynamic array storing prescribed dof val in this rank (for acceleration)*/
    PetscScalar *m_arrayZero;                         /**< ptr to dynamic array of the same size arrayConstrainedRows as storing zero val in this rank (for acceleration)*/
    PetscScalar  m_penalty;                           /**< penalty for apply penalty pivot*/ 
protected:
    inline void setMeshSysPtr(MeshSystem *t_meshSysPtr){
        m_meshSysPtr=t_meshSysPtr;
        m_ifSetMeshSysPtr=true;
    }    
public:
    BCsSystem():
        m_ifSetMeshSysPtr(false),m_ifHasReadBCDes(false),m_ifHasInit(false),
        m_meshSysPtr(nullptr),m_bcDesPtr(nullptr),
        m_arrayConstrainedRows(nullptr),m_arrayPresetVals(nullptr),m_arrayZero(nullptr),
        m_penalty(0.0){};

    BCsSystem(BCDescription *t_bcDesPtr):
        m_ifSetMeshSysPtr(false),m_ifHasReadBCDes(true),m_ifHasInit(false),
        m_meshSysPtr(nullptr),m_bcDesPtr(t_bcDesPtr),
        m_arrayConstrainedRows(nullptr),m_arrayPresetVals(nullptr),m_arrayZero(nullptr),
        m_penalty(0.0){};

    BCsSystem(BCDescription *t_bcDesPtr, MeshSystem *t_meshSysPtr):
        m_ifHasReadBCDes(true),m_ifHasInit(false),m_meshSysPtr(nullptr),
        m_bcDesPtr(t_bcDesPtr),
        m_arrayConstrainedRows(nullptr),m_arrayPresetVals(nullptr),m_arrayZero(nullptr),
        m_penalty(0.0){
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
     * set solution initial guess, it's a must do not just for efficiency
    */
    virtual PetscErrorCode setInitialSolution(PetscScalar facInc)=0;
    /**
     * apply boundary condition to global residual Vec
     * @param residualPtr > ptr to global residual Vec
    */
    virtual PetscErrorCode applyResidualBoundaryCondition(Vec *residualPtr)=0;
    /**
     * note: this function will zero t_b first
    */
    virtual PetscErrorCode applyBoundaryConditionArc(Mat *AMatrixPtr, Vec *t_b)=0;
    /**
     * apply boundary condition to global jacobian Mat
     * @param AMatrixPtr > ptr to global Jacobian Mat
    */
    virtual PetscErrorCode applyJacobianBoundaryCondition(Mat *AMatrixPtr)=0;
    virtual PetscErrorCode update_penalty(Mat *AMatrixPtr)=0;
public:
    DirichletMethod m_drclt_method=DirichletMethod::SETUNIT;
    Vec m_uIncInitial;                                /**< Vec for initial incremental dof*/ 
};
