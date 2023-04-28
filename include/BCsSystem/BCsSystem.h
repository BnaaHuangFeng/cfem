# pragma once
#include "MeshSystem/MeshSystem.h"
#include <vector>
#include <map>
using namespace std;
class BCsSystem
{
public:
    BCsSystem():m_meshPtr(nullptr),m_bcDesPtr(nullptr),m_ifHasReadBCDes(false),m_ifHasInit(false){};
    BCsSystem(BCDescription *t_bcDesPtr):m_meshPtr(nullptr),m_bcDesPtr(t_bcDesPtr),m_ifHasReadBCDes(true),m_ifHasInit(false){};
    virtual PetscErrorCode init()=0;
    virtual PetscErrorCode init(BCDescription *t_bcDesPtr)=0;
    virtual ~BCsSystem(){};
    /**
     * apply boundary condition to global jacobian Mat, and global residual Vec
     * @param factor > boundary condition scaling factor
     * @param residualPtr > ptr to global residual Vec
     * @param AMatrixPtr > ptr to global Jacobian Mat
    */
    virtual PetscErrorCode applyBoundaryCondition(PetscScalar factor,Vec *residualPtr, Mat *AMatrixPtr)=0;
public:
    MeshSystem *m_meshPtr;                          /**< ptr to mesh system it relies on*/
    BCDescription *m_bcDesPtr;                      /**< ptr to boundary condition description*/
    map<PetscInt,PetscScalar> m_bcValsMap;          /**< (constrained dof id in this rank, prescribed dof val)*/
    bool m_ifHasReadBCDes;                            /**< if has read the boundary condition description*/
    bool m_ifHasInit;                                 /**< if has completed inition*/
};
