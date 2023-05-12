# pragma once
#include "BCsSystem/BCsSystem.h"
using namespace std;
class BCsSysStructured2d: public BCsSystem
{
private:
    /**
     * set arrayPresetVals
     * @param facInc > Incremental boundary condition scaling factor
    */
    PetscErrorCode setArrayPresetVals(PetscScalar facInc);
public:
    BCsSysStructured2d(){};
    BCsSysStructured2d(BCDescription *t_bcDesPtr):BCsSystem(t_bcDesPtr){};
    BCsSysStructured2d(BCDescription *t_bcDesPtr, MeshSystem *t_meshSysPtr):BCsSystem(t_bcDesPtr, t_meshSysPtr){};
    virtual PetscErrorCode init();
    virtual PetscErrorCode init(BCDescription *t_bcDesPtr);
    virtual PetscErrorCode init(BCDescription *t_bcDesPtr, MeshSystem *t_meshSysPtr);
    virtual PetscErrorCode checkInit();
    virtual ~BCsSysStructured2d();

    /**
     * set solution initial guess, it's a must do not just for efficiency
     * @param uInc1Ptr > ptr to Vec to store the solution guess
    */
    virtual PetscErrorCode setInitialSolution(PetscScalar facInc);
    /**
     * apply boundary condition to global residual Vec
     * @param residualPtr > ptr to global residual Vec
    */
    virtual PetscErrorCode applyResidualBoundaryCondition(Vec *residualPtr);
    virtual PetscErrorCode applyBoundaryConditionArc(Mat *AMatrixPtr, Vec *t_b);
    /**
     * apply boundary condition to global jacobian Mat
     * @param AMatrixPtr > ptr to global Jacobian Mat
    */
    virtual PetscErrorCode applyJacobianBoundaryCondition(Mat *AMatrixPtr);
    virtual PetscErrorCode update_penalty(Mat *AMatrixPtr);
};
