# pragma once
#include "BCsSystem/BCsSystem.h"
using namespace std;
class BCsSysStructured2d: public BCsSystem
{
private:
    Vec m_u;    /**< global Vec for apply bcs*/
    Vec m_b;    /**< global Vec for apply bcs*/
public:
    BCsSysStructured2d(){};
    BCsSysStructured2d(BCDescription *t_bcDesPtr):BCsSystem(t_bcDesPtr){};
    virtual PetscErrorCode init();
    virtual PetscErrorCode init(BCDescription *t_bcDesPtr);
    virtual ~BCsSysStructured2d();
    /**
     * apply boundary condition to global jacobian Mat, and global residual Vec
     * @param factor > boundary condition scaling factor
     * @param residualPtr > ptr to global residual Vec
     * @param AMatrixPtr > ptr to global Jacobian Mat
    */
    virtual PetscErrorCode applyBoundaryCondition(PetscScalar factor,Vec *residualPtr, Mat *AMatrixPtr);
};
