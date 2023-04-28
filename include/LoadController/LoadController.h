#pragma once
#include "petsc.h"
#include "MeshSystem/MeshSystem.h" 
#include "InputSystem/DescriptionInfo.h"
class LoadController
{
private:
    /* data */
public:
    LoadController():m_stepDesPtr(nullptr){};
    LoadController(StepDescriptiom *t_stepDesPtr);
    ~LoadController(){};
private:
    StepDescriptiom *m_stepDesPtr;  /**< ptr to step description*/
    double m_factorFinal;   /**< final factor to end simulation*/
    double m_growthRatio;   /**< factor growth ratio*/
    double m_cutbackRatio;  /**< factor cut back ratio*/
    double m_factorIncMin;  /**< min incremental factor*/
    double m_factorIncMax;  /**< max incremental factor*/
public:
    int m_mGrowth;          /**< factor will grow if there are m_mGrowth consecutive convergence*/
    int m_convergeCnt;      /**< consecutive convergence count*/
    int m_iterCnt2;         /**< iteration count of latest convergence*/
    int m_iterCnt1;         /**< iteration count of latest iteration*/
    double m_factor2;       /**< total factor of latest convergence*/
    double m_factorInc2;    /**< incremental factor of latest convergence*/
    double m_factor1;       /**< current total factor*/
    double m_factorInc1;    /**< current incremental factor*/
    double m_l0;            /**< arc length of latest convergence*/
    MeshSystem *m_meshPtr;  /**< ptr to its relied mesh system*/
    Vec *m_loadVecPtr;      /**< global load Vec ptr*/
public:
    /**
     * apply outer load to global residual Vec
     * @param factor > total load scaling factor
     * @param residualPtr > ptr to global residual Vec
    */
    PetscErrorCode applyLoad(PetscScalar factor,Vec *residualPtr);
    /**
     * update load controller's memeber, including factor, count, ...
     * @param ifConverged > if the latest increment converged.
     * @param mIter > iteration num of the latest increment.
     * @param return < if the final factor reached.
    */
    bool update(bool ifConverged, int mIter);
};