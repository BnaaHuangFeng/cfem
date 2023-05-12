#pragma once
#include "petsc.h"
#include "MeshSystem/MeshSystem.h" 
#include "InputSystem/DescriptionInfo.h"
class LoadController
{
private:
    bool m_ifReadStepDes;   /**< if has read step desctiption*/
    bool m_ifSetMeshSys;    /**< if has set mesh system ptr*/
    bool m_ifHasNodeLoad;     /**< if appliy outer load, set by load controller*/
    bool m_ifInitialArcLenSet;
    StepDescriptiom *m_stepDesPtr;  /**< ptr to step description*/
    MeshSystem *m_meshPtr;  /**< ptr to its relied mesh system*/
    double m_factorFinal;   /**< final factor to end simulation*/
    double m_growthRatio;   /**< factor growth ratio*/
    double m_cutbackRatio;  /**< factor cut back ratio*/
    double m_factorIncMin;  /**< min incremental factor*/
    double m_factorIncMax;  /**< max incremental factor*/
    double m_expIters;      /**< expected iterations num (for arc length method)*/
    double m_arcLenMaxParam;/**< arc length max paramater*/
    double m_arcLenMax;     /**< max arc length limit*/
    double m_arcLen1;       /**< arc len of current incre*/
    double m_arcLen2;       /**< arc len of last converged incre*/
public:
    int m_mGrowth;          /**< factor will grow if there are m_mGrowth consecutive convergence*/
    int m_convergeCnt;      /**< consecutive convergence count*/
    int m_iterExpCnt;       /**< expected iteration num*/   
    double m_factor2;       /**< total factor of latest convergence*/
    double m_factorInc2;    /**< incremental factor of latest convergence*/
    double m_factor1;       /**< current total factor*/
    double m_factorInc1;    /**< current incremental factor*/
    Vec *m_loadVecPtr;      /**< global load Vec ptr*/
private:
    /**
     * get the arc len for increment except the 1st increment
     * @param t_lastArcLen > arc len of the last converged increment
     * @param t_mIters > iteration num of the last converged increment
     * @return < the arc len of the current increment (except the 1st increment)
    */
    double nextArcLen(int t_mIters);
public:
    LoadController();
    LoadController(StepDescriptiom *t_stepDesPtr);
    LoadController(StepDescriptiom *t_stepDesPtr, MeshSystem *t_meshSysPtr);
    void readStepDes(StepDescriptiom *t_stepDesPtr);
    void setMeshSysPtr(MeshSystem *t_meshSysPtr);
    /**
     * check if the load controller's init completed
    */
    PetscErrorCode checkInit();
    ~LoadController(){};
    /**
     * apply outer load to global residual Vec
     * @param factor > total load scaling factor
     * @param residualPtr > ptr to global residual Vec
    */
    PetscErrorCode applyLoad(PetscScalar factor,Vec *residualPtr);
    /**
     * update load controller's memeber, including factor, count, ...
     * @param ifConverged > if the latest increment converged.
     * @param return < if the final factor reached.
    */
    bool update(bool ifConverged);
    /**
     * for arc length method, update load controller's memeber, including factor, count, ...
     * @param ifConverged > if the latest increment converged.
     * @param factorInc > increnmental factor get by arc length method
     * @param miters > num of inter of last convergence
     * @param return < if the final factor reached.
    */
    bool update(bool ifConverged, double factorInc, int miters);
    inline double factorFinal(){return m_factorFinal;}
    inline void setInitialArc(double initialArcLen){
        m_arcLen1=initialArcLen;
        m_arcLenMax=initialArcLen*m_arcLenMaxParam;
        m_ifInitialArcLenSet=true;
    };
    inline Vec *getLoadVecPtr(){return m_loadVecPtr;};
    inline PetscScalar getArcLen(){return m_arcLen1;};
};