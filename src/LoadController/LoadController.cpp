#include "LoadController/LoadController.h"
LoadController::LoadController(StepDescriptiom *t_stepDesPtr):
m_stepDesPtr(t_stepDesPtr),m_mGrowth(5),m_convergeCnt(0),m_iterCnt2(0),m_iterCnt1(0),
m_factor2(0.0),m_factorInc2(0.0),m_factor1(0.0),m_factorInc1(0.0){
    m_growthRatio=t_stepDesPtr->s_growFactor;
    m_cutbackRatio=t_stepDesPtr->s_cutbackFactor;
    m_factorIncMin=t_stepDesPtr->s_dtmin;
    m_factorIncMax=t_stepDesPtr->s_dtmax;
    m_factorInc2=t_stepDesPtr->s_dt0;

}
bool LoadController::update(bool ifConverged, int mIter){
    if(ifConverged){
        m_iterCnt2=m_iterCnt1;
        m_iterCnt1=mIter;
        m_factorInc2=m_factorInc1;
        m_factor2=m_factor1;
        if(m_factor2==m_factorFinal)return true;
        if(m_convergeCnt>=m_mGrowth){
            m_factorInc1=min(m_growthRatio*m_factorInc2,m_factorIncMax);
            m_factorInc1=min(m_factorInc1,m_factorFinal-m_factor2);
            m_convergeCnt=0;
        }
        else{
            m_factorInc1=min(m_factorInc2,m_factorFinal-m_factor2);
            ++m_convergeCnt;
        }
        m_factor1=m_factor2+m_factorInc1;
    }
    else{
        m_factorInc1=max(m_factorInc1*m_cutbackRatio,m_factorIncMin);
        m_factorInc1=min(m_factorInc1,m_factorFinal-m_factor2);
        m_factor1=m_factor2+m_factorInc1;
    }
    return false;
}
PetscErrorCode LoadController::applyLoad(PetscScalar factor,Vec *residualPtr){
    PetscCall(VecAXPY(*residualPtr,-factor,*m_loadVecPtr));
    return 0;
}