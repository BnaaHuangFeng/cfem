#include "LoadController/LoadController.h"
LoadController::LoadController():
m_stepDesPtr(nullptr),m_meshPtr(nullptr),
m_factorFinal(1.0),m_growthRatio(1.1),m_cutbackRatio(0.25),m_factorIncMin(0.001),m_factorIncMax(1.0),
m_mGrowth(5),m_convergeCnt(0),
m_factor2(0.0),m_factorInc2(0.0),m_factor1(0.0),m_factorInc1(0.0){
    m_arcLen1=0.0;
    m_arcLen2=0.0;
    m_ifReadStepDes=false;
    m_ifSetMeshSys=false;
    m_ifHasNodeLoad=false;
    m_ifInitialArcLenSet=false;
}
LoadController::LoadController(StepDescriptiom *t_stepDesPtr):
m_meshPtr(nullptr),
m_mGrowth(5),m_convergeCnt(0),
m_factor2(0.0),m_factorInc2(0.0),m_factor1(0.0),m_factorInc1(0.0){
    m_arcLen1=0.0;
    m_arcLen2=0.0;
    m_ifReadStepDes=false;
    readStepDes(t_stepDesPtr);
    m_ifSetMeshSys=false;
    m_ifHasNodeLoad=false;
    m_ifInitialArcLenSet=false;
}
LoadController::LoadController(StepDescriptiom *t_stepDesPtr, MeshSystem *t_meshSysPtr):
m_mGrowth(5),m_convergeCnt(0),
m_factor2(0.0),m_factorInc2(0.0),m_factor1(0.0),m_factorInc1(0.0){
    m_arcLen1=0.0;
    m_arcLen2=0.0;
    m_ifReadStepDes=false;
    m_ifSetMeshSys=false;
    readStepDes(t_stepDesPtr);
    setMeshSysPtr(t_meshSysPtr);
    m_ifHasNodeLoad=false;
    m_ifInitialArcLenSet=false;
}
bool LoadController::update(bool ifConverged){
    if(ifConverged){
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
bool LoadController::update(bool ifConverged, double factorInc, int miters){
    if(ifConverged){
        m_factorInc2=m_factorInc1;
        m_factorInc1=factorInc;
        m_factor2=m_factor1;
        m_factor1=m_factor1+factorInc;
        m_arcLen2=m_arcLen1;
        m_arcLen1=nextArcLen(miters);
        if(m_factor1>=m_factorFinal)return true;
    }
    else{
        m_arcLen1=m_arcLen1*m_cutbackRatio;
    }
    return false;    
}
PetscErrorCode LoadController::checkInit(){
    if(!m_ifReadStepDes){
        MessagePrinter::printErrorTxt("LoadController: it has not read step description");
        MessagePrinter::exitcfem();
    }
    if(!m_ifSetMeshSys){
        MessagePrinter::printErrorTxt("LoadController: its relied mesh system has not yet been set.");
        MessagePrinter::exitcfem();        
    }
    return 0;
}
PetscErrorCode LoadController::applyLoad(PetscScalar factor,Vec *residualPtr){
    if(!m_ifHasNodeLoad)return 0;
    PetscCall(VecAXPY(*residualPtr,-factor,*m_loadVecPtr));
    return 0;
}
void LoadController::readStepDes(StepDescriptiom *t_stepDesPtr){
    m_stepDesPtr=t_stepDesPtr;
    m_growthRatio=t_stepDesPtr->s_growFactor;
    m_cutbackRatio=t_stepDesPtr->s_cutbackFactor;
    m_factorIncMin=t_stepDesPtr->s_dtmin;
    m_factorIncMax=t_stepDesPtr->s_dtmax;
    m_factorInc1=t_stepDesPtr->s_dt0;
    m_factorInc2=0.0;
    m_factorFinal=t_stepDesPtr->s_t;
    m_expIters=t_stepDesPtr->s_expIters;
    m_arcLenMaxParam=t_stepDesPtr->s_arcLenMaxParam;
    m_arcLenMax=0.0;
    m_ifReadStepDes=true;
}
void LoadController::setMeshSysPtr(MeshSystem *t_meshSysPtr){
    m_meshPtr=t_meshSysPtr;
    m_loadVecPtr=&t_meshSysPtr->m_node_load;
    m_ifSetMeshSys=true;
}
double LoadController::nextArcLen(int t_mIters){
    return min(m_arcLen2*sqrt(m_expIters/t_mIters),m_arcLenMax);
}