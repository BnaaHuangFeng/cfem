#include "LoadController/LoadController.h"
LoadController::LoadController():
m_stepDesPtr(nullptr),m_meshPtr(nullptr),
m_factorFinal(1.0),m_growthRatio(1.1),m_cutbackRatio(0.25),m_factorIncMin(0.001),m_factorIncMax(1.0),
m_mGrowth(5),m_convergeCnt(0),m_iterCnt2(0),m_iterCnt1(0),
m_factor2(0.0),m_factorInc2(0.0),m_factor1(0.0),m_factorInc1(0.0),m_l0(0.1){
    m_ifReadStepDes=false;
    m_ifSetMeshSys=false;
    m_ifHasNodeLoad=false;
}
LoadController::LoadController(StepDescriptiom *t_stepDesPtr):
m_meshPtr(nullptr),
m_mGrowth(5),m_convergeCnt(0),m_iterCnt2(0),m_iterCnt1(0),
m_factor2(0.0),m_factorInc2(0.0),m_factor1(0.0),m_factorInc1(0.0){
    readStepDes(t_stepDesPtr);
    m_ifSetMeshSys=false;
    m_ifHasNodeLoad=false;
}
LoadController::LoadController(StepDescriptiom *t_stepDesPtr, MeshSystem *t_meshSysPtr):
m_mGrowth(5),m_convergeCnt(0),m_iterCnt2(0),m_iterCnt1(0),
m_factor2(0.0),m_factorInc2(0.0),m_factor1(0.0),m_factorInc1(0.0){
    readStepDes(t_stepDesPtr);
    setMeshSysPtr(t_meshSysPtr);
    m_ifHasNodeLoad=false;
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
    m_factorInc2=t_stepDesPtr->s_dt0;
    m_factorFinal=t_stepDesPtr->s_t;
    m_ifReadStepDes=true;
}
void LoadController::setMeshSysPtr(MeshSystem *t_meshSysPtr){
    m_meshPtr=t_meshSysPtr;
    m_loadVecPtr=&t_meshSysPtr->m_node_load;
    m_ifSetMeshSys=true;
}