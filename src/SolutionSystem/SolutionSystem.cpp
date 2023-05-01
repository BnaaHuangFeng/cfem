#include "SolutionSystem/SolutionSystem.h"
SolutionSystem::SolutionSystem():
                m_stepDesPtr(nullptr),m_meshSysPtr(nullptr),
                m_increI(0),m_iterI(0),m_mIter(0),
                m_rnorm(0.0),m_duNorm(0.0),m_rnorm0(0.0),m_uNorm(0.0),m_ifShowIterInfo(true)
{
    m_solutionCtx.s_bcsSysPtr=nullptr;
    m_solutionCtx.s_elmtSysPtr=nullptr;
    m_solutionCtx.s_loadCtrlPtr=nullptr;
    m_ifStepDesRead=false;
    m_ifSetMeshSysPtr=false;
    m_ifSolerInit=false;
    m_ifSolutionCtxInit=false;
}
SolutionSystem::SolutionSystem(StepDescriptiom *t_stepDesPtr):
                m_meshSysPtr(nullptr),
                m_increI(0),m_iterI(0),m_mIter(0),
                m_rnorm(0.0),m_duNorm(0.0),m_rnorm0(0.0),m_uNorm(0.0),m_ifShowIterInfo(true)
{
    m_solutionCtx.s_bcsSysPtr=nullptr;
    m_solutionCtx.s_elmtSysPtr=nullptr;
    m_solutionCtx.s_loadCtrlPtr=nullptr;
    m_ifStepDesRead=false;
    readStepDes(t_stepDesPtr);
    m_ifSetMeshSysPtr=false;
    m_ifSolerInit=false;
    m_ifSolutionCtxInit=false;
}
SolutionSystem::SolutionSystem(StepDescriptiom *t_stepDesPtr, MeshSystem *t_meshSysPtr):
                m_increI(0),m_iterI(0),m_mIter(0),
                m_rnorm(0.0),m_duNorm(0.0),m_rnorm0(0.0),m_uNorm(0.0),m_ifShowIterInfo(true)
{
    m_solutionCtx.s_bcsSysPtr=nullptr;
    m_solutionCtx.s_elmtSysPtr=nullptr;
    m_solutionCtx.s_loadCtrlPtr=nullptr;
    m_ifStepDesRead=false;
    m_ifSetMeshSysPtr=false;
    readStepDes(t_stepDesPtr);
    setMeshSysPtr(t_meshSysPtr);
    m_ifSolerInit=false;
    m_ifSolutionCtxInit=false;   
}
SolutionSystem::~SolutionSystem(){
    SNESDestroy(&m_snes);
}
PetscErrorCode SolutionSystem::init(StepDescriptiom *t_stepDesPtr,MeshSystem *t_meshSysPtr,ElementSystem *t_elmtSysPtr,BCsSystem *t_bcsSysPtr, LoadController *t_loadCtrlPtr){
    readStepDes(t_stepDesPtr);
    setMeshSysPtr(t_meshSysPtr);
    initStep(t_stepDesPtr);
    initSolutionCtx(t_elmtSysPtr, t_bcsSysPtr, t_loadCtrlPtr);
    initMonitorCtx();
    if(m_ifShowIterInfo){
        PetscCall(SNESMonitorSet(m_snes,monitorFunction,&m_monitorCtx,nullptr));
    }
    bindCallBack();
    return 0;
};
PetscErrorCode SolutionSystem::checkInit(){
    if(!m_ifStepDesRead){
        MessagePrinter::printErrorTxt("SolutionSystem: it has not read step description.");
        MessagePrinter::exitcfem();        
    }
    if(!m_ifSetMeshSysPtr){
        MessagePrinter::printErrorTxt("SolutionSystem: its relied mesh system has not yet been set.");
        MessagePrinter::exitcfem();
    }
    if(!m_ifSolutionCtxInit){
        MessagePrinter::printErrorTxt("SolutionSystem: its solution context has not yet been init.");
        MessagePrinter::exitcfem();        
    }
    if(!m_ifSolerInit){
        MessagePrinter::printErrorTxt("SolutionSystem: its Petsc solver has not yet been init.");
        MessagePrinter::exitcfem();          
    }
    return 0;
}
PetscErrorCode SolutionSystem::initStep(StepDescriptiom *t_stepDesPtr){
    if(!m_ifStepDesRead)
        readStepDes(t_stepDesPtr);
    initStep();
    return 0;
}
PetscErrorCode SolutionSystem::initStep(){
    if(m_ifSolerInit)return 0;
    if(!m_ifSetMeshSysPtr){
        MessagePrinter::printErrorTxt("SolutionSystem: its relied mesh system has not yet been set.");
        MessagePrinter::exitcfem();
    }
    PetscCall(SNESCreate(PETSC_COMM_WORLD,&m_snes));
    PetscCall(SNESSetDM(m_snes,m_meshSysPtr->m_dm));
    PetscCall(SNESGetKSP(m_snes,&m_ksp));
    PetscCall(KSPGMRESSetRestart(m_ksp,2500));
    PetscCall(KSPGetPC(m_ksp,&m_pc));
    PetscCall(PCSetType(m_pc,m_PCType));
    PetscCall(KSPSetType(m_ksp,m_KSPType));
    PetscCall(KSPSetFromOptions(m_ksp));
    PetscCall(PCSetFromOptions(m_pc));
    // extra setting for PC*******/
    /*****************************/
    if(strcmp(m_PCType,PCLU)==0)
        PetscCall(PCFactorSetMatSolverType(m_pc,MATSOLVERSUPERLU_DIST));
    PetscCall(PCFactorSetReuseOrdering(m_pc,PETSC_TRUE)); // ???
    // basic setting for SNES*****/
    //****************************/
    PetscCall(SNESSetTolerances(m_snes,m_absTol,m_relTol,m_uIncTol,m_maxIter,-1));
    PetscCall(SNESSetDivergenceTolerance(m_snes,-1));
    // for different types of SNES solver**/
    /**************************************/
    if(strcmp(m_SNESType,SNESNEWTONLS)==0){
        PetscCall(SNESSetType(m_snes,SNESNEWTONLS));
        PetscCall(SNESGetLineSearch(m_snes,&m_snesLinesearch));
        PetscCall(SNESLineSearchSetType(m_snesLinesearch,SNESLINESEARCHBT));
        PetscCall(SNESLineSearchSetOrder(m_snesLinesearch,3));
    }
    else if(strcmp(m_SNESType,SNESNEWTONTR)==0){
        PetscCall(SNESSetType(m_snes,SNESNEWTONTR));
    }
    else if(strcmp(m_SNESType,SNESNRICHARDSON)==0){
        PetscCall(SNESSetType(m_snes,SNESNRICHARDSON));
    }
    else if(strcmp(m_SNESType,SNESKSPONLY)==0){
        PetscCall(SNESSetType(m_snes,SNESKSPONLY));
    }
    else if(strcmp(m_SNESType,SNESNGMRES)==0){
        PetscCall(SNESSetType(m_snes,SNESNGMRES));
    }
    PetscCall(SNESSetFromOptions(m_snes));
    m_ifSolerInit=true;
    return 0;
}
void SolutionSystem::initSolutionCtx(ElementSystem *t_elmtSysPtr,BCsSystem *t_bcsSysPtr, LoadController *t_loadCtrlPtr){
    m_solutionCtx.s_elmtSysPtr=t_elmtSysPtr;
    m_solutionCtx.s_bcsSysPtr=t_bcsSysPtr;
    m_solutionCtx.s_loadCtrlPtr=t_loadCtrlPtr;
    m_ifSolutionCtxInit=true;
}
void SolutionSystem::initMonitorCtx(){
    m_monitorCtx.s_mIterPtr=&m_mIter;
    m_monitorCtx.s_duNormPtr=&m_duNorm;
    m_monitorCtx.s_iterIPtr=&m_iterI;
    m_monitorCtx.s_rnormPtr=&m_rnorm;
    m_monitorCtx.s_rnorm0Ptr=&m_rnorm0;
    m_monitorCtx.s_uNormPtr=&m_uNorm;
}
PetscErrorCode SolutionSystem::bindCallBack(){
    PetscCall(SNESSetFunction(m_snes,m_meshSysPtr->m_node_residual2,formFunction,&m_solutionCtx));
    PetscCall(SNESSetJacobian(m_snes,m_meshSysPtr->m_AMatrix2,m_meshSysPtr->m_AMatrix2,formJacobian,&m_solutionCtx));
    return 0;
}
void SolutionSystem::readStepDes(StepDescriptiom *t_stepDesPtr){
    m_stepDesPtr=t_stepDesPtr;
    m_algorithm=m_stepDesPtr->s_algorithm;
    m_SNESType=m_stepDesPtr->s_SNESType;
    m_KSPType=m_stepDesPtr->s_KSPType;
    m_PCType=m_stepDesPtr->s_PCType;
    m_maxIter=m_stepDesPtr->s_maxIterNum;
    m_absTol=m_stepDesPtr->s_absTol;
    m_relTol=m_stepDesPtr->s_relTol;
    m_uIncTol=m_stepDesPtr->s_duTol;
    m_ifStepDesRead=true;
}
PetscErrorCode SolutionSystem::showIterInfo(bool t_ifShowIterInfo){
    if(t_ifShowIterInfo){
        if(m_ifShowIterInfo)return 0;
        m_ifShowIterInfo=t_ifShowIterInfo;
        PetscCall(SNESMonitorSet(m_snes,monitorFunction,&m_monitorCtx,nullptr));
    }
    else{
        if(!m_ifShowIterInfo)return 0;
        m_ifShowIterInfo=t_ifShowIterInfo;
        PetscCall(SNESMonitorCancel(m_snes));
    }
    return 0;
}
PetscErrorCode SolutionSystem::run(bool t_ifLastConverged, bool *t_ifConverged, bool *t_ifcompleted){
    const int buffLen=100;
    char charBuff[buffLen];
    string increInfo;
    m_iterI=0;
    m_mIter=0;
    m_uNorm=0.0;
    m_duNorm=0.0;
    m_rnorm0=0.0;
    m_rnorm=0.0;
    if(m_solutionCtx.s_loadCtrlPtr->update(t_ifLastConverged,m_mIter)){
        *t_ifConverged=true;
        *t_ifcompleted=true;
        return 0;
    }
    m_solutionCtx.s_bcsSysPtr->setInitialSolution();
    PetscCall(SNESSolve(m_snes,NULL,m_solutionCtx.s_bcsSysPtr->m_uIncInitial));
    SNESConvergedReason converReason;
    PetscCall(SNESGetConvergedReason(m_snes,&converReason));
    switch (converReason)
    {
    case SNES_CONVERGED_FNORM_ABS:
        snprintf(charBuff,buffLen,"  Converged for |R| < abs tolerance, final iters=%3d", m_mIter);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        snprintf(charBuff,buffLen,"  incre %3d: |R|=%12.5e, abs-tol=%12.5e",m_increI,m_rnorm,m_absTol);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        *t_ifConverged=true;
        *t_ifcompleted=false;
        ++m_increI;
        return 0;
        break;
    case SNES_CONVERGED_FNORM_RELATIVE:
        snprintf(charBuff,buffLen,"  Converged for |R|/|R0| < rel tolerance, final iters=%3d", m_mIter);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        snprintf(charBuff,buffLen,"  incre %3d: |R|=%12.5e, |R0|=%12.5e, rel-tol=%12.5e",m_increI,m_rnorm,m_rnorm0,m_absTol);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        *t_ifConverged=true;
        *t_ifcompleted=false;
        ++m_increI;
        return 0;
        break;
    case SNES_CONVERGED_SNORM_RELATIVE:
        snprintf(charBuff,buffLen,"  Converged for |du|/|u| < stol, final iters=%3d", m_mIter);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        snprintf(charBuff,buffLen,"  incre %3d: |du|=%12.5e, |u|=%12.5e, stol=%12.5e",m_increI,m_duNorm,m_uNorm,m_uIncTol);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        *t_ifConverged=true;
        *t_ifcompleted=false;
        ++m_increI;
        return 0;
        break;        
    default:
        snprintf(charBuff,buffLen,"  Divergent, SNES nolinear solver diverged, incre=%3d, iters=%3d", m_increI, m_mIter);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        *t_ifConverged=false;
        *t_ifcompleted=false;
        return 0;        
        break;
    }
    return 0;
}
PetscErrorCode formFunction(SNES t_snes, Vec t_uInc, Vec t_function, void *ctx){
    if(ctx||t_snes){}
    SolutionCtx *ctxPtr=(SolutionCtx *)ctx;
    ctxPtr->s_elmtSysPtr->assemblRVec(&t_uInc,&t_function);
    ctxPtr->s_loadCtrlPtr->applyLoad(ctxPtr->s_loadCtrlPtr->m_factor1,&t_function);
    ctxPtr->s_bcsSysPtr->applyResidualBoundaryCondition(&t_function);
    return 0;
}
PetscErrorCode formJacobian(SNES t_snes, Vec t_uInc, Mat t_AMat, Mat t_PMat, void *ctx){
    if(ctx||t_snes){}
    SolutionCtx *ctxPtr=(SolutionCtx *)ctx;
    ctxPtr->s_elmtSysPtr->assembleAMatrix(&t_uInc,&t_PMat);
    ctxPtr->s_bcsSysPtr->applyJacobianBoundaryCondition(&t_PMat);
    if(t_AMat!=t_PMat){
        PetscCall(MatAssemblyBegin(t_AMat,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(t_AMat,MAT_FINAL_ASSEMBLY));
    }
    return 0;
}
PetscErrorCode monitorFunction(SNES snes, PetscInt its, PetscScalar norm, void *mctx){
    char buff[70];
    string str;
    MoniterCtx *mctxPtr=(MoniterCtx *)mctx;
    *mctxPtr->s_iterIPtr=its;
    *mctxPtr->s_mIterPtr=its;
    *mctxPtr->s_rnormPtr=norm;
    PetscCall(SNESGetSolutionNorm(snes,mctxPtr->s_uNormPtr));
    if(its==0){
        *mctxPtr->s_duNormPtr=*mctxPtr->s_uNormPtr;
    }
    else{
        PetscCall(SNESGetUpdateNorm(snes,mctxPtr->s_duNormPtr));
    }
    if(its==0){
        *mctxPtr->s_rnorm0Ptr=norm;
    }
    snprintf(buff,70," SNES solver: iters=%4d, |R|=%12.5e, |dU|=%12.5e",its,norm,*mctxPtr->s_duNormPtr);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    return 0;
}