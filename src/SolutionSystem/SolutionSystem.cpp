#include "SolutionSystem/SolutionSystem.h"
SolutionSystem::SolutionSystem():
                m_stepDesPtr(nullptr),m_meshSysPtr(nullptr),
                m_mDiverged(0),m_iterI(-1),m_mIter(0),
                m_rnorm(0.0),m_duNorm(0.0),m_rnorm0(0.0),m_uNorm(0.0),m_ifShowIterInfo(true)
{
    m_increI=0;
    m_solutionCtx.s_bcsSysPtr=nullptr;
    m_solutionCtx.s_elmtSysPtr=nullptr;
    m_solutionCtx.s_loadCtrlPtr=nullptr;
    m_ifStepDesRead=false;
    m_ifSetMeshSysPtr=false;
    m_ifSolverInit=false;
    m_ifSolutionCtxInit=false;
}
SolutionSystem::SolutionSystem(StepDescriptiom *t_stepDesPtr):
                m_meshSysPtr(nullptr),
                m_mDiverged(0),m_iterI(-1),m_mIter(0),
                m_rnorm(0.0),m_duNorm(0.0),m_rnorm0(0.0),m_uNorm(0.0),m_ifShowIterInfo(true)
{
    m_increI=0;
    m_solutionCtx.s_bcsSysPtr=nullptr;
    m_solutionCtx.s_elmtSysPtr=nullptr;
    m_solutionCtx.s_loadCtrlPtr=nullptr;
    m_ifStepDesRead=false;
    readStepDes(t_stepDesPtr);
    m_ifSetMeshSysPtr=false;
    m_ifSolverInit=false;
    m_ifSolutionCtxInit=false;
}
SolutionSystem::SolutionSystem(StepDescriptiom *t_stepDesPtr, MeshSystem *t_meshSysPtr):
                m_mDiverged(0),m_iterI(-1),m_mIter(0),
                m_rnorm(0.0),m_duNorm(0.0),m_rnorm0(0.0),m_uNorm(0.0),m_ifShowIterInfo(true)
{
    m_increI=0;
    m_solutionCtx.s_bcsSysPtr=nullptr;
    m_solutionCtx.s_elmtSysPtr=nullptr;
    m_solutionCtx.s_loadCtrlPtr=nullptr;
    m_ifStepDesRead=false;
    m_ifSetMeshSysPtr=false;
    readStepDes(t_stepDesPtr);
    setMeshSysPtr(t_meshSysPtr);
    m_ifSolverInit=false;
    m_ifSolutionCtxInit=false;   
}
SolutionSystem::~SolutionSystem(){
    SNESDestroy(&m_snes);
    if(m_arcLenSolverPtr)delete m_arcLenSolverPtr;
}
PetscErrorCode SolutionSystem::init(StepDescriptiom *t_stepDesPtr,MeshSystem *t_meshSysPtr,ElementSystem *t_elmtSysPtr,BCsSystem *t_bcsSysPtr, LoadController *t_loadCtrlPtr){
    readStepDes(t_stepDesPtr);
    setMeshSysPtr(t_meshSysPtr);
    initSolutionCtx(t_elmtSysPtr, t_bcsSysPtr, t_loadCtrlPtr);
    initMonitorCtx();
    initStep(t_stepDesPtr);
    m_ifShowIterInfo=true;
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
    if(!m_ifSolverInit){
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
    if(m_ifSolverInit)return 0;
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
    // PetscCall(PCFactorSetReuseOrdering(m_pc,PETSC_TRUE)); // ???
    // basic setting for SNES*****/
    //****************************/
    PetscCall(SNESSetTolerances(m_snes,m_absTol,m_relTol,m_uIncTol,m_maxIter,-1));
    PetscCall(SNESSetDivergenceTolerance(m_snes,-1));
    m_div_tol=-1;
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
    // for arc length method solver inition**/
    /****************************************/
    if(m_algorithm==AlgorithmType::ARCLENGTH_CYLENDER){
        m_arcLenSolverPtr=new ArcLengthSolver(m_stepDesPtr);
        m_arcLenSolverPtr->init(m_stepDesPtr,m_meshSysPtr,&m_ksp,&m_pc,&m_snesLinesearch);
        m_arcLenSolverPtr->setLoadVecPtr(m_solutionCtx.s_loadCtrlPtr->getLoadVecPtr());
    }
    m_ifSolverInit=true;
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
    switch(m_algorithm){
        case AlgorithmType::STANDARD:
            PetscCall(SNESSetFunction(m_snes,m_meshSysPtr->m_node_residual2,formFunction,&m_solutionCtx));
            PetscCall(SNESSetJacobian(m_snes,m_meshSysPtr->m_AMatrix2,m_meshSysPtr->m_AMatrix2,formJacobian,&m_solutionCtx));
            PetscCall(SNESMonitorSet(m_snes,monitorFunction,&m_monitorCtx,nullptr));
            break;
        case AlgorithmType::ARCLENGTH_CYLENDER:
            m_arcLenSolverPtr->setFunction(&m_meshSysPtr->m_node_residual2,formFunctionArcLen,&m_solutionCtx);
            m_arcLenSolverPtr->setJacobian(&m_meshSysPtr->m_AMatrix2,&m_meshSysPtr->m_AMatrix2,formJacobianArcLen,&m_solutionCtx);
            m_arcLenSolverPtr->setLoad(applyLoadArcLen,&m_solutionCtx);
            m_arcLenSolverPtr->monitorSet(monitorFunctionArcLen,&m_monitorCtx);
            break;
    }
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
    string increInfo;
    m_uNorm=0.0;
    m_duNorm=0.0;
    m_rnorm0=0.0;
    m_rnorm=0.0;
    switch(m_algorithm){
        case AlgorithmType::STANDARD:
            if(m_solutionCtx.s_loadCtrlPtr->update(t_ifLastConverged)){
                *t_ifConverged=true;
                *t_ifcompleted=true;
                return 0;
            }
            else{
                *t_ifcompleted=false;
            }
            break;
        default:
            break;
    }
    if(m_algorithm==AlgorithmType::STANDARD){
        m_solutionCtx.s_bcsSysPtr->setInitialSolution(m_solutionCtx.s_loadCtrlPtr->m_factorInc1);
        PetscCall(SNESSolve(m_snes,NULL,m_solutionCtx.s_bcsSysPtr->m_uIncInitial));
    }
    else if(m_algorithm==AlgorithmType::ARCLENGTH_CYLENDER){
        m_solutionCtx.s_bcsSysPtr->setInitialSolution(1.0);
        PetscScalar arcLen=m_solutionCtx.s_loadCtrlPtr->getArcLen();
        if(m_increI==1){
            arcLen=m_arcLenSolverPtr->getIntialArcLen(m_solutionCtx.s_loadCtrlPtr->m_factorInc1);
            m_solutionCtx.s_loadCtrlPtr->m_factorInc1=0.0;
            m_solutionCtx.s_loadCtrlPtr->setInitialArc(arcLen);
        }
        else if(m_increI==2){
            m_arcLenSolverPtr->setLastSolutionPtr(&m_meshSysPtr->m_nodes_uInc2);
        }
        m_arcLenSolverPtr->solve(arcLen);
    }
    SNESConvergedReason converReason;
    if(m_algorithm==AlgorithmType::STANDARD)
        PetscCall(SNESGetConvergedReason(m_snes,&converReason));
    else if(m_algorithm==AlgorithmType::ARCLENGTH_CYLENDER){
        PetscScalar facInc=m_arcLenSolverPtr->getFactorInc();
        int iterNum=m_arcLenSolverPtr->getIterNum();
        if(m_solutionCtx.s_loadCtrlPtr->update(t_ifLastConverged,facInc,iterNum)){
            *t_ifcompleted=true;
        }           
        else{
            *t_ifcompleted=false;
        }
        m_arcLenSolverPtr->getConvergedReason(&converReason);     
    }
    printConvergedReason(converReason);
    if(converReason>0){// for converged case
        *t_ifConverged=true;
        m_mDiverged=0;        
    }
    else if(converReason<0){// for diverged reason
        *t_ifConverged=false;
        ++m_mDiverged;
        if(m_mDiverged>=m_maxAttempt){
            MessagePrinter::printErrorTxt("too manay attempt to get converged!");
            MessagePrinter::exitcfem();
        }        
    }
    else{
        MessagePrinter::printErrorTxt("unexpected SNESConvergedReason value = 0!");
        MessagePrinter::exitcfem();
    }
    return 0;
}

PetscErrorCode formFunction(SNES t_snes, Vec t_uInc, Vec t_function, void *ctx){
    if(t_snes){}
    SolutionCtx *ctxPtr=(SolutionCtx *)ctx;
    // for debug
    // MessagePrinter::printTxt("incremental u before assemble residual:");
    // PetscCall(VecView(t_uInc,PETSC_VIEWER_STDOUT_WORLD));
    ctxPtr->s_elmtSysPtr->assemblRVec(&t_uInc,&t_function);
    // for debug
    // MessagePrinter::printTxt("function after assemble bcs:");
    // PetscCall(VecView(t_function,PETSC_VIEWER_STDOUT_WORLD));
    ctxPtr->s_loadCtrlPtr->applyLoad(ctxPtr->s_loadCtrlPtr->m_factor1,&t_function);
    ctxPtr->s_bcsSysPtr->applyResidualBoundaryCondition(&t_function);
    // for debug
    // MessagePrinter::printTxt("function after apply bcs:");
    // PetscCall(VecView(t_function,PETSC_VIEWER_STDOUT_WORLD));
    return 0;
}
PetscErrorCode formFunctionArcLen(ArcLengthSolver *t_solverPtr, Vec *t_uInc, Vec *t_function, void *ctx){
    SolutionCtx *ctxPtr=(SolutionCtx *)ctx;
    ctxPtr->s_elmtSysPtr->assemblRVec(t_uInc,t_function);
    double loadFactor=ctxPtr->s_loadCtrlPtr->m_factor1+t_solverPtr->getFactorInc();
    ctxPtr->s_loadCtrlPtr->applyLoad(loadFactor,t_function);
    ctxPtr->s_bcsSysPtr->applyResidualBoundaryCondition(t_function);
    PetscCall(VecScale(*t_function,-1.0));
    // PetscCall(VecView(*t_uInc,PETSC_VIEWER_STDOUT_WORLD));
    return 0;
}
PetscErrorCode applyLoadArcLen(Mat *AMatrixPtr, Vec *t_b, void *ctx){
    SolutionCtx *ctxPtr=(SolutionCtx *)ctx;
    ctxPtr->s_bcsSysPtr->applyBoundaryConditionArc(AMatrixPtr,t_b);
    ctxPtr->s_loadCtrlPtr->applyLoad(-1.0,t_b);
    return 0;
}

PetscErrorCode formJacobian(SNES t_snes, Vec t_uInc, Mat t_AMat, Mat t_PMat, void *ctx){
    if(t_snes){}
    SolutionCtx *ctxPtr=(SolutionCtx *)ctx;
    // for debug
    // MessagePrinter::printTxt("incremental u before assemble Jacobian:");
    // PetscCall(VecView(t_uInc,PETSC_VIEWER_STDOUT_WORLD));
    ctxPtr->s_elmtSysPtr->assembleAMatrix(&t_uInc,&t_PMat);
    if(ctxPtr->s_bcsSysPtr->m_drclt_method==DirichletMethod::SETLARGE)
        ctxPtr->s_bcsSysPtr->update_penalty(&t_PMat);
    // for debug:
    // MessagePrinter::printTxt("J(u) after assembleAMatrix:");
    // PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE));
    // PetscCall(MatView(t_AMat,PETSC_VIEWER_STDOUT_WORLD));
    ctxPtr->s_bcsSysPtr->applyJacobianBoundaryCondition(&t_PMat);
    if(t_AMat!=t_PMat){
        PetscCall(MatAssemblyBegin(t_AMat,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(t_AMat,MAT_FINAL_ASSEMBLY));
    }
    // MessagePrinter::printTxt("J(u) after applyJacobianBoundaryCondition:");
    // MatView(t_PMat,PETSC_VIEWER_STDOUT_WORLD);
    // MatView(t_PMat,PETSC_VIEWER_STDOUT_WORLD);
    // PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB));
    // PetscCall(MatView(t_AMat,PETSC_VIEWER_STDOUT_WORLD));
    return 0;
}
PetscErrorCode formJacobianArcLen(ArcLengthSolver *t_solverPtr, Vec *t_uInc, Mat *t_AMat, Mat *t_PMat, void *ctx){
    if(t_solverPtr){}
    SolutionCtx *ctxPtr=(SolutionCtx *)ctx;
    ctxPtr->s_elmtSysPtr->assembleAMatrix(t_uInc,t_PMat);
    if(ctxPtr->s_bcsSysPtr->m_drclt_method==DirichletMethod::SETLARGE)
        ctxPtr->s_bcsSysPtr->update_penalty(t_PMat);
    if(t_AMat!=t_PMat){
        PetscCall(MatAssemblyBegin(*t_AMat,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(*t_AMat,MAT_FINAL_ASSEMBLY));
    }
    return 0;
}
PetscErrorCode monitorFunction(SNES snes, PetscInt its, PetscScalar norm, void *mctx){
    // char buff[70];
    // string str;
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
    // snprintf(buff,70," SNES solver: iters=%4d, |R|=%12.5e, |dU|=%12.5e",its,norm,*mctxPtr->s_duNormPtr);
    // str=buff;
    // MessagePrinter::printNormalTxt(str);
    return 0;
}
PetscErrorCode monitorFunctionArcLen(ArcLengthSolver *t_solverPtr, PetscInt its, PetscScalar norm, void *mctx){
    // char buff[70];
    // string str;
    if(t_solverPtr){}
    MoniterCtx *mctxPtr=(MoniterCtx *)mctx;
    *mctxPtr->s_iterIPtr=its;
    *mctxPtr->s_mIterPtr=its;
    *mctxPtr->s_rnormPtr=norm;
    if(its==0){
        *mctxPtr->s_rnorm0Ptr=norm;
    }
    return 0;
}
PetscErrorCode SolutionSystem::printConvergedReason(SNESConvergedReason converReason){
    const int buffLen=200;
    char charBuff[buffLen];
    string increInfo;
    snprintf(charBuff,buffLen,"Increment %4d: attempt=%3d, incre-t=%12.5e, total-t=%12.5e",m_increI,m_mDiverged,
            converReason>0?m_solutionCtx.s_loadCtrlPtr->m_factorInc1:0.0,
            converReason>0?m_solutionCtx.s_loadCtrlPtr->m_factor1:m_solutionCtx.s_loadCtrlPtr->m_factor2);
    increInfo=charBuff;
    MessagePrinter::printNormalTxt(increInfo);    
    switch (converReason)
    {
/*************************/
/* converged reason print*/
/*************************/
    case SNES_CONVERGED_FNORM_ABS:
        snprintf(charBuff,buffLen,"Converged for |F| < abs tolerance, final iters=%3d",m_mIter);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        snprintf(charBuff,buffLen,"  |F|=%12.5e, abs-tol=%12.5e",m_rnorm,m_absTol);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;
    case SNES_CONVERGED_FNORM_RELATIVE:
        snprintf(charBuff,buffLen,"Converged for |F|/|F0| < rel tolerance, final iters=%3d", m_mIter);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        snprintf(charBuff,buffLen,"  |F|=%12.5e, |F0|=%12.5e, rel-tol=%12.5e",m_rnorm,m_rnorm0,m_absTol);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;
    case SNES_CONVERGED_SNORM_RELATIVE:
        snprintf(charBuff,buffLen,"Converged for |du|/|u| < stol, final iters=%3d", m_mIter);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        snprintf(charBuff,buffLen,"  |du|=%12.5e, |u|=%12.5e, stol=%12.5e",m_duNorm,m_uNorm,m_uIncTol);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break; 
    case SNES_CONVERGED_ITS:
        snprintf(charBuff,buffLen,"Maximum iterations reached");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;      
    case SNES_BREAKOUT_INNER_ITER:
        snprintf(charBuff,buffLen,"Flag to break out of inner loop after checking custom convergence (it is used in multi-phase flow when state changes)");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;  
/*************************/
/* diverged reason print*/
/*************************/    
    case SNES_DIVERGED_FUNCTION_DOMAIN:
        snprintf(charBuff,buffLen,"the new u location passed the function is not in the domain of F");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;              
    case SNES_DIVERGED_FUNCTION_COUNT:
        snprintf(charBuff,buffLen,"The user provided function has been called more times than the maximum");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;    
    case SNES_DIVERGED_LINEAR_SOLVE:
        snprintf(charBuff,buffLen,"the linear solve failed");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;
    case SNES_DIVERGED_FNORM_NAN:
        snprintf(charBuff,buffLen,"the 2-norm of the current function evaluation is not-a-number (NaN), this is usually caused by a division of 0 by 0.");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;              
    case SNES_DIVERGED_MAX_IT:
        snprintf(charBuff,buffLen,"the solver reached the maximum number of iterations (%d) without satisfying any convergence criteria.", m_maxIter);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);        
        break;
    case SNES_DIVERGED_LINE_SEARCH:
        snprintf(charBuff,buffLen,"the line search failed.");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;   
    case SNES_DIVERGED_INNER:
        snprintf(charBuff,buffLen,"Inner solve failed.");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;         
    case SNES_DIVERGED_LOCAL_MIN:
        snprintf(charBuff,buffLen,"|| J^T b || is small, implies converged to local minimum of F().");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;     
    case SNES_DIVERGED_DTOL:
        snprintf(charBuff,buffLen,"|F| > div-tol * |F0|.");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        snprintf(charBuff,buffLen,"  |F|=%12.5e, |F0|=%12.5e, div-tol=%12.5e",m_rnorm,m_rnorm0,m_div_tol);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);
        break;   
    case SNES_DIVERGED_JACOBIAN_DOMAIN:
        snprintf(charBuff,buffLen,"Jacobian calculation does not make sense.");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);    
        break;
    case SNES_DIVERGED_TR_DELTA:
        snprintf(charBuff,buffLen,"SNESConvergedReason is SNES_DIVERGED_TR_DELTA.");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);    
        break;        
    case SNES_CONVERGED_ITERATING:
        snprintf(charBuff,buffLen,"Iteration is continuing.");
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);        
        break;
    default:
        snprintf(charBuff,buffLen,"Unknown SNESConvergedReason value = %d.", converReason);
        increInfo=charBuff;
        MessagePrinter::printNormalTxt(increInfo);   
        break;
    }
    MessagePrinter::printDashLine();
    return 0;
}