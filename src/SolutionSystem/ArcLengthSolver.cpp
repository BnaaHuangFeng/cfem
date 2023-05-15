#include "SolutionSystem/ArcLengthSolver.h"
#include "MathUtils/PetsExtension.h"
ArcLengthSolver::ArcLengthSolver(){
    m_stepDesPtr=nullptr;
    m_div_tol=-1;
    m_arcLen=0;
    signum=1.0;
    m_factorInc=0.0;
    m_factorIter=0.0;
    m_kspPtr=nullptr; m_pcPtr=nullptr; m_snesLinesearchPtr=nullptr;
    m_mDiverged=0;
    m_mIter=0;
    m_rnorm=0.0;
    m_duNorm=0.0;
    m_rnorm0=0.0;
    m_uNorm=0.0;
    m_ifStepDesRead=false; m_ifSolverInit=false; m_ifMeshSysSet=false; m_ifFunctionSet=false; m_ifJacobianSet=false;
    m_rPtr=nullptr; m_loadPtr=nullptr; m_AMatPtr=nullptr; m_PMatPtr=nullptr; m_uInc2Ptr=nullptr;
    m_converReason = SNES_CONVERGED_ITERATING;
    m_functionCtx=nullptr; m_jacobianCtx=nullptr;
}
ArcLengthSolver::ArcLengthSolver(StepDescriptiom *t_stepDesPtr){
    m_stepDesPtr=nullptr;
    m_div_tol=-1;
    m_arcLen=0;
    signum=1.0;
    m_factorInc=0.0;
    m_factorIter=0.0;
    m_kspPtr=nullptr; m_pcPtr=nullptr; m_snesLinesearchPtr=nullptr; m_uInc2Ptr=nullptr;
    m_mDiverged=0;
    m_mIter=0;
    m_rnorm=0.0;
    m_duNorm=0.0;
    m_rnorm0=0.0;
    m_uNorm=0.0;
    m_ifStepDesRead=false; m_ifSolverInit=false; m_ifMeshSysSet=false; m_ifFunctionSet=false; m_ifJacobianSet=false;
    m_ifLoadVecPtrSet=false; m_ifLastSolutionPtrSet=false;
    m_rPtr=nullptr; m_loadPtr=nullptr; m_AMatPtr=nullptr; m_PMatPtr=nullptr;
    m_converReason = SNES_CONVERGED_ITERATING;
    readStepDes(t_stepDesPtr);
    m_functionCtx=nullptr; m_jacobianCtx=nullptr;
}
ArcLengthSolver::~ArcLengthSolver(){
    if(m_ifMeshSysSet){
        m_meshSysPtr->destroyGlobalVec(&m_uIterRes);
        m_meshSysPtr->destroyGlobalVec(&m_uIterTangent);
        m_meshSysPtr->destroyGlobalVec(&m_uInc);
        m_meshSysPtr->destroyGlobalVec(&m_b);
    }
}
PetscErrorCode ArcLengthSolver::init(StepDescriptiom *t_stepDesPtr,MeshSystem *t_meshSysPtr, KSP *t_kspPtr, PC *t_pcPtr, SNESLineSearch *t_snesLinesearchPtr){
    if(!m_ifStepDesRead)readStepDes(t_stepDesPtr);
    setMeshSystem(t_meshSysPtr);
    m_kspPtr=t_kspPtr;
    m_pcPtr=t_pcPtr;
    m_snesLinesearchPtr=t_snesLinesearchPtr;
    m_ifSolverInit=true;
    return 0;
}
PetscErrorCode ArcLengthSolver::setFunction(Vec *rPtr, PetscErrorCode (*f)(ArcLengthSolver *t_solverPtr, Vec *t_uInc, Vec *t_function, void *ctx),void *t_ctx){
    m_rPtr=rPtr;
    m_functionCal=f;
    m_functionCtx=t_ctx;
    m_ifFunctionSet=true;
    return 0;
}
PetscErrorCode ArcLengthSolver::setJacobian(Mat *AmatPtr, Mat *PmatPtr, PetscErrorCode (*J)(ArcLengthSolver *t_solverPtr, Vec *t_uInc, Mat *t_AMat, Mat *t_PMat, void *ctx),void *t_ctx){
    m_AMatPtr=AmatPtr;
    m_PMatPtr=PmatPtr;
    m_jacobianCal=J;
    m_jacobianCtx=t_ctx;
    m_ifJacobianSet=true;
    return 0;
}

PetscErrorCode ArcLengthSolver::solve(double t_arcLen){
    checkInit();
    m_arcLen=t_arcLen;
    PetscCall(VecZeroEntries(m_uInc));
    // for debug
    // PetscCall(VecView(m_uInc,PETSC_VIEWER_STDOUT_WORLD));
    m_factorInc=0.0;
    m_mIter=0;
    while(m_mIter<=m_maxIter){
        (*m_functionCal)(this,&m_uInc,m_rPtr,m_functionCtx);
        if(m_mIter>0){
            updateConvergenceReason();
            (*m_monitor)(this,m_mIter,m_rnorm,m_monitorCtx);
            if(m_converReason>0) return 0;
            else if(m_converReason<0) return 0;
        }
        (*m_jacobianCal)(this,&m_uInc,m_AMatPtr,m_PMatPtr,m_jacobianCtx);
        (*m_applyLoad)(m_AMatPtr,&m_b,m_LoadCtx);
        KSPConvergedReason kspReason;
        calTangentIter(&kspReason);
        if(kspReason<0){
            m_converReason=SNES_DIVERGED_LINEAR_SOLVE;
            break;
        }
        calResidualIter(&kspReason);
        if(kspReason<0){
            m_converReason=SNES_DIVERGED_LINEAR_SOLVE;
            break;
        }
        // PetscCall(VecView(m_uIterRes,PETSC_VIEWER_STDOUT_WORLD));
        // PetscCall(VecView(m_uIterTangent,PETSC_VIEWER_STDOUT_WORLD));
        PetscScalar sigNum;
        if(m_mIter==0)
            getSigNum(&sigNum);
        updateFactor(sigNum);
        if(m_converReason==SNES_DIVERGED_INNER)break;
        PetscCall(VecAXPBYPCZ(m_uInc,1.0,m_factorIter,1.0,m_uIterRes,m_uIterTangent));
        ++m_mIter;
    }
    return 0;
}
PetscErrorCode ArcLengthSolver::getSigNum(double *t_sigNum){
    if(m_mIter==0){
        if(!m_ifLastSolutionPtrSet){ // it means this is the first increment
            *t_sigNum=1.0;
            return 0;
        }
        PetscScalar vecDot;
        PetscCall(VecDot(*m_uInc2Ptr,m_uIterTangent,&vecDot));
        if(vecDot>0) *t_sigNum=1.0;
        else *t_sigNum=-1.0;
        return 0;
    }
    else{
        MessagePrinter::printErrorTxt("sign predictor can only be called in the first iteration");
        MessagePrinter::exitcfem();
    }
    return 0;
}
PetscErrorCode ArcLengthSolver::updateFactor(double sigNum){
    if(m_mIter==0){
        PetscScalar iterTangentNorm;
        PetscCall(VecNorm(m_uIterTangent,NORM_2,&iterTangentNorm));
        m_factorIter=sigNum*m_arcLen/iterTangentNorm;
        m_factorInc=m_factorIter;
    }
    else{
        PetscScalar a,b,c;
        PetscScalar tdt, idt, rdt, idi, idr, rdr;
        PetscCall(VecDot(m_uIterTangent,m_uIterTangent,&tdt));
        PetscCall(VecDot(m_uInc,m_uIterTangent,&idt));
        PetscCall(VecDot(m_uIterRes,m_uIterTangent,&rdt));
        PetscCall(VecDot(m_uInc,m_uInc,&idi));
        PetscCall(VecDot(m_uInc,m_uIterRes,&idr));
        PetscCall(VecDot(m_uIterRes,m_uIterRes,&rdr));
        a=tdt;
        b=2.0*(idt+rdt);
        c=idi+2*idr+rdr-m_arcLen*m_arcLen;
        PetscScalar root1, root2;
        PetscInt rootNum=0;
        rootNum = PetscExtension::solQua(a,b,c,&root1, &root2);
        if(rootNum==0){
            m_converReason=SNES_DIVERGED_INNER;
            return 0;
        }
        else if(rootNum==1){
            m_factorIter=root1;
        }
        else if(rootNum==2){
            PetscScalar cos1=0.0, cos2=0.0;
            cos1=idi+idr+idt*root1;
            cos2=idi+idr+idt*root2;
            m_factorIter=cos1>cos2?root1:root2;
        }
        m_factorInc+=m_factorIter;
    }
    return 0;
}
double ArcLengthSolver::getIntialArcLen(double factorInc0){
    (*m_jacobianCal)(this,&m_uInc,m_AMatPtr,m_PMatPtr,m_jacobianCtx);
    (*m_applyLoad)(m_AMatPtr,&m_b,m_LoadCtx);
    KSPConvergedReason kspReason;
    calTangentIter(&kspReason);
    PetscCall(VecAXPY(m_uInc,factorInc0,m_uIterTangent));
    double initialArcLen=0.0;
    VecNorm(m_uInc,NORM_2,&initialArcLen);
    PetscCall(VecZeroEntries(m_uInc));
    return initialArcLen;    
}
PetscErrorCode ArcLengthSolver::calTangentIter(KSPConvergedReason *reasonPtr){
    PetscCall(KSPSetOperators(*m_kspPtr,*m_AMatPtr,*m_PMatPtr));
    PetscCall(KSPSolve(*m_kspPtr,m_b,m_uIterTangent));
    PetscCall(KSPGetConvergedReason(*m_kspPtr,reasonPtr));
    if(*reasonPtr<0)return 0;
    return 0;
}
PetscErrorCode ArcLengthSolver::calResidualIter(KSPConvergedReason *reasonPtr){
    PetscCall(KSPSetOperators(*m_kspPtr,*m_AMatPtr,*m_PMatPtr));
    PetscCall(KSPSolve(*m_kspPtr,*m_rPtr,m_uIterRes));
    PetscCall(KSPGetConvergedReason(*m_kspPtr,reasonPtr));
    return 0;
}
PetscErrorCode ArcLengthSolver::getSolution(Vec *t_solution){
    PetscCall(VecCopy(m_uInc,*t_solution));
    return 0;
}
PetscErrorCode ArcLengthSolver::getConvergedReason(SNESConvergedReason *t_converReasonPtr){
    *t_converReasonPtr=m_converReason;
    return 0;
}
void ArcLengthSolver::readStepDes(StepDescriptiom *t_stepDesPtr){
    m_stepDesPtr=t_stepDesPtr;
    m_algorithm=m_stepDesPtr->s_algorithm;
    m_maxIter=m_stepDesPtr->s_maxIterNum;
    m_absTol=m_stepDesPtr->s_absTol;
    m_relTol=m_stepDesPtr->s_relTol;
    m_uIncTol=m_stepDesPtr->s_duTol;
    m_ifStepDesRead=true;
}
PetscErrorCode ArcLengthSolver::checkInit(){
    if(!m_ifStepDesRead){
        MessagePrinter::printErrorTxt("ArcLengthSolver: step description was not read.");
        MessagePrinter::exitcfem();
    }
    if(!m_ifSolverInit){
        MessagePrinter::printErrorTxt("ArcLengthSolver: solver was not inited.");
        MessagePrinter::exitcfem();        
    }
    if(!m_ifMeshSysSet){
        MessagePrinter::printErrorTxt("ArcLengthSolver: mesh system was not set.");
        MessagePrinter::exitcfem();        
    }
    if(!m_ifFunctionSet){
        MessagePrinter::printErrorTxt("ArcLengthSolver: callback function for residual calculation was not set.");
        MessagePrinter::exitcfem();        
    }
    if(!m_ifJacobianSet){
        MessagePrinter::printErrorTxt("ArcLengthSolver: callback function for jacobian calculation was not set.");
        MessagePrinter::exitcfem();        
    }    
    if(!m_loadPtr){
        MessagePrinter::printErrorTxt("ArcLengthSolver: load Vec Pointer was not set.");
        MessagePrinter::exitcfem();            
    }
    return 0;
}
void ArcLengthSolver::setMeshSystem(MeshSystem *t_meshSysPtr){
    if(m_ifMeshSysSet){
        MessagePrinter::printErrorTxt("can not set mesh system pointer of arc length solver repeatedly.");
        MessagePrinter::exitcfem();
    }
    m_meshSysPtr=t_meshSysPtr;
    m_meshSysPtr->createGlobalVec(&m_uIterRes);
    m_meshSysPtr->createGlobalVec(&m_uIterTangent);
    m_meshSysPtr->createGlobalVec(&m_uInc);
    m_meshSysPtr->createGlobalVec(&m_b);
    // for debug
    // VecView(m_uIterRes,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(m_uIterTangent,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(m_uInc,PETSC_VIEWER_STDOUT_WORLD);
    m_ifMeshSysSet=true;
}
PetscErrorCode ArcLengthSolver::updateConvergenceReason(){
    if(m_mIter==1){
        PetscCall(VecNorm(*m_rPtr,NORM_2,&m_rnorm0));
    }
    PetscCall(VecNorm(*m_rPtr,NORM_2,&m_rnorm));
    if(m_rnorm<m_absTol){
        m_converReason=SNES_CONVERGED_FNORM_ABS;
        return 0;
    }
    if(m_rnorm<m_relTol*m_rnorm0){
        m_converReason=SNES_CONVERGED_FNORM_RELATIVE;
        return 0;
    }
    if(m_div_tol>0&&m_rnorm>m_div_tol*m_rnorm0&&m_mIter>1){
        m_converReason=SNES_DIVERGED_DTOL;
        return 0;
    }
    if(m_mIter>=m_maxIter){
        m_converReason=SNES_DIVERGED_MAX_IT;
        return 0;
    }
    m_converReason=SNES_CONVERGED_ITERATING;
    return 0;
}