#pragma once
#include "InputSystem/DescriptionInfo.h"
#include "SolutionSystem/SolutionCtx.h"
/**
 * form the global function Vec for SNES solver
 * @param t_uInc > state at which to evaluate residual
 * @param t_ctx > user-defined function context
 * @param f < vector to put residual (function value)
*/
PetscErrorCode formFunction(SNES t_snes, Vec t_uInc, Vec t_function, void *ctx);
/**
 * form the global jacobian Mat for SNES solver
 * @param t_uInc > input vector, the jacobian is to be computed at this value
 * @param t_ctx > user-defined jacobian context
 * @param t_AMat < the matrix that defines the jacobian
 * @param t_PMat < the matrix to be used in constructing the preconditioner, usually the same as t_AMat
*/
PetscErrorCode formJacobian(SNES t_snes, Vec t_uInc, Mat t_AMat, Mat t_PMat, void *ctx);
/**
 * functional form passed to SNESMonitorSet() to monitor convergence of nolinear solver
 * @param snes > the SNES context
 * @param its > iteration number
 * @param norm > 2 norm function value (may be estimated)
 * @param mctx > [optional] monitor context
*/
PetscErrorCode monitorFunction(SNES snes, PetscInt its, PetscScalar norm, void *mctx);
struct MoniterCtx{
    int *s_mIterPtr;                    /**< ptr to iteration num of the latest convergence*/
    int *s_iterIPtr;                    /**< ptr to current iteration id*/
    PetscScalar *s_rnormPtr;            /**< ptr to 2 norm function value (may be estimated)*/
    PetscScalar *s_duNormPtr;           /**< ptr to |du| of current iter*/
    PetscScalar *s_rnorm0Ptr;           /**< ptr to 2 norm function value of iteration 0*/
    PetscScalar *s_uNormPtr;            /**< ptr to 2 norm of solution*/
};

class SolutionSystem
{
private:
    StepDescriptiom *m_stepDesPtr;      /**< ptr to step description*/
    MeshSystem *m_meshSysPtr;           /**< mesh system it rely on*/
    AlgorithmType m_algorithm;
    SNESType m_SNESType;
    KSPType m_KSPType;
    PCType m_PCType;
    int m_maxIter;                  /**< max iteration num limit*/
    PetscScalar m_absTol;           /**< absolute tolerance*/
    PetscScalar m_relTol;           /**< relative tolerance*/
    PetscScalar m_uIncTol;          /**< incremental dof tolerance*/
    PetscScalar m_div_tol;            /**< the divergence tolerance. Use -1 to deactivate the test, default is 1e4*/
    KSP m_ksp;                      /**< ksp solver*/
    PC m_pc;                        /**< preconditioner*/
    SNESLineSearch m_snesLinesearch;/**< snes line search solver*/
    bool m_ifSolerInit;             /**< if the Petsc solver be inited*/
    bool m_ifStepDesRead;           /**< if has read step description*/
    bool m_ifSetMeshSysPtr;         /**< if has set mesh system ptr*/
    bool m_ifSolutionCtxInit;       /**< if solution context inited*/

    const int m_maxAttempt=10;      /**< max diverged increments tolerance*/
    int m_mDiverged;             /**< consecutive divergence num*/
    int m_increI;                   /**< current increment id*/
    int m_iterI;                    /**< current iteration id*/
    int m_mIter;                    /**< iteration num of the latest convergence*/
    PetscScalar m_rnorm;            /**< 2 norm function value (may be estimated)*/
    PetscScalar m_duNorm;           /**< |du| of current iter*/
    PetscScalar m_rnorm0;           /**< 2 norm function value of iteration 0*/
    PetscScalar m_uNorm;            /**< 2 norm of solution*/
private:
    void readStepDes(StepDescriptiom *t_stepDesPtr);
    PetscErrorCode initStep();
    PetscErrorCode initStep(StepDescriptiom *t_stepDesPtr);
    void initSolutionCtx(ElementSystem *t_elmtSysPtr,BCsSystem *t_bcsSysPtr, LoadController *t_loadCtrlPtr);
    void initMonitorCtx();
    inline void setMeshSysPtr(MeshSystem *t_meshSysPtr){
        m_meshSysPtr=t_meshSysPtr;
        m_ifSetMeshSysPtr=true;
    }
    /**
     * bind call back function for SNES solver
    */
    PetscErrorCode bindCallBack();
    /**
     * print the increment's SNES calculation's convergence(divergence) reason
    */
    PetscErrorCode printConvergedReason(SNESConvergedReason converReason);
public:
    SolutionSystem();
    SolutionSystem(StepDescriptiom *t_stepDesPtr);
    SolutionSystem(StepDescriptiom *t_stepDesPtr, MeshSystem *t_meshSysPtr);
    ~SolutionSystem();
    PetscErrorCode init(StepDescriptiom *t_stepDesPtr,MeshSystem *t_meshSysPtr,ElementSystem *t_elmtSysPtr,BCsSystem *t_bcsSysPtr, LoadController *t_loadCtrlPtr);
    /**
     * check solution system's inition
    */
    PetscErrorCode checkInit();
    /**
     * set if show every iteration's information
     * @param t_ifShowIterInfo > true for show, false for close
    */
    PetscErrorCode showIterInfo(bool t_ifShowIterInfo);
    /**
     * run a single increment
     * @param t_ifLastConverged > if last increment converged
     * @param t_ifConverged > if this increment converged
     * @param t_ifCompleted > if the load factor reach the final load factor
    */
    PetscErrorCode run(bool t_ifLastConverged,bool *t_ifConverged, bool *t_ifcompleted);
public:
    bool m_ifShowIterInfo;          /**< if show every iteration's information*/
    SolutionCtx m_solutionCtx;      /**< solution context*/
    MoniterCtx m_monitorCtx;        /**< monitor context*/
    SNES m_snes;                    /**< snes solver*/
};
