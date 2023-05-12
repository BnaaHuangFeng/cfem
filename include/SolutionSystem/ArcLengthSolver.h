#pragma once
#include "petsc.h"
#include "InputSystem/DescriptionInfo.h"
#include "MeshSystem/MeshSystem.h"
enum class PredictorMethod{
    STIFFNESSDET,
    INCREWORK,
    SECANTPATH
};
class ArcLengthSolver
{
private:
    AlgorithmType m_algorithm;
    int         m_maxIter;          /**< max iteration num limit*/      
    PetscScalar m_absTol;           /**< absolute tolerance*/
    PetscScalar m_relTol;           /**< relative tolerance*/
    PetscScalar m_uIncTol;          /**< incremental dof tolerance*/
    PetscScalar m_div_tol;          /**< the divergence tolerance. Use -1 to deactivate the test, default is 1e4*/
    PetscScalar m_arcLen;
    double      signum;             /**< stand for sign of arc len, 1.0 for positive, -1.0 for negative*/
    PetscScalar m_factorInc;        /**< incremental factor*/
    PetscScalar m_factorIter;        /**< incremental factor*/
    KSP         *m_kspPtr;          /**< ksp solver*/
    PC          *m_pcPtr;           /**< preconditioner*/
    SNESLineSearch *m_snesLinesearchPtr;/**< snes line search solver*/
    bool m_ifStepDesRead;           /**< if has read step description*/
    bool m_ifSolverInit;            /**< if the solver inited*/
    bool m_ifMeshSysSet;            /**< if the mesh system it relt on has been set*/
    bool m_ifFunctionSet;
    bool m_ifJacobianSet;
    bool m_ifLoadVecPtrSet;
    bool m_ifLastSolutionPtrSet;

    const int m_maxAttempt=10;      /**< max diverged increments tolerance*/
    int m_mDiverged;                /**< consecutive divergence num*/
    int m_mIter;                    /**< iteration num*/

    PetscScalar m_rnorm;            /**< 2 norm function value (may be estimated)*/
    PetscScalar m_duNorm;           /**< |du| of current iter*/
    PetscScalar m_rnorm0;           /**< 2 norm function value of iteration 0*/
    PetscScalar m_uNorm;            /**< 2 norm of solution*/
private:
    void readStepDes(StepDescriptiom *t_stepDesPtr);
    PetscErrorCode checkInit();
    /**
     * set mesh systemptr, including create global Vec
     * @param t_meshSysPtr > ptr to mesh system
    */
    void setMeshSystem(MeshSystem *t_meshSysPtr);
    /**
     * cal the sign of the incremental factor
     * @param t_uInc2 < ptr to incremental solution of last converged
     * @return > double val stand for sign of arc len, 1.0 for positive, -1.0 for negative
    */
    PetscErrorCode getSigNum(double *t_sigNum);
    PetscErrorCode updateFactor(double sigNum);
    PetscErrorCode updateConvergenceReason();
public:
    ArcLengthSolver();
    ArcLengthSolver(StepDescriptiom *t_stepDesPtr);
    ~ArcLengthSolver();
    PetscErrorCode init(StepDescriptiom *t_stepDesPtr,MeshSystem *t_meshSysPtr, KSP *t_kspPtr, PC *t_pcPtr, SNESLineSearch *t_snesLinesearchPtr);
    /**
     * set residual ptr (the Vec need preallocation), and residual cal callback function ptr
     * @param rPtr > ptr to residual (need preallocation, dont't destroy the Vec until this solver isn's needed)
    */
    PetscErrorCode setFunction(Vec *rPtr, PetscErrorCode (*f)(ArcLengthSolver *t_solverPtr, Vec *t_uInc, Vec *t_function, void *ctx),void *t_ctx);
    PetscErrorCode setLoad(PetscErrorCode (*t_applyPresetDof)(Mat *AMatrixPtr, Vec *t_b, void *ctx),void *t_ctx){
        m_applyLoad=t_applyPresetDof;
        m_LoadCtx=t_ctx;
        return 0;
    }
    /**
     * set  ptr to Jacobian matrix, ptr to material for precondition of jacobian and residual cal callback function ptr
     * @param AmatPtr > ptr to Jacobian matrix (need preallocation, dont't destroy the Mat until this solver isn's needed)
     * @param PmatPtr > ptr to material for precondition of jacobian (need preallocation, dont't destroy the Mat until this solver isn's needed)
    */    
    PetscErrorCode setJacobian(Mat *AmatPtr, Mat *PmatPtr, PetscErrorCode (*J)(ArcLengthSolver *t_solverPtr, Vec *t_uInc, Mat *t_AMat, Mat *t_PMat, void *ctx),void *t_ctx);
    PetscErrorCode monitorSet(PetscErrorCode (*t_monitor)(ArcLengthSolver *t_solverPtr, PetscInt its, PetscScalar norm, void *mctx),void *t_ctx){
        m_monitor=t_monitor;
        m_monitorCtx=t_ctx;
        return 0;
    }
    /**
     * Set m_uInc2Ptr ponit to solution of the last converged increment
     * @param t_solution > ptr to last converged solution
    */
    void setLastSolutionPtr(Vec *t_solutionPtr){m_uInc2Ptr=t_solutionPtr; m_ifLastSolutionPtrSet=true;}
    /**
     * Set m_loadPtr point to load Vec
    */
    void setLoadVecPtr(Vec *t_loadVecPtr){m_loadPtr=t_loadVecPtr; m_ifLoadVecPtrSet=true;}
    /**
     * start solution iteration
     * @param uIncInitialPtr > ptr to Vec used for initial solution (need preallocation, can be destroyed after the call)
    */
    PetscErrorCode solve(double t_arcLen);
    double getIntialArcLen(double factorInc0);
    PetscErrorCode calTangentIter(KSPConvergedReason *reasonPtr);
    PetscErrorCode calResidualIter(KSPConvergedReason *reasonPtr);

    PetscErrorCode getConvergedReason(SNESConvergedReason *t_converReasonPtr);
    /**
     * get solution Vec (need preallocation, dont't destroy the Vec until this solver isn's needed)
     * @param t_solution > the Vec (need preallocation)
    */
    PetscErrorCode getSolution(Vec *t_solution);
    inline double getFactorInc(){return m_factorInc;};
    inline int getIterNum(){return m_mIter;};
public:
    StepDescriptiom         *m_stepDesPtr;      /**< ptr to step description*/
    Vec                     *m_rPtr;            /**< ptr to residual (need preallocation, dont't destroy the Vec until this solver isn's needed)*/
    Vec                     *m_loadPtr;         /**< ptr to global load Vec*/
    Vec                     m_uIterRes;         /**< residual iterational solution*/
    Vec                     m_uIterTangent;     /**< tangent iterational solution*/
    Vec                     m_uInc;             /**< incremental solution*/
    Vec                     m_b;                /**< RHS of equation about m_uIterTangent*/
    Vec                     *m_uInc2Ptr;        /**< ptr to last converged solution (used for load direction predictor)*/
    Mat                     *m_AMatPtr;         /**< ptr to Jacobian matrix (need preallocation, dont't destroy the Mat until this solver isn's needed)*/
    Mat                     *m_PMatPtr;         /**< ptr to material for precondition of jacobian (need preallocation, dont't destroy the Mat until this solver isn's needed)*/
    SNESConvergedReason     m_converReason;
    MeshSystem              *m_meshSysPtr;      /**< ptr to mesh system it relies on*/
    void                    *m_functionCtx;
    void                    *m_jacobianCtx;
    void                    *m_LoadCtx;
    void                    *m_monitorCtx;
    PetscErrorCode (*m_functionCal)(ArcLengthSolver *t_solverPtr, Vec *t_uInc, Vec *t_function, void *ctx);
    PetscErrorCode (*m_jacobianCal)(ArcLengthSolver *t_solverPtr, Vec *t_uInc, Mat *t_AMat, Mat *t_PMat, void *ctx);
    PetscErrorCode (*m_applyLoad)(Mat *AMatrixPtr, Vec *t_b, void *ctx);
    PetscErrorCode (*m_monitor)(ArcLengthSolver *t_solverPtr, PetscInt its, PetscScalar norm, void *mctx);
};
