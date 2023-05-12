#include "MathUtils/PetsExtension.h"
#include "petscmat.h"
PetscErrorCode PetscExtension::MatDet(Mat *matPtr, PetscScalar *detPtr, PetscInt rank){
    KSP         ksp; /* linear solver context */
    PC          pc;    
    PetscInt  icntl;
    Mat F;
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE));
    PetscCall(MatView(*matPtr,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, *matPtr, *matPtr));
    PetscCall(KSPSetType(ksp, KSPPREONLY));
    PetscCall(KSPGetPC(ksp, &pc));    
    PetscCall(PCSetType(pc, PCLU));
    PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    PetscCall(PCFactorSetUpMatSolverType(pc)); 
    PetscCall(PCFactorGetMatrix(pc, &F));
    PetscCall(MatMumpsSetIcntl(F, 33, 1));
    PetscCall(PCSetType(pc, PCLU));
    PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERPETSC));
    PetscCall(PCFactorSetUpMatSolverType(pc)); /* call MatGetFactor() to create F */
    PetscCall(PCFactorGetMatrix(pc, &F));
    PetscCall(KSPSetUp(ksp));
    PetscInt  infog34;
    PetscReal cntl, rinfo12, rinfo13;
    if (rank == 0){
        icntl = 3;
        PetscCall(MatMumpsGetCntl(F, icntl, &cntl));
        PetscCall(MatMumpsGetInfog(F, 34, &infog34));
        PetscCall(MatMumpsGetRinfog(F, 12, &rinfo12));
        PetscCall(MatMumpsGetRinfog(F, 13, &rinfo13));
        // MatMumpsGetNullPivots(F, &num_null_pivots, &null_pivots);
        PetscCall(PetscPrintf(PETSC_COMM_SELF, "  Mumps row pivot threshold = %g\n", cntl));
        PetscCall(PetscPrintf(PETSC_COMM_SELF, "  Mumps determinant = (%g, %g) * 2^%" PetscInt_FMT " \n", (double)rinfo12, (double)rinfo13, infog34));
    }
    *detPtr=rinfo12;
    return 0;
}
int PetscExtension::solQua(PetscScalar a, PetscScalar b, PetscScalar c, PetscScalar *r1Ptr, PetscScalar *r2Ptr){
    const PetscScalar small=1.0e-20;
    PetscScalar sigB=0.0, b2, r4AC, squar, q;
    if(a!=0.0){
        if(b!=0.0){
            if(b!=0.0){
                sigB=b/abs(b);
            }
            else{
                sigB=1.0;
            }
            b2=b*b;
            r4AC=4.0*a*c;
            squar=b2-r4AC;
            if(squar>small){
                squar=sqrt(squar);
                q=-(b+sigB*squar)/2.0;
                *r1Ptr=q/a;
                *r2Ptr=c/q;
                return 2;
            }
            else if(squar>=0){
                *r1Ptr=-b/(2.0*a);
                return 1;
            }
        }
    }
    else{
        if(b!=0.0){
            *r1Ptr=-c/b;
            *r2Ptr=*r1Ptr;
            return 1;
        }
    }
    return 0;
}