#pragma once
#include "petsc.h"
namespace PetscExtension{
/**
 * cal a matrxi's determinant
 * @param matPtr > ptr to the matrix to calculate
 * @param detPtr > address to receive the determinant result 
*/
PetscErrorCode MatDet(Mat *matPtr, PetscScalar *detPtr, PetscInt rank); 
int solQua(PetscScalar a, PetscScalar b, PetscScalar c, PetscScalar *r1Ptr, PetscScalar *r2Ptr);   
};
