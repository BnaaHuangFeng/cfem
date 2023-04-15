#include <iostream>
#include "petsc.h"

int main(int args,char *argv[]){
    PetscErrorCode ierr;
    ierr=PetscInitialize(&args,&argv,NULL,NULL);if (ierr) return ierr;
   

    ierr=PetscFinalize();CHKERRQ(ierr);
    return ierr;
}
