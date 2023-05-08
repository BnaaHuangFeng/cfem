#pragma once
#include "petsc.h"
#include "MeshSystem/MeshSystem.h"
#include "InputSystem/InputSystem.h"
#include "ElementSystem/ElementSystem.h"
#include "BCsSystem/BCsSystem.h"
#include "LoadController/LoadController.h"
#include "SolutionSystem/SolutionSystem.h"
#include "PostProcessSystem/PostProcessSystem.h"
/**
 * Init Mesh System
 * @param meshDesPtr > ptr to mesh description
 * @param meshSysPtr > ptr to dynamic created mesh system
*/
PetscErrorCode MeshSystemInit(Timer *timerPtr, MeshDescription *meshDesPtr,MeshSystem **meshSysPtrAdr);
/**
 * Init element system
*/
PetscErrorCode ElmtSystemInit(Timer *timerPtr, ElementSystem **elmtSysPtrAdr, ElementDescription *elmtDesPtr,
                            MaterialDescription *matDesPtr, MeshSystem *meshSysPtr);
/**
 * Init boundary condition system
*/
PetscErrorCode BCsSystemInit(BCsSystem **BCsSysPtrAdr, BCDescription *BCDesPtr, MeshSystem *meshSysPtr);
/**
 * Init load controller
*/
PetscErrorCode LoadCtrolInit(LoadController **loadCtrlPtrAdr, StepDescriptiom *stepDesPtr, MeshSystem *meshSysPtr);
/**
 * Init solution system
*/
PetscErrorCode SolutionSysInit(SolutionSystem **solutionSysPtrAdr, StepDescriptiom *stepDesPtr, MeshSystem *meshSysPtr,
                            ElementSystem *elmtSysPtr, BCsSystem *BCsSysPtr, LoadController *loadCtrlPtr);
/**
 * Init postprocess system
*/
PetscErrorCode PostSysInit(PostProcessSystem **postSysPtrAdr,OutputDescription *t_outputDesPtr,MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr, LoadController *t_loadCtrlPtr);