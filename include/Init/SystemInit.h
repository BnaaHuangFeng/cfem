#pragma once
#include"petsc.h"
#include"MeshSystem/MeshSystem.h"
#include"InputSystem/InputSystem.h"
/**
 * Init Mesh System
 * @param MeshDesPtr > ptr to mesh description
 * @param meshSysPtr > ptr to dynamic created mesh system
*/
PetscErrorCode MeshSystemInit(Timer *timerPtr, MeshDescription *MeshDesPtr,MeshSystem *meshSysPtr);
