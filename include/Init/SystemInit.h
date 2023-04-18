#pragma once
#include"petsc.h"
#include"MeshSystem/StructuredMesh2D.h"
#include"InputSystem/InputSystem.h"
#include"InputSystem/DescriptionInfo.h"
/**
 * Init Mesh System
 * @param 
*/
PetscErrorCode MeshSystemInit(MeshDescription *MeshDesPtr,MeshSystem *meshSystemPtr);
