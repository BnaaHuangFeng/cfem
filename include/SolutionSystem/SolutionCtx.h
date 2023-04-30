#pragma once
#include "MeshSystem/MeshSystem.h"
#include "ElementSystem/ElementSystem.h"
#include "BCsSystem/BCsSystem.h"
#include "LoadController/LoadController.h"
/**
 * context structure that nolinear solver need
*/
struct SolutionCtx
{
    ElementSystem *s_elmtSysPtr;      /**< ptr to the elmt system it relied on*/
    BCsSystem *s_bcsSysPtr;           /**< ptr to the boundary conditon system it relied on*/
    LoadController *s_loadCtrlPtr;    /**< ptr to the load controller it relied on*/
};
