#pragma once
#include "petsc.h"
#include <vector>
#include "InputSystem/DescriptionInfo.h"
#include "MeshSystem/MeshSystem.h"
#include "MathUtils/Vector2d.h"
#include "Utils/Timer.h"
class Vector2d;
/**
 * this class store the topnology structure of 2D mesh,including node's ID. coords, element's connectivity
 */
using namespace std;
class MeshDescription;
class MeshSystem;
class MeshSystem2D:public MeshSystem{
public:
   inline MeshSystem2D(){m_dim=2;m_mDof_node=2;}
   inline MeshSystem2D(Timer *timerPtr):MeshSystem(timerPtr){m_dim=2;m_mDof_node=2;};
};