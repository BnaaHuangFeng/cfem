#pragma once
#include "MaterialSystem/Material.h"
#include "MathUtils/ViogtRank2Tensor2D.h"
#include "MathUtils/ViogtRank4Tensor2D.h"
#include "MathUtils/Rank2Tensor2d.h"
#include "MathUtils/TensorConst2D.h"
#include "InputSystem/EnumDataType.h"
#include "nlohmann/json.hpp"
class ViogtRank2Tensor2D;
class ViogtRank4Tensor2D;
class Rank2Tensor2d;
class Material2D:public Material{
public:
    /**
    * construction
    */
    Material2D(bool nLarge):m_nLarge(nLarge),
                            m_F(TensorConst2D::I),m_F0(TensorConst2D::I),m_S(0.0),m_J(1.0){}
    Material2D():m_nLarge(false),
                m_F(TensorConst2D::I),m_F0(TensorConst2D::I),m_S(0.0),m_J(1.0){}
    /**
    * destruction
    */
    virtual ~Material2D(){};
public:
    bool                m_nLarge;   /**< true for large strain,false for small strain*/
public:
    Rank2Tensor2d       m_F;        /**< deformation Tensor F*/
    Rank2Tensor2d       m_F0;       /**< deformation Tensor F of last converged*/
    ViogtRank2Tensor2D  m_S;        /**< cauchy stress S*/
    double              m_J;        /**< det(m_F)*/

};