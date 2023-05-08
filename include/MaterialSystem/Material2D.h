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
    Material2D(bool nLarge,double t_det_dx0dr):m_nLarge(nLarge),
                            m_F(Rank2Tensor2d::InitMethod::IDENTITY),
                            m_F0(Rank2Tensor2d::InitMethod::IDENTITY),
                            m_S(),m_J(1.0),m_det_dx0dr(t_det_dx0dr){}
    Material2D():m_nLarge(false),
                m_F(Rank2Tensor2d::InitMethod::IDENTITY),
                m_F0(Rank2Tensor2d::InitMethod::IDENTITY),
                m_S(),m_J(1.0),m_det_dx0dr(0.0){}
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
    double              m_det_dx0dr;/**< det(dx0/dr), x0 is elmt's coords in ref config*/
};