#include "MathUtils/TensorConst2D.h"
#include "MathUtils/Rank2Tensor2d.h"
#include "MathUtils/ViogtRank2Tensor2D.h"
#include "MathUtils/ViogtRank4Tensor2D.h"
using namespace TensorConst2D;
const double TensorConst2D::_I_d[3]={1.0,1.0,0.0};
const double TensorConst2D::_IXI_d[9]={
    1.0,    1.0,    0.0,
    1.0,    1.0,    0.0,
    0.0,    0.0,    0.0
};
const double TensorConst2D::_IISym_d[9]={
    1.0,    0.0,    0.0,
    0.0,    1.0,    0.0,
    0.0,    0.0,    0.5       
};
ViogtRank2Tensor2D TensorConst2D::I=ViogtRank2Tensor2D(_I_d);       /**< 2 rank unit tensor del_ij*/
// ViogtRank4Tensor2D II=ViogtRank4Tensor2D(unitRank4);      /**< 4 rank unit tensor del_ik*del_jl*/
// ViogtRank4Tensor2D IIT;     /**< transpositional 4 rank unit tensor del_il*del_jk*/
ViogtRank4Tensor2D TensorConst2D::IXI=ViogtRank4Tensor2D(_IXI_d);     /**< spherical 4 rank unit tensor del_ij*del_kl*/
ViogtRank4Tensor2D TensorConst2D::IISym=ViogtRank4Tensor2D(_IISym_d);   /**< symmetric 4 rank unit tensor 1/2*(II+IIT)*/
// ViogtRank4Tensor2D IISkew;  /**< skew-symmetric 4 rank tensor 1/2*(II-IIT)*/
ViogtRank4Tensor2D TensorConst2D::IIVol=IXI*(1.0/3.0);   /**< volumetric 4 rank unit tensor 1/3*IXI*/
ViogtRank4Tensor2D TensorConst2D::IIDev=IISym-IIVol;   /**< deviatoric 4 rank unit tensor IISym-IIVol*/