#pragma once
class Rank2Tensor2d;
class ViogtRank2Tensor2D;
class ViogtRank4Tensor2D;
namespace TensorConst2D{
    /*****************************************
     * declare of const tensor            ****
    *****************************************/
    extern const double _I_d[3];
    // const double unitRank4[9]={
    //     1.0,    0.0,    0.0,
    //     0.0,    1.0,    0.0,
    //     0.0,    0.0,    1.0
    // };
    extern const double _IXI_d[9];
    extern const double _IISym_d[9];
    extern ViogtRank2Tensor2D I;       /**< 2 rank unit tensor del_ij*/
    // ViogtRank4Tensor2D II=ViogtRank4Tensor2D(unitRank4);      /**< 4 rank unit tensor del_ik*del_jl*/
    // ViogtRank4Tensor2D IIT;     /**< transpositional 4 rank unit tensor del_il*del_jk*/
    extern ViogtRank4Tensor2D IXI;     /**< spherical 4 rank unit tensor del_ij*del_kl*/
    extern ViogtRank4Tensor2D IISym;   /**< symmetric 4 rank unit tensor 1/2*(II+IIT)*/
    // ViogtRank4Tensor2D IISkew;  /**< skew-symmetric 4 rank tensor 1/2*(II-IIT)*/
    extern ViogtRank4Tensor2D IIVol;   /**< volumetric 4 rank unit tensor 1/3*IXI*/
    extern ViogtRank4Tensor2D IIDev;   /**< deviatoric 4 rank unit tensor IISym-IIVol*/
}