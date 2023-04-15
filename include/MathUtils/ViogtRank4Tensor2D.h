//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : HuangFeng
//+++ Date   : 2023.04.11
//+++ Purpose: defines the tensor of 2D and rank2 with viogt symmetry
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# pragma once
# include "MathUtils/MatrixXd.h"
# include "MathUtils/ViogtRank2Tensor2D.h"
# include "MathUtils/Rank2Tensor2d.h"
# include "MathUtils/Rank4Tensor.h"
class ViogtRank2Tensor2D;
class Rank2Tensor2d;
class Rank4Tensor;
class ViogtRank4Tensor2D:public MatrixXd{
    public:
    /**
     * different initial method for rank-2 tensor
    */
    enum InitMethod{
        ZERO,
        IDENTITY,   // A_ijkl=delt_ik*delt_jl
        RANDOM
    };
    public:
    ViogtRank4Tensor2D():MatrixXd(3,3){};
    ViogtRank4Tensor2D(const double &val):MatrixXd(3,3,val){};
    ViogtRank4Tensor2D(InitMethod initmethod);
    ViogtRank4Tensor2D(Rank4Tensor rank4Tensor);
    /**
     * return the corresponding Rank4Tensor
    */
    Rank4Tensor toRank4Tensor();
    /**
     * @param indij viogt index sf 1
     * @param indkl viogt index sf 1
    */
    inline double & operator()(const int indij,const int indkl){
        return MatrixXd::operator()(indij,indkl);
    }
    /**
     * @param indij viogt index sf 1
     * @param indkl viogt index sf 1
    */
    inline double operator()(const int indij,const int indkl)const{
        return MatrixXd::operator()(indij,indkl);
    }   
    inline double & operator()(const int i,const int j,const int k,const int l){
        return MatrixXd::operator()(ij2ind[i-1][j-1],ij2ind[k-1][l-1]);
    }
    double operator()(const int i,const int j,const int k,const int l)const{
        return MatrixXd::operator()(ij2ind[i-1][j-1],ij2ind[k-1][l-1]);
    };
    /**
     * return tmp_ij=a_ijkl*b_kl
    */
    ViogtRank2Tensor2D operator*(const ViogtRank2Tensor2D &b)const;
    /**
     * return tmp_kl=b_ij*a_ijkl
    */
    friend ViogtRank2Tensor2D operator*(const ViogtRank2Tensor2D &b,const ViogtRank4Tensor2D &a);
    /**
     * return tmp_ij=L_ijkl*R_kl
    */
    ViogtRank2Tensor2D operator*(const Rank2Tensor2d &b)const;
    /**
     * return tmp_kl=L_ij*R_ijkl
    */
    friend ViogtRank2Tensor2D operator*(const Rank2Tensor2d &L,const ViogtRank4Tensor2D & R);
    /**
     * renturn tmp_ijkl=L_ijmn*R_mnkl
    */
    ViogtRank4Tensor2D operator*(const ViogtRank4Tensor2D &R);
    void print() const;
    public:
    static const int dim=2; /**< tensor dimension*/
    static const int NViogt=3; /**< number of viogt index*/
    constexpr static const int ij2ind[2][2]={{1,3},{3,2}};
    constexpr static const int ind2ij[3][2]={{1,1},{2,2},{1,2}};
};