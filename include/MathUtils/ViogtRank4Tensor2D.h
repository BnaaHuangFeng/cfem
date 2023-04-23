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
class ViogtRank2Tensor2D;
class Rank2Tensor2d;
class Rank4Tensor3d;
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
    ViogtRank4Tensor2D();
    ViogtRank4Tensor2D(const double &val);
    ViogtRank4Tensor2D(const double *vals);
    ViogtRank4Tensor2D(const MatrixXd &matrix);
    ViogtRank4Tensor2D(InitMethod initmethod);
    ViogtRank4Tensor2D(Rank4Tensor3d rank4Tensor);
    /**
     * return the corresponding Rank4Tensor3d
    */
    Rank4Tensor3d toRank4Tensor();
    /**
     * @param indij viogt index sf 1
     * @param indkl viogt index sf 1
    */
    double & operator()(const int indij,const int indkl);
    /**
     * @param indij viogt index sf 1
     * @param indkl viogt index sf 1
    */
    double operator()(const int indij,const int indkl)const;  
    double & operator()(const int i,const int j,const int k,const int l);
    double operator()(const int i,const int j,const int k,const int l)const;
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
    ViogtRank4Tensor2D operator+(ViogtRank4Tensor2D &R);
    ViogtRank4Tensor2D operator+(ViogtRank4Tensor2D R);
    ViogtRank4Tensor2D operator-(const ViogtRank4Tensor2D &R);
    /**
     * return it's full matrix form, index order: 11->1 21->2 12->3 22->4
    */
    MatrixXd toFullMatrix()const;
    void print() const;
    public:
    static const int dim=2; /**< tensor dimension*/
    static const int NViogt=3; /**< number of viogt index*/
    constexpr static const int ij2ind[2][2]={{0,2},{2,1}};
    constexpr static const int ind2ij[3][2]={{0,0},{1,1},{0,1}};
};