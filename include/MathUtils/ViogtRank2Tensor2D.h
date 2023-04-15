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
#pragma once
#include"Vector3d.h"
#include"Vector2d.h"
#include"Rank2Tensor.h"
#include"Rank2Tensor2d.h"
#include"MathUtils/ViogtRank4Tensor2D.h"
#include <fstream>
#include <iostream>
class Rank2Tensor;
class ViogtRank2Tensor2D:public Vector3d{
    public:
    /**
     * different initial method for rank-2 tensor
    */
    enum InitMethod{
        ZERO,
        IDENTITY,
        RANDOM
    };
    public:
        /**
         * constructor
        */
        ViogtRank2Tensor2D():Vector3d(){}
        /**
         * @param val the scalar value
        */
        ViogtRank2Tensor2D(const double val):Vector3d(val){};
        /**
         * @param a the right hand vector
        */
        ViogtRank2Tensor2D(const ViogtRank2Tensor2D &a):Vector3d((Vector3d &)a){};
        ViogtRank2Tensor2D(InitMethod initmethod);
        /**
         * () operator
         * @param i viogt index sf 1
         */
        inline double& operator()(const int i){
            return Vector3d::operator()(i);
        }  
        /**
         * () operator
         * @param i viogt index sf 1
         */
        inline double operator()(const int i)const{
            return Vector3d::operator()(i);
        }       
        /**
         * @param i the first index sf 1
         * @param j the second index sf 1
        */
        double& operator()(const int i,const int j);
        /**
         * @param i the first index sf 1
         * @param j the second index sf 1
        */
        double operator()(const int i,const int j)const;
        /**
         * * operator
         * @param a right hand side Vector2d value
         */
        Vector2d operator*(const Vector2d &a)const;
         /**
         * * operator
         * @param a left hand side Vector2d value
         * @param b right hand side ViogtRank2Tensor2D value
         */
        friend Vector2d operator*(const Vector2d &a,const ViogtRank2Tensor2D &b){
            return b*a;
        }
        /**
         * * operator return L_ik*R_kj
        */
        Rank2Tensor2d operator*(const Rank2Tensor2d &R);
        Rank2Tensor2d toRank2Tensor2d()const;
        Rank2Tensor toRank2Tensor3d()const;
        inline void print() const{
            PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ***\n",(*this)(1,1),(*this)(1,2));
            PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ***\n",(*this)(2,1),(*this)(2,2));
        };
    public:
    static const int dim=2; /**< tensor dimension*/
    static const int NViogt=3; /**< number of viogt index*/
    constexpr static const int ij2ind[2][2]={{1,3},{3,2}};
    constexpr static const int ind2ij[3][2]={{1,1},{2,2},{1,2}};
};