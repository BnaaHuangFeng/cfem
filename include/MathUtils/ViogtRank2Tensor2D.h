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
#include"Eigen/Eigen"
#include <fstream>
#include <iostream>
#include "Utils/MessagePrinter.h"
#include "MathUtils/Vector3d.h"
#include <cmath>
class Vector2d;
class Rank2Tensor3d;
class Rank2Tensor2d;
class ViogtRank4Tensor2D;
class MatrixXd;
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
    ViogtRank2Tensor2D(){ViogtRank2Tensor2D(ViogtRank2Tensor2D::InitMethod::ZERO);};
    /**
     * @param a the right hand vector
    */
    ViogtRank2Tensor2D(const double vals[3]):Vector3d(vals){}
    ViogtRank2Tensor2D(Vector3d valvec):Vector3d(valvec){}
    ViogtRank2Tensor2D(InitMethod initmethod);
    void setFromRank2Tensor2D(const Rank2Tensor2d &R);
    ~ViogtRank2Tensor2D(){}
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
    ViogtRank2Tensor2D operator*(const double a){
        ViogtRank2Tensor2D tmp;
        tmp=ViogtRank2Tensor2D(Vector3d::operator*(a));
        return tmp;
    }
    friend ViogtRank2Tensor2D operator*(const double L,ViogtRank2Tensor2D &R){
        return R*L;
    }
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
    friend Vector2d operator*(const Vector2d &a,const ViogtRank2Tensor2D &b);
    /**
     * operation /
    */
    ViogtRank2Tensor2D operator/(const double R);
    inline double trace(){return (*this)(0)+(*this)(1);}
    /**
     * operation %: tensor product:a_ijkl=L_ij*R_Kl
    */
    ViogtRank4Tensor2D operator%(const ViogtRank2Tensor2D &R);
    /**
     * * operator return L_ik*R_kj
    */
    Rank2Tensor2d operator*(const Rank2Tensor2d &R);
    Rank2Tensor2d toRank2Tensor2d()const;
    Rank2Tensor3d toRank2Tensor3d()const;
    inline void print() const{
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ***\n",(*this)(0,0),(*this)(0,1));
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ***\n",(*this)(1,0),(*this)(1,1));
    };
    //**********************************************************************
    /**
     * calculate the eigen value and eigen vector for current rank-2 tensor
     * @param eigvalPtr the double array, which stores the eigen value
     * @param eigvecPtr the related eigen vector2d
     */
    void calcEigenValueAndEigenVectors(double eigvalPtr[2],Vector2d eigvecPtr[2]);
    //**********************************************************************
    //*** for decomposition
    //**********************************************************************
    /**
     * spectral decomposition for symmtry matrix
     * @param eigvalPtr < ref to eigen values
     * @param eigProj < ref to corresponding eigenvalues projection
     * @param repeated < if the eigen values repeated
    */
    void spectralDecomposition(double eigvalPtr[2],ViogtRank2Tensor2D eigProjPtr[2],bool &repeated);
    /**
     * isotropic tensor function calculation for this tensor
     * @param func > ptr to conrresponding scalar function
    */
    ViogtRank2Tensor2D isotropicFunc(double (*func)(double));
    /**
     * refer to Box A.3 of CMFP, Cal isotropic symmetric rank 2 tensor function's derivate to symmetric rank 2 tensor for this tensor
     * @param func > ptr to conrresponding isotropic scalar function
     * @param funcDeriv > ptr to conrresponding isotropic scalar function's derivate (also dfunc(t)/dt)
    */
    ViogtRank4Tensor2D iostropicFuncDeriv(double (*func)(double),double (*funcDeriv)(double));
    /**
     * return rank 4 tensor's matrix form M_IJ=L_il*R_jk, ordering 11 21 12 22
    */
    MatrixXd iljk(const ViogtRank2Tensor2D &R)const;
    /**
     * return rank 4 tensor's matrix form M_IJ=L_ik*R_jl, ordering 11 21 12 22
    */
    MatrixXd ikjl(const ViogtRank2Tensor2D &R)const;
    private:
    constexpr static const int matInd2ij[4][2]={{0,0},{1,0},{0,1},{1,1}}; /**< Rank 4 tensor matrix form index matInd -> (i,j)*/
    constexpr static const int matInd2viogtInd[4]={0,2,2,1};
    constexpr static const int ij2ind[2][2]={{0,2},{2,1}};
    constexpr static const int ind2ij[3][2]={{0,0},{1,1},{0,1}};
    public:
    static const int dim=2; /**< tensor dimension*/
    static const int NViogt=3; /**< number of viogt index*/

};