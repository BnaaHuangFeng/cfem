//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.10.17
//+++ Update : 2022.07.24->re-write & re-design @Yang Bai
//+++ Purpose: Implement rank-2 tensor class for the common
//+++          tensor manipulation in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/Rank2Tensor2d.h"
#include "Eigen/Eigen"

Rank2Tensor2d::Rank2Tensor2d():m_vals{}{}
Rank2Tensor2d::Rank2Tensor2d(const double val):m_vals{val,val,val,val}{}
Rank2Tensor2d::Rank2Tensor2d(const Rank2Tensor2d &a){
    for(int i=0;i<N2;i++) m_vals[i]=a.m_vals[i];
}
Rank2Tensor2d::Rank2Tensor2d(const InitMethod &initmethod){
    if(initmethod==InitMethod::ZERO){
        for(int i=0;i<N2;i++)m_vals[i]=0.0;
    }
    else if(initmethod==InitMethod::IDENTITY){
        setToIdentity();
    }
    else if(initmethod==InitMethod::RANDOM){
        setToRandom();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported initialize method in rank-2 tensor");
        MessagePrinter::exitcfem();
    }
}
Rank2Tensor2d::~Rank2Tensor2d(){
}
Rank2Tensor Rank2Tensor2d::toRank2Tensor3d()const{
    Rank2Tensor tmp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            tmp(i,j)=(*this)(i,j);
        }
    }
    return tmp;
}
//**********************************************************************
Rank2Tensor2d operator*(const double lhs,const Rank2Tensor2d &a){
    Rank2Tensor2d temp(0.0);
    for(int i=0;i<a.N2;i++) temp.m_vals[i]=lhs*a.m_vals[i];
    return temp;
}
Vector2d operator*(const Vector2d &lhs,const Rank2Tensor2d &a){
    Vector2d temp(0.0);
    for(int j=1;j<=a.N;j++){
        for(int i=1;i<=a.N;i++)
            temp(j)+=lhs(i)*a(i,j);
    }
    return temp;
}
Rank2Tensor2d Rank2Tensor2d::operator*(const ViogtRank2Tensor2D &a) const{
    // return A*B(still rank-2 tensor)
    Rank2Tensor2d temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            temp(i,j)=(*this)(i,1)*a(1,j)+(*this)(i,2)*a(2,j);
        }
    }
    return temp;
}
double Rank2Tensor2d::doubledot(const ViogtRank2Tensor2D &a) const{
    // return A:B calculation
    double sum=0.0;
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            // You may see A:B=A_ijB_ji in other books/literature, here we use A_ijB_ij
            // in Rank4Tensor, we follow the same definition!
            sum+=(*this)(i,j)*a(i,j);// use this to get the positive definite case!!!
        }
    }
    return sum;
}
//**************************************************************
//*** For eigen value and eigen vectors and other 
//*** stress and strain decomposition related functions
//**************************************************************
void Rank2Tensor2d::calcEigenValueAndEigenVectors(double (&eigval)[2],Rank2Tensor2d &eigvec) const{
    Eigen::Matrix2d _M;

    _M<<(*this)(1,1),(*this)(1,2),
        (*this)(2,1),(*this)(2,2);
    
    Eigen::EigenSolver<Eigen::Matrix2d> _eigen_solver;
    _eigen_solver.compute(_M);
    //
    eigval[0]=_eigen_solver.eigenvalues()(0).real();
    eigval[1]=_eigen_solver.eigenvalues()(1).real();
    for(int i=0;i<N;i++){
        eigvec(1,i+1)=_eigen_solver.eigenvectors()(0,i).real();
        eigvec(2,i+1)=_eigen_solver.eigenvectors()(1,i).real();
    }
}