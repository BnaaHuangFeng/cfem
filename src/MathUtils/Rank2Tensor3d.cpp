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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "MathUtils/Rank2Tensor3d.h"
#include "Eigen/Eigen"
#include "MathUtils/Rank4Tensor3d.h"
Rank2Tensor3d::Rank2Tensor3d(){
    m_vals.resize(9,0.0);
}
Rank2Tensor3d::Rank2Tensor3d(const Rank2Tensor3d &a){
    m_vals.resize(9,0.0);
    for(int i=0;i<N2;i++) m_vals[i]=a.m_vals[i];
}
Rank2Tensor3d::Rank2Tensor3d(const InitMethod &initmethod){
    if(initmethod==InitMethod::ZERO){
        m_vals.resize(9,0.0);
    }
    else if(initmethod==InitMethod::IDENTITY){
        m_vals.resize(9,0.0);
        setToIdentity();
    }
    else if(initmethod==InitMethod::RANDOM){
        m_vals.resize(9,0.0);
        setToRandom();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported initialize method in rank-2 tensor");
        MessagePrinter::exitcfem();
    }
}
Rank2Tensor3d::~Rank2Tensor3d(){
    m_vals.clear();
}
//**********************************************************************
Rank2Tensor3d operator*(const double &lhs,const Rank2Tensor3d &a){
    Rank2Tensor3d temp(Rank2Tensor3d::InitMethod::ZERO);
    for(int i=0;i<a.N2;i++) temp.m_vals[i]=lhs*a.m_vals[i];
    return temp;
}
Vector3d operator*(const Vector3d &lhs,const Rank2Tensor3d &a){
    Vector3d temp(Vector3d::InitMethod::ZERO);
    for(int j=0;j<a.N;j++){
        temp(j)=lhs(0)*a(0,j)+lhs(1)*a(1,j)+lhs(2)*a(2,j);
    }
    return temp;
}
Rank2Tensor3d Rank2Tensor3d::doubledot(const Rank4Tensor3d &a) const{
        // return A:B calculation
        Rank2Tensor3d temp(Rank2Tensor3d::InitMethod::ZERO);
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                for(int k=0;k<N;k++){
                    for(int l=0;l<N;l++){
                        temp(k,l)+=(*this)(i,j)*a(i,j,k,l);
                    }
                }
            }
        }
        return temp;
    }
//*******************************************************************
//*** for advanced math operators
//*******************************************************************
Rank2Tensor3d exp(const Rank2Tensor3d &a){
    Rank2Tensor3d I;
    I.setToIdentity();
    return I
          +a
          +a*a*(1.0/(1.0*2.0))
          +a*a*a*(1.0/(1.0*2.0*3.0))
          +a*a*a*a*(1.0/(1.0*2.0*3.0*4.0))
          +a*a*a*a*a*(1.0/(1.0*2.0*3.0*4.0*5.0))
          +a*a*a*a*a*a*(1.0/(1.0*2.0*3.0*4.0*5.0*6.0));
}
Rank2Tensor3d dexp(const double &a,const Rank2Tensor3d &b){
    // return dexp(ab)/db
    return b
          +b*b*a*(1.0/1.0)
          +b*b*b*a*a*(1.0/(1.0*2.0))
          +b*b*b*b*a*a*a*(1.0/(1.0*2.0*3.0))
          +b*b*b*b*b*a*a*a*a*(1.0/(1.0*2.0*3.0*4.0))
          +b*b*b*b+b*b*a*a*a*a*a*(1.0/(1.0*2.0*3.0*4.0*5.0));
}
//**************************************************************
//*** For eigen value and eigen vectors and other 
//*** stress and strain decomposition related functions
//**************************************************************
void Rank2Tensor3d::calcEigenValueAndEigenVectors(double (&eigval)[3],Rank2Tensor3d &eigvec) const{
    Eigen::Matrix3d _M;

    _M<<(*this)(0,0),(*this)(0,1),(*this)(0,2),
        (*this)(1,0),(*this)(1,1),(*this)(1,2),
        (*this)(2,0),(*this)(2,1),(*this)(2,2);
    
    Eigen::EigenSolver<Eigen::Matrix3d> _eigen_solver;
    _eigen_solver.compute(_M);
    //
    eigval[0]=_eigen_solver.eigenvalues()(0).real();
    eigval[1]=_eigen_solver.eigenvalues()(1).real();
    eigval[2]=_eigen_solver.eigenvalues()(2).real();
    for(int i=0;i<N;i++){
        eigvec(0,i)=_eigen_solver.eigenvectors()(0,i).real();
        eigvec(1,i)=_eigen_solver.eigenvectors()(1,i).real();
        eigvec(2,i)=_eigen_solver.eigenvectors()(2,i).real();
    }
}