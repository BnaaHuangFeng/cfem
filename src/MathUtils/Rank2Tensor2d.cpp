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
#include "MathUtils/Vector2d.h"
#include "MathUtils/Rank2Tensor3d.h"
#include "MathUtils/Rank2Tensor2d.h"
#include "Eigen/Eigen"
#include <algorithm>
Rank2Tensor2d::Rank2Tensor2d():m_vals{}{
    Rank2Tensor2d(Rank2Tensor2d::InitMethod::ZERO);
}
Rank2Tensor2d::Rank2Tensor2d(const Rank2Tensor2d &a){
    for(int i=0;i<N2;i++) m_vals[i]=a.m_vals[i];
}
Rank2Tensor2d::Rank2Tensor2d(const ViogtRank2Tensor2D &a){
    for(int mI=0;mI<N;mI++){
        for(int nI=0;nI<N;nI++){
            (*this)(mI,nI)=a(mI,nI);
        }
    }
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
Rank2Tensor3d Rank2Tensor2d::toRank2Tensor3d()const{
    Rank2Tensor3d tmp(Rank2Tensor3d::InitMethod::ZERO);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            tmp(i,j)=(*this)(i,j);
        }
    }
    return tmp;
}
//**********************************************************************
Rank2Tensor2d operator*(const double lhs,const Rank2Tensor2d &a){
    Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
    for(int i=0;i<a.N2;i++) temp.m_vals[i]=lhs*a.m_vals[i];
    return temp;
}
Vector2d operator*(const Vector2d &lhs,const Rank2Tensor2d &a){
    Vector2d temp(Vector2d::InitMethod::ZERO);
    for(int j=0;j<a.N;j++){
        for(int i=0;i<a.N;i++)
            temp(j)+=lhs(i)*a(i,j);
    }
    return temp;
}
Rank2Tensor2d Rank2Tensor2d::operator*(const ViogtRank2Tensor2D &a) const{
    // return A*B(still rank-2 tensor)
    Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            temp(i,j)=(*this)(i,0)*a(0,j)+(*this)(i,1)*a(1,j);
        }
    }
    return temp;
}
double Rank2Tensor2d::doubledot(const ViogtRank2Tensor2D &a) const{
    // return A:B calculation
    double sum=0.0;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            // You may see A:B=A_ijB_ji in other books/literature, here we use A_ijB_ij
            // in Rank4Tensor3d, we follow the same definition!
            sum+=(*this)(i,j)*a(i,j);// use this to get the positive definite case!!!
        }
    }
    return sum;
}
//**************************************************************
//*** For eigen value and eigen vectors and other 
//*** stress and strain decomposition related functions
//**************************************************************
void Rank2Tensor2d::calcEigenValueAndEigenVectors(double eigvalPtr[2],Vector2d eigvecPtr[2]) const{
    Eigen::Matrix2d _M;

    _M<<(*this)(0,0),(*this)(0,1),
        (*this)(1,0),(*this)(1,1);
    Eigen::EigenSolver<Eigen::Matrix2d> m_eigen_solver;    /**< solver for eigen solve*/
    m_eigen_solver.compute(_M);
    //
    eigvalPtr[0]=m_eigen_solver.eigenvalues()(0).real();
    eigvalPtr[1]=m_eigen_solver.eigenvalues()(1).real();
    for(int i=0;i<N;i++){
        eigvecPtr[i](0)=m_eigen_solver.eigenvectors()(0,i).real();
        eigvecPtr[i](1)=m_eigen_solver.eigenvectors()(1,i).real();
    }
}
void Rank2Tensor2d::spectralDecomposition(double eigvalPtr[2],ViogtRank2Tensor2D eigprojPtr[2],bool &repeated)const{
    static const double small=1E-5;
    Vector2d eigenvec[2];
    calcEigenValueAndEigenVectors(eigvalPtr,eigenvec);
    double differ=abs(eigvalPtr[0]-eigvalPtr[1]);
    double maxEigen=std::max(abs(eigvalPtr[0]),abs(eigvalPtr[1]));
    if(maxEigen!=0.0)
        differ=differ/maxEigen;
    repeated=differ<small;
    for(int eigI=0;eigI<N;eigI++){
        eigprojPtr[eigI](0)=eigenvec[eigI](0)*eigenvec[eigI](0);
        eigprojPtr[eigI](1)=eigenvec[eigI](1)*eigenvec[eigI](1);
        eigprojPtr[eigI](2)=eigenvec[eigI](0)*eigenvec[eigI](1);
    }
}

Vector2d Rank2Tensor2d::getIthRow(const int i)const{
    Vector2d temp(Vector2d::InitMethod::ZERO);
    temp(0)=(*this)(i,0);
    temp(1)=(*this)(i,1);
    return temp;
}
Vector2d Rank2Tensor2d::getIthCol(const int i)const{
    Vector2d temp(Vector2d::InitMethod::ZERO);
    temp(0)=(*this)(0,i);
    temp(1)=(*this)(1,i);
    return temp;
}
Vector2d Rank2Tensor2d::operator*(const Vector2d &a) const{
    Vector2d temp(Vector2d::InitMethod::ZERO);
    for(int i=0;i<N;i++){
        temp(i)=(*this)(i,0)*a(0)+(*this)(i,1)*a(1);
    }
    return temp;
}
void Rank2Tensor2d::setFromVectorDyad(const Vector2d &a,const Vector2d &b){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            (*this)(i,j)=a(i)*b(j);
        }
    }
}