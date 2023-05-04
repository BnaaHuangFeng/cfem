//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author   : Yang Bai
//+++ Date     : 2020.10.18
//+++ Reviewer : Xiaoyuan @ 2021.08.20
//+++ Purpose  : Define the general Matrix  in AsFem
//+++            we mainly use this for the calculation of jacobian
//+++            If one wants to use Eigen's MatrixXd, please use
//+++            Eigen::MatrixXd, which is different with ours !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "MathUtils/MatrixXd.h"
#include "MathUtils/VectorXd.h"
#include "MathUtils/ViogtRank2Tensor2D.h"
MatrixXd::MatrixXd(){
    m_m=m_n=m_mn=0;
    m_vals.clear();
}
MatrixXd::MatrixXd(const MatrixXd &a){
    m_vals=a.m_vals;
    m_m=a.m_m;m_n=a.m_n;
    m_mn=m_m*m_n;
}
MatrixXd::MatrixXd(const int  m,const int  n){
    m_m=m;m_n=n;
    m_mn=m*n;
    m_vals.resize(m_mn,0.0);
}
MatrixXd::MatrixXd(const int  m,const int  n,const double val){
    m_m=m;m_n=n;
    m_mn=m*n;
    m_vals.resize(m_mn,val);
}
MatrixXd::MatrixXd(const int  m,const int  n,const double *vals){
    m_m=m;m_n=n;
    m_mn=m*n;
    m_vals.resize(m_mn,0.0);
    for(int rowI=0;rowI<m_m;rowI++){// loop over row
        for(int colI=0;colI<m_n;colI++){// loop over col
            m_vals[rowI*m_n+colI]=*vals++;
        }
    }
}
void MatrixXd::solve(const VectorXd &b,VectorXd &x) const{
    if(b.getM()!=getM()){
        MessagePrinter::printErrorTxt("size of rhs vector b is not equal to the row number of your matrix, can\'t execute the solve function");
        MessagePrinter::exitcfem();
    }
    if(x.getM()!=getN()){
        MessagePrinter::printErrorTxt("size of solution vector x is not equal to the column number of your matrix, can\'t execute the solve function");
        MessagePrinter::exitcfem();
    }
    Eigen::MatrixXd A(getM(),getN());
    Eigen::VectorXd B(getM()),X(getN());
    for(int i=0;i<getM();i++){
        for(int j=0;j<getN();j++){
            A.coeffRef(i,j)=(*this)(i,j);
        }
        B.coeffRef(i)=b(i);
    }
    X=A.fullPivLu().solve(B);
    for(int i=0;i<getM();i++){
        x(i)=X.coeff(i);
    }
}
VectorXd MatrixXd::solve(const VectorXd &b) const{
    if(b.getM()!=getM()){
        MessagePrinter::printErrorTxt("size of rhs vector b is not equal to the row number of your matrix, can\'t execute the solve function");
        MessagePrinter::exitcfem();
    }
    VectorXd x(getM(),0.0);
    Eigen::MatrixXd A(getM(),getN());
    Eigen::VectorXd B(getM()),X(getN());
    for(int i=0;i<getM();i++){
        for(int j=0;j<getN();j++){
            A.coeffRef(i,j)=(*this)(i,j);
        }
        B.coeffRef(i)=b(i);
    }
    X=A.fullPivLu().solve(B);
    for(int i=0;i<getM();i++){
        x(i)=X.coeff(i);
    }
    return x;
}
VectorXd MatrixXd::operator*(const VectorXd &a)const{
    VectorXd temp(m_m,0.0);
    if(m_n!=a.getM()){
        MessagePrinter::printErrorTxt("A*b should be applied to A matrix with the same cols as b vector!");
        MessagePrinter::exitcfem();
    }
    else{
        for(int i=0;i<m_m;i++){
            temp(i)=0.0;
            for(int j=0;j<m_n;j++){
                temp(i)+=(*this)(i,j)*a(j);
            }
        }
        return temp;
    }
    return temp;
}
VectorXd MatrixXd::operator*(const ViogtRank2Tensor2D &a)const{
    VectorXd temp(m_m,0.0);
    if(m_n!=3){
        MessagePrinter::printErrorTxt("A*b should be applied to A matrix with the same cols as b vector!");
        MessagePrinter::exitcfem();
    }
    else{
        for(int i=0;i<m_m;i++){
            temp(i)=0.0;
            for(int j=0;j<m_n;j++){
                temp(i)+=(*this)(i,j)*a(j);
            }
        }
        return temp;
    }
    return temp;       
}
void MatrixXd::print(){
    for(int rowI=0;rowI<m_m;++rowI){
        for(int colI=0;colI<m_n;++colI){
            printf("%12.5e",(*this)(rowI,colI));
        }
        printf("\n");
    }
}