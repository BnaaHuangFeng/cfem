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
//+++ Purpose: Implement rank-4 tensor class for the common
//+++          tensor manipulation in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "MathUtils/Rank4Tensor3d.h"
Rank4Tensor3d::Rank4Tensor3d(){
    m_vals.resize(81,0.0);
}
Rank4Tensor3d::Rank4Tensor3d(const double val){
    m_vals.resize(81,val);
}
Rank4Tensor3d::Rank4Tensor3d(const Rank4Tensor3d &a){
    m_vals.resize(81,0.0);
    for(int i=0;i<N4;i++) m_vals[i]=a.m_vals[i];
}
Rank4Tensor3d::Rank4Tensor3d(const InitMethod &method){
    if(method==InitMethod::ZERO){
        m_vals.resize(81,0.0);
    }
    else if(method==InitMethod::IDENTITY){
        m_vals.resize(81,0.0);
        setToIdentity();
    }
    else if(method==InitMethod::IDENTITY4){
        m_vals.resize(81,0.0);
        setToIdentity4();
    }
    else if(method==InitMethod::IDENTITY4SYMMETRIC){
        m_vals.resize(81,0.0);
        setToIdentity4Symmetric();
    }
    else if(method==InitMethod::IDENTITY4TRANS){
        m_vals.resize(81,0.0);
        setIdentity4Transpose();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported initialize method in rank-4 tensor");
        MessagePrinter::exitcfem();
    }
}
Rank4Tensor3d::~Rank4Tensor3d(){
    m_vals.clear();
}
//************************************************************************
double Rank4Tensor3d::getVoigtComponent(const int &i,const int &j)const{
    if(i<0||i>=6||j<0||j>=6){
        MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
        MessagePrinter::exitcfem();
    }
    else{
        return (*this)(viogtInd2ij[i][0],viogtInd2ij[i][1],viogtInd2ij[j][0],viogtInd2ij[j][1]);
    }
    return (*this)(0,0,0,0);
}
//************************************************************************
double& Rank4Tensor3d::voigtComponent(const int &i,const int &j){
   if(i<0||i>=6||j<0||j>=6){
        MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
        MessagePrinter::exitcfem();
    }
    else{
        return (*this)(viogtInd2ij[i][0],viogtInd2ij[i][1],viogtInd2ij[j][0],viogtInd2ij[j][1]);
    }
    return (*this)(0,0,0,0);
}
//************************************************************************
Rank4Tensor3d Rank4Tensor3d::operator*(const Rank2Tensor3d &a) const{
    Rank4Tensor3d temp(0.0);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                for(int l=0;l<N;l++){
                    temp(i,j,k,l)=0.0;
                    for(int m=0;m<N;m++){
                        temp(i,j,k,l)+=(*this)(i,j,k,m)*a(m,l);
                    }
                }
            }
        }
    }
    return temp;
}
Rank2Tensor3d Rank4Tensor3d::doubledot(const Rank2Tensor3d &a) const{
    Rank2Tensor3d temp(0.0);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            temp(i,j)=0.0;
            for(int k=0;k<N;k++){
                for(int l=0;l<N;l++){
                    temp(i,j)+=(*this)(i,j,k,l)*a(k,l);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor3d Rank4Tensor3d::doubledot(const Rank4Tensor3d &a) const{
    Rank4Tensor3d temp(0.0);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                for(int l=0;l<N;l++){
                    temp(i,j,k,l)=0.0;
                    for(int m=0;m<N;m++){
                        for(int n=0;n<N;n++){
                            temp(i,j,k,l)+=(*this)(i,j,m,n)*a(m,n,k,l);
                        }
                    }
                }
            }
        }
    }
    return temp;
}
//*******************************************************************
//*** for left hand side manipulation
//*******************************************************************
Rank4Tensor3d operator*(const double &lhs,const Rank4Tensor3d &a){
    Rank4Tensor3d temp(0.0);
    for(int i=0;i<a.N4;i++) temp.m_vals[i]=lhs*a.m_vals[i];
    return temp;
}
Rank4Tensor3d operator*(const Rank2Tensor3d &lhs,const Rank4Tensor3d &a){
    // C_ijkl=B_ip*A_pjkl
    Rank4Tensor3d temp(0.0);
    for(int i=0;i<a.N;i++){
        for(int j=0;j<a.N;j++){
            for(int k=0;k<a.N;k++){
                for(int l=0;l<a.N;l++){
                    temp(i,j,k,l)=0.0;
                    for(int p=0;p<a.N;p++){
                        temp(i,j,k,l)+=lhs(i,p)*a(p,j,k,l);
                    }
                }
            }
        }
    }
    return temp;
}