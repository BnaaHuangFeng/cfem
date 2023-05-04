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

#pragma once

#include <iostream>
#include <cmath>
#include "MathUtils/Vector3d.h"
#include "MathUtils/Rank2Tensor3d.h"
#include "Utils/MessagePrinter.h"

class Rank2Tensor3d;

using std::fill;
using std::sqrt;
using std::abs;

/**
 * This class implement the general manipulation for rank-2 tensor.
 */
class Rank4Tensor3d{
public:
    /**
     * enum for different initializing methods
     */
    enum InitMethod{
        ZERO,
        IDENTITY,
        IDENTITY4,
        IDENTITY4TRANS,
        IDENTITY4SYMMETRIC,
        RANDOM
    };
public:
    /**
     * constructor
     */
    Rank4Tensor3d();
    Rank4Tensor3d(const Rank4Tensor3d &a);
    Rank4Tensor3d(const InitMethod &method);
    ~Rank4Tensor3d();

    //**********************************************************************
    //*** for operator override
    //**********************************************************************
    //*******************************
    //*** for ()  operator
    //*******************************
    /**
     * () operator for rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index, start from 1 instead of 0 !!!
     * @param j j index, start from 1 instead of 0 !!!
     * @param k k index, start from 1 instead of 0 !!!
     * @param l l index, start from 1 instead of 0 !!!
     */
    inline double operator()(const int &i,const int &j,const int &k,const int &l) const{
        if(i<0||i>2 || j<0||j>2 || k<0||k>2 || l<0||l>2){
            MessagePrinter::printErrorTxt("your i or j or k or l is out of range when you call a rank-4 tensor");
            MessagePrinter::exitcfem();
        }
        return m_vals[((i*N+j)*N+k)*N+l];
    }
    /**
     * () operator for rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index, start from 1 instead of 0 !!!
     * @param j j index, start from 1 instead of 0 !!!
     * @param k k index, start from 1 instead of 0 !!!
     * @param l l index, start from 1 instead of 0 !!!
     */
    inline double& operator()(const int &i,const int &j,const int &k,const int &l){
        if(i<0||i>2 || j<0||j>2 || k<0||k>2 || l<0||l>2){
            MessagePrinter::printErrorTxt("your i or j or k or l is out of range when you call a rank-4 tensor");
            MessagePrinter::exitcfem();
        }
        return m_vals[((i*N+j)*N+k)*N+l];
    }
    /**
     * get the voigt component of current rank-4 tensor
     * @param i i-index for 1st dimension
     * @param j j-index for 2nd dimension
     */
    double getVoigtComponent(const int &i,const int &j)const;

    /**
     * return the reference via the voigt component of current rank-4 tensor
     * Hf: K_ik=N_pi,q*D_pqrs*N_rk,s
     * @param i i-index for 1st dimension
     * @param j j-index for 2nd dimension
     */
    double& voigtComponent(const int &i,const int &j);
    //*******************************
    //*** for []  operator
    //*******************************
    /**
     * [] operator to get the component of a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index
     */
    inline double  operator[](const int &i) const{
        return m_vals[i];
    }
    /**
     * [] operator to get the component of a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index
     */
    inline double& operator[](const int &i){
        return m_vals[i];
    }
    //*******************************
    //*** for =  operator
    //*******************************
    /**
     * = operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d& operator=(const double &a){
        fill(m_vals.begin(),m_vals.end(),a);
        return *this;
    }
    /**
     * = operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand rank-4 tensor
     */
    inline Rank4Tensor3d& operator=(const Rank4Tensor3d &a){
        for(int i=0;i<N4;i++) m_vals[i]=a.m_vals[i];
        return *this;
    }
    //*******************************
    //*** for + operator
    //*******************************
    /**
     * + operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d operator+(const double &a) const{
        Rank4Tensor3d temp(Rank4Tensor3d::InitMethod::ZERO);
        for(int i=0;i<N4;i++) temp.m_vals[i]=m_vals[i]+a;
        return temp;
    }
    /**
     * + operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand rank-4 tensor
     */
    inline Rank4Tensor3d operator+(const Rank4Tensor3d &a) const{
        Rank4Tensor3d temp(Rank4Tensor3d::InitMethod::ZERO);
        for(int i=0;i<N4;i++) temp.m_vals[i]=m_vals[i]+a.m_vals[i];
        return temp;
    }
    //*******************************
    //*** for += operator
    //*******************************
    /**
     * += operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d& operator+=(const double &a){
        for(int i=0;i<N4;i++) m_vals[i]=m_vals[i]+a;
        return *this;
    }
    /**
     * += operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side rank-4 tensor
     */
    inline Rank4Tensor3d& operator+=(const Rank4Tensor3d &a){
        for(int i=0;i<N4;i++) m_vals[i]=m_vals[i]+a.m_vals[i];
        return *this;
    }
    //*******************************
    //*** for - operator
    //*******************************
    /**
     * - operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d operator-(const double &a) const{
        Rank4Tensor3d temp(Rank4Tensor3d::InitMethod::ZERO);
        for(int i=0;i<N4;i++) temp.m_vals[i]=m_vals[i]-a;
        return temp;
    }
    /**
     * - operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand rank-4 tensor
     */
    inline Rank4Tensor3d operator-(const Rank4Tensor3d &a) const{
        Rank4Tensor3d temp(Rank4Tensor3d::InitMethod::ZERO);
        for(int i=0;i<N4;i++) temp.m_vals[i]=m_vals[i]-a.m_vals[i];
        return temp;
    }
    //*******************************
    //*** for -= operator
    //*******************************
    /**
     * -= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d& operator-=(const double &a){
        for(int i=0;i<N4;i++) m_vals[i]=m_vals[i]-a;
        return *this;
    }
    /**
     * -= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side rank-4 tensor
     */
    inline Rank4Tensor3d& operator-=(const Rank4Tensor3d &a){
        for(int i=0;i<N4;i++) m_vals[i]=m_vals[i]-a.m_vals[i];
        return *this;
    }
    //*******************************
    //*** for * operator
    //*******************************
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d operator*(const double &a) const{
        Rank4Tensor3d temp(Rank4Tensor3d::InitMethod::ZERO);
        for(int i=0;i<N4;i++) temp.m_vals[i]=m_vals[i]*a;
        return temp;
    }
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}=\mathbb{C}_{ijkm}\mathbf{A}_{ml}\f$
     * @param a right hand side ran-2 tensor
     */
    Rank4Tensor3d operator*(const Rank2Tensor3d &a) const;
    /**
     * double dot : between a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$ and rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    Rank2Tensor3d doubledot(const Rank2Tensor3d &a) const;
    /**
     * double dot : between a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$ and another rank-4 tensor
     * return \f$\mathbb{C}_{ijkl}=\mathbb{A}_{ijmn}\mathbb{B}_{mnkl}\f$
     * @param a right hand side rank-4 tensor
     */
    Rank4Tensor3d doubledot(const Rank4Tensor3d &a) const;
    //*******************************
    //*** for *= operator
    //*******************************
    /**
     * *= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d& operator*=(const double &a){
        for(int i=0;i<N4;++i) m_vals[i]=m_vals[i]*a;
        return *this;
    }
    //*******************************
    //*** for / operator
    //*******************************
    /**
     * / operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d operator/(const double &a) const{
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for / operator in rank-4 tensor");
            MessagePrinter::exitcfem();
        }
        Rank4Tensor3d temp(Rank4Tensor3d::InitMethod::ZERO);
        for(int i=0;i<N4;i++) temp.m_vals[i]=m_vals[i]/a;
        return temp;
    }
    //*******************************
    //*** for /= operator
    //*******************************
    /**
     * /= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor3d& operator/=(const double &a){
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for /= operator in rank-4 tensor");
            MessagePrinter::exitcfem();
        }
        for(int i=0;i<N4;++i) m_vals[i]=m_vals[i]*a;
        return *this;
    }
    //*******************************************************************
    //*** for left hand side manipulation
    //*******************************************************************
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param lhs left hand side scalar
     * @param a right hand side rank-4 tensor
     */
    friend Rank4Tensor3d operator*(const double &lhs,const Rank4Tensor3d &a);
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param lhs left hand side rank-2 tensor
     * @param a right hand side rank-4 tensor
     */
    friend Rank4Tensor3d operator*(const Rank2Tensor3d &lhs,const Rank4Tensor3d &a);
    //**********************************************************************
    //*** for general settings
    //**********************************************************************
    /**
     * set current rank-4 tensor to 0
     */
    inline void setToZeros(){
        fill(m_vals.begin(),m_vals.end(),0.0);
    }
    /**
     * set current rank-4 tensor to idendity
     */
    inline void setToIdentity(){
        setToZeros();
        for(int i=0;i<N;i++) (*this)(i,i,i,i)=1.0;
    }
    /**
     * set current rank-4 tensor to rank-4 identity
     */
    inline void setToIdentity4(){
        // maps a rank2 tensor to itself(no symmetric consideriation here),i.e. Iden4:rank2=rank2
        setToZeros();
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                for(int k=0;k<N;k++){
                    for(int l=0;l<N;l++){
                        (*this)(i,j,k,l)=1.0*((i==k)&&(j==l));
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to rank-4 transposed identity
     */
    inline void setIdentity4Transpose(){
        // maps a rank-2 tensor to its transpose, A^{T}=I4:A
        setToZeros();
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                for(int k=0;k<N;k++){
                    for(int l=0;l<N;l++){
                        (*this)(i,j,k,l)=1.0*((j==k)&&(i==l));
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to symmetric identity rank-4 tensor
     */
    inline void setToIdentity4Symmetric(){
        // symmetric fourth-order tensor
        setToZeros();
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                for(int k=0;k<N;k++){
                    for(int l=0;l<N;l++){
                        (*this)(i,j,k,l)=0.5*((i==k)&&(j==l))
                                        +0.5*((i==l)&&(j==k));
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to random values
     */
    inline void setToRandom(){
        srand(time(0));
        for(int i=0;i<N;++i){
            for(int j=0;j<N;++j){
                for(int k=0;k<N;++k){
                    for(int l=0;l<N;++l){
                        (*this)(i,j,k,l)=static_cast<double>(1.0*rand()/RAND_MAX);
                    }
                }
            }
        }
    }

private:
    constexpr static const int viogtInd2ij[6][2]={{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}};
    const int N=3;/**< the dimension of current tensor */
    const int N4=81;/**< the total length of current tensor */
    vector<double> m_vals;/**< the tensor components */

};