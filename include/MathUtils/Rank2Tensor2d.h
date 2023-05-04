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
#pragma once
#include <iostream>
#include <cmath>
#include "MathUtils/ViogtRank2Tensor2D.h"
#include "Utils/MessagePrinter.h"
#include "Eigen/Eigen"
using std::fill;
using std::sqrt;
using std::abs;
class Vector2d;
class ViogtRank2Tensor2D;
class Rank2Tensor3d;
/**
 * This class implement the general manipulation for rank-2 tensor.
 */
class Rank2Tensor2d{
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
    Rank2Tensor2d();
    Rank2Tensor2d(const Rank2Tensor2d &a);
    Rank2Tensor2d(const ViogtRank2Tensor2D &a);
    Rank2Tensor2d(const InitMethod &initmethod);
    ~Rank2Tensor2d();
    Rank2Tensor3d toRank2Tensor3d()const;
    //**********************************************************************
    //*** for row, col and other elements access
    //**********************************************************************
    /**
     * get the ith row of the rank-2 tensor
     * @param i i-th row number, start from 1 to 3
     */
    Vector2d getIthRow(const int i)const;
    /**
     * get the ith column of the rank-2 tensor
     * @param i i-th col number, start from 1 to 3
     */
    Vector2d getIthCol(const int i)const;
    //**********************************************************************
    //*** for operator override
    //**********************************************************************
    //*******************************
    //*** for ()  operator
    //*******************************
    /** for index based access(start from 1, instead of zero !!!)
     * @param i i index of the rank-2 tensor, start from 1
     * @param j j index of the rank-2 tensor, start from 1
     */
    inline double operator()(const int i,const int j) const{
        if(i<0||i>=N || j<0||j>=N ){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" or j="+to_string(j)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitcfem();
        }
        return m_vals[i*N+j];
    }
    /** for index based access(start from 1, instead of zero !!!)
     * @param i i index of the rank-2 tensor, start from 1
     * @param j j index of the rank-2 tensor, start from 1
     */
    inline double& operator()(const int i,const int j){
        if(i<0||i>=N || j<0||j>=N ){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" or j="+to_string(j)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitcfem();
        }
        return m_vals[i*N+j];
    }
    //*******************************
    //*** for []  operator
    //*******************************
    /**
     * [] operator for element access
     * @param i global index, range from 1~9
     */
    inline double operator[](const int i) const{
        if(i<0||i>=N2){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitcfem();
        }
        return m_vals[i];
    }
    /**
     * [] operator for element access
     * @param i global index, range from 1~9
     */
    inline double& operator[](const int i){
        if(i<0||i>=N2){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitcfem();
        }
        return m_vals[i];
    }
    //*******************************
    //*** for =  operator
    //*******************************
    /**
     * '=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d& operator=(const double a){
        for(int i=0;i<N2;i++) m_vals[i]=a;
        return *this;
    }
    /**
     * '=' operator for rank-2 tensor
     * @param a the right hand side rank-2 tensor
     */
    inline Rank2Tensor2d& operator=(const Rank2Tensor2d &a){
        for(int i=0;i<N2;i++) m_vals[i]=a[i];
        return *this;
    }
    //*******************************
    //*** for +  operator
    //*******************************
    /**
     * '+' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d operator+(const double &a) const{
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]+a;
        return temp;
    }
    /**
     * '+' operator for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor2d operator+(const Rank2Tensor2d &a) const{
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]+a[i];
        return temp;
    }
    //*******************************
    //*** for +=  operator
    //*******************************
    /**
     * '+=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d& operator+=(const double a) {
        for(int i=0;i<N2;i++) m_vals[i]+=a;
        return *this;
    }
    /**
     * '+=' for rank-2 tensor
     * @param a right hand side tensor
     */
    inline Rank2Tensor2d& operator+=(const Rank2Tensor2d &a){
        for(int i=0;i<N2;i++) m_vals[i]+=a[i];
        return *this;
    }
    /**
     * add scalar value to diagnal element, the current tensor=old-tensor+I*a, where I is the identity tensor
     * @param val the scalar value to be added
     */
    inline Rank2Tensor2d& addIa(const double val){
        (*this)(0,0)+=val;
        (*this)(1,1)+=val;
        return *this;
    }
    //*******************************
    //*** for -  operator
    //*******************************
    /**
     * '-' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d operator-(const double a) const{
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]-a;
        return temp;
    }
    /**
     * '-' operator for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor2d operator-(const Rank2Tensor2d &a) const{
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]-a[i];
        return temp;
    }
    //*******************************
    //*** for -=  operator
    //*******************************
    /**
     * '-=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d& operator-=(const double a) {
        for(int i=0;i<N2;i++) m_vals[i]-=a;
        return *this;
    }
    /**
     * '-=' for rank-2 tensor
     * @param a right hand side tensor
     */
    inline Rank2Tensor2d& operator-=(const Rank2Tensor2d &a){
        for(int i=0;i<N2;i++) m_vals[i]-=a[i];
        return *this;
    }
    //*******************************
    //*** for *  operator
    //*******************************
    /**
     * '*' for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d operator*(const double a) const{
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        for(int i=0;i<N2;++i) temp.m_vals[i]=m_vals[i]*a;
        return temp;
    }
    /**
     * '*' operator for vector2d
     * @param a right hand side vector2d
     */
    Vector2d operator*(const Vector2d &a) const;
    /**
     * '*' operator for rank-2 tensor, return \f$a_{ij}=b_{ik}c_{kj}\f$
     * @param a the right hand side rank-2 tensor
     */
    inline Rank2Tensor2d operator*(const Rank2Tensor2d &a) const{
        // return A*B(still rank-2 tensor)
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                temp(i,j)=(*this)(i,0)*a(0,j)+(*this)(i,1)*a(1,j);
            }
        }
        return temp;
    }
    /**
     * '*' operator for vioget rank-2 tensor, return \f$a_{ij}=b_{ik}c_{kj}\f$
     * @param a the right hand side rank-2 tensor
     */
    Rank2Tensor2d operator*(const ViogtRank2Tensor2D &a) const;
    /**
     * double dot(:, the double contruction) between two rank-2 tensor, the result is \f$\sum a_{ij}b_{ij}\f$, which is a scalar
     * @param a the right hand side rank-2 tensor
     */
    inline double doubledot(const Rank2Tensor2d &a) const{
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
    /**
     * double dot(:, the double contruction) between two rank-2 tensor, the result is \f$\sum a_{ij}b_{ij}\f$, which is a scalar
     * @param a the right hand side rank-2 tensor
     */
    double doubledot(const ViogtRank2Tensor2D &a) const;
    //*******************************
    //*** for *=  operator
    //*******************************
    /**
     * '*=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d& operator*=(const double &a) {
        for(int i=0;i<N2;i++) m_vals[i]*=a;
        return *this;
    }
    /**
     * '*=' for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor2d& operator*=(const Rank2Tensor2d &a){
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        temp=(*this)*a;
        (*this)=temp;
        return *this;
    }
    //*******************************
    //*** for /  operator
    //*******************************
    /**
     * '/' for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d operator/(const double a) const{
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for / operator in rank-2 tensor");
            MessagePrinter::exitcfem();
        }
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]/a;
        return temp;
    }
    //*******************************
    //*** for /=  operator
    //*******************************
    /**
     * '/=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor2d& operator/=(const double a){
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for /= operator in rank-2 tensor");
            MessagePrinter::exitcfem();
        }
        for(int i=0;i<N2;i++) m_vals[i]/=a;
        return *this;
    }
    //*******************************************************************
    //*** for left hand side manipulation
    //*******************************************************************
    /**
     * '*' for left hand side scalar
     * @param lhs left hand side scalar value
     * @param a right hand side rank-2 tensor
     */
    friend Rank2Tensor2d operator*(const double lhs,const Rank2Tensor2d &a);
    //*** for left hand vector times rank-2 tensor
    /**
     * '*' operator for left hand side vector3d
     * @param lhs the left hand side vector3d
     * @param a the right hand side rank-2 tensor
     */
    friend Vector2d operator*(const Vector2d &lhs,const Rank2Tensor2d &a);
    //**********************************************************************
    //*** for general settings
    //**********************************************************************
    /**
     * set the elements of current rank-2 tensor to be zero
     */
    inline void setToZeros(){
        for(int i=0;i<N2;i++) m_vals[i]=0.0;
    }
    /**
     * set current rannk-2 tensor to be an identitiy tensor, where \f$a_{ij}=\delta_{ij}\f$.
     */
    inline void setToIdentity(){
        for(int i=0;i<N;++i){
            for(int j=0;j<N;++j){
                if(i==j){
                    (*this)(i,j)=1.0;
                }
                else{
                    (*this)(i,j)=0.0;
                }
            }
        }
    }
    /**
     * set current rank-2 tensor to be a random one
     */
    inline void setToRandom(){
        srand(time(0));
        for(int i=0;i<N;++i){
            for(int j=0;j<N;++j){
                (*this)(i,j)=static_cast<double>(1.0*rand()/RAND_MAX);
            }
        }
    }
    /**
     * cross dot \f$\otimes\f$ for two vector(from STL), set current one to \f$c_{ij}=a_{i}\times b_{j}\f$.
     * @param a vector<double> for 1st dimension
     * @param b vector<double> for 2nd dimension
     */
    inline void setFromVectorDyad(const vector<double> &a,const vector<double> &b){
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                (*this)(i,j)=a[i]*b[j];
            }
        }
    }
    /**
     * fill up current rank-2 tensor from the dyad of two vector
     * @param a the first vector
     * @param b the second vector
     */ 
    void setFromVectorDyad(const Vector2d &a,const Vector2d &b);
    //**********************************************************************
    //*** for some common mathematic manipulations
    //**********************************************************************
    /**
     * return the trace of a rank-2 tensor, result is \f$\sum a_{ii}\f$
     */
    inline double trace() const{
        return (*this)(0,0)+(*this)(1,1);
    }
    /**
     * return the determinant of the current rank-2 tensor
     */
    inline double det() const{
        // taken from http://mathworld.wolfram.com/Determinant.html
        // Eq.10
        return (*this)(0,0)*(*this)(1,1)-(*this)(0,1)*(*this)(1,0);
    }
    /**
     * return the \f$L_{2}\f$ norm of current rank-2 tensor, result is \f$\sqrt{\sum a_{ij}^{2}}\f$
     */
    inline double norm() const{
        double sum=0.0;
        for(int i=0;i<N2;i++) sum+=static_cast<double>(m_vals[i]*m_vals[i]);
        return sqrt(sum);
    }
    
    /**
     * return the \f$L_{2}\f$ norm^2 of current rank-2 tensor, result is \f$\sum a_{ij}^{2}\f$
     */
    inline double normsq() const{
        double sum=0.0;
        for(int i=0;i<N2;i++) sum+=static_cast<double>(m_vals[i]*m_vals[i]);
        return sum;
    }
    //*** for the different invariants of stress(strain)
    /**
     * return the first invariant of current rank-2 tensor, namely, \f$I_{1}\f$.
     */
    inline double firstInvariant() const{
        return trace();
    }
    /**
     * return the second invariant of current tensor, namely, \f$I_{2}\f$.
     */
    inline double secondInvariant() const{
        return det();
    }
    //*** for inverse
    /**
     * return the inverse tensor of current one, namely, \f$\mathbf{A}^{-1}\f$.
     */
    inline Rank2Tensor2d inverse() const{
        double J=det();
        if(abs(J)<1.0e-16){
            MessagePrinter::printErrorTxt("inverse operation failed for a singular rank-2 tensor !");
            MessagePrinter::exitcfem();
        }
        Rank2Tensor2d inv(Rank2Tensor2d::InitMethod::ZERO);
        inv(0,0)=(*this)(1,1)/J;
        inv(0,1)=-(*this)(1,0)/J;
        inv(1,0)=-(*this)(0,1)/J;
        inv(1,1)=(*this)(0,0)/J;
        return inv;
    }
    /**
     * return the transpose tensor of current one, namely, \f$\mathbf{A}^{T}\f$.
     */
    inline Rank2Tensor2d transpose() const{
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                temp(i,j)=(*this)(j,i);
            }
        }
        return temp;
    }
    /**
     * transpose current tensor, and overwrite the original one
     */
    inline void transposed(){
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO);
        temp=(*this).transpose();
        (*this)=temp;
    }
    /**
     * return the deviatoric part of current rank-2 tensor
     */
    inline Rank2Tensor2d dev()const{
        Rank2Tensor2d temp(Rank2Tensor2d::InitMethod::ZERO),I(Rank2Tensor2d::InitMethod::ZERO);
        I.setToIdentity();
        temp=(*this)-I*(this->trace()/this->N);
        return temp;
    }
    //**********************************************************************
    //*** for decomposition
    //**********************************************************************
    /**
     * calculate the eigen value and eigen vector for current rank-2 tensor
     * @param eigvalPtr the double array, which stores the eigen value
     * @param eigvecPtr the related eigen vector2d
     */
    void calcEigenValueAndEigenVectors(double eigvalPtr[2],Vector2d eigvecPtr[2]) const;
    //**********************************************************************
    //*** for decomposition
    //**********************************************************************
    /**
     * spectral decomposition for symmtry matrix
     * @param eigvalPtr < ref to eigen values
     * @param eigprojPtr < ref to corresponding eigenvalues projection
     * @param repeated < if the eigen values repeated
    */
    void spectralDecomposition(double eigvalPtr[2],ViogtRank2Tensor2D eigprojPtr[2],bool &repeated)const;
    /**
     * print all the elements to the terminal
     */
    inline void print() const{
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ***\n",(*this)(0,0),(*this)(0,1));
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ***\n",(*this)(1,0),(*this)(1,1));
    }
public:
    static const int N=2;                                               /**< the dimension of current rank-2 tensor */
    static const int N2=4;                                              /**< the total length of current rank-2 tensor */
    double m_vals[N2];                                                  /**< vector for the elements of rank-2 tensor */
};