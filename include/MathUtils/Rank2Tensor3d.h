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
#include "MathUtils/Vector3d.h"
#include "Utils/MessagePrinter.h"

class Rank4Tensor3d;

using std::fill;
using std::sqrt;
using std::abs;

/**
 * This class implement the general manipulation for rank-2 tensor.
 */
class Rank2Tensor3d{
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
    Rank2Tensor3d();
    Rank2Tensor3d(const double val);
    Rank2Tensor3d(const Rank2Tensor3d &a);
    Rank2Tensor3d(const InitMethod &initmethod);
    ~Rank2Tensor3d();
    //**********************************************************************
    //*** for row, col and other elements access
    //**********************************************************************
    /**
     * get the ith row of the rank-2 tensor
     * @param i i-th row number, start from 1 to 3
     */
    inline Vector3d getIthRow(const int &i)const{
        Vector3d temp(0.0);
        temp(0)=(*this)(i,0);
        temp(1)=(*this)(i,1);
        temp(2)=(*this)(i,2);
        return temp;
    }
    /**
     * get the ith column of the rank-2 tensor
     * @param i i-th col number, start from 1 to 3
     */
    inline Vector3d getIthCol(const int &i)const{
        Vector3d temp(0.0);
        temp(0)=(*this)(0,i);
        temp(1)=(*this)(1,i);
        temp(2)=(*this)(2,i);
        return temp;
    }
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
    inline double operator()(const int &i,const int &j) const{
        if(i<0||i>2 || j<0||j>2 ){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" or j="+to_string(j)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitcfem();
        }
        return m_vals[i*N+j];
    }
    /** for index based access(start from 1, instead of zero !!!)
     * @param i i index of the rank-2 tensor, start from 1
     * @param j j index of the rank-2 tensor, start from 1
     */
    inline double& operator()(const int &i,const int &j){
        if(i<0||i>2 || j<0||j>2 ){
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
    inline double operator[](const int &i) const{
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
    inline double& operator[](const int &i){
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
    inline Rank2Tensor3d& operator=(const double &a){
        for(int i=0;i<N2;i++) m_vals[i]=a;
        return *this;
    }
    /**
     * '=' operator for rank-2 tensor
     * @param a the right hand side rank-2 tensor
     */
    inline Rank2Tensor3d& operator=(const Rank2Tensor3d &a){
        for(int i=0;i<N2;i++) m_vals[i]=a.m_vals[i];
        return *this;
    }
    //*******************************
    //*** for +  operator
    //*******************************
    /**
     * '+' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor3d operator+(const double &a) const{
        Rank2Tensor3d temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]+a;
        return temp;
    }
    /**
     * '+' operator for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor3d operator+(const Rank2Tensor3d &a) const{
        Rank2Tensor3d temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]+a.m_vals[i];
        return temp;
    }
    //*******************************
    //*** for +=  operator
    //*******************************
    /**
     * '+=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor3d& operator+=(const double &a) {
        for(int i=0;i<N2;i++) m_vals[i]+=a;
        return *this;
    }
    /**
     * '+=' for rank-2 tensor
     * @param a right hand side tensor
     */
    inline Rank2Tensor3d& operator+=(const Rank2Tensor3d &a){
        for(int i=0;i<N2;i++) m_vals[i]+=a.m_vals[i];
        return *this;
    }
    /**
     * add scalar value to diagnal element, the current tensor=old-tensor+I*a, where I is the identity tensor
     * @param val the scalar value to be added
     */
    inline Rank2Tensor3d& addIa(const double &val){
        (*this)(0,0)+=val;
        (*this)(1,1)+=val;
        (*this)(2,2)+=val;
        return *this;
    }
    //*******************************
    //*** for -  operator
    //*******************************
    /**
     * '-' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor3d operator-(const double &a) const{
        Rank2Tensor3d temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]-a;
        return temp;
    }
    /**
     * '-' operator for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor3d operator-(const Rank2Tensor3d &a) const{
        Rank2Tensor3d temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]-a.m_vals[i];
        return temp;
    }
    //*******************************
    //*** for -=  operator
    //*******************************
    /**
     * '-=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor3d& operator-=(const double &a) {
        for(int i=0;i<N2;i++) m_vals[i]-=a;
        return *this;
    }
    /**
     * '-=' for rank-2 tensor
     * @param a right hand side tensor
     */
    inline Rank2Tensor3d& operator-=(const Rank2Tensor3d &a){
        for(int i=0;i<N2;i++) m_vals[i]-=a.m_vals[i];
        return *this;
    }
    //*******************************
    //*** for *  operator
    //*******************************
    /**
     * '*' for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor3d operator*(const double &a) const{
        Rank2Tensor3d temp(0.0);
        for(int i=0;i<N2;++i) temp.m_vals[i]=m_vals[i]*a;
        return temp;
    }
    /**
     * '*' operator for vector3d
     * @param a right hand side vector3d
     */
    inline Vector3d operator*(const Vector3d &a) const{
        Vector3d temp(0.0);
        for(int i=0;i<N;i++){
            temp(i)=(*this)(i,0)*a(0)+(*this)(i,1)*a(1)+(*this)(i,2)*a(2);
        }
        return temp;
    }
    /**
     * '*' operator for rank-2 tensor, return \f$a_{ij}=b_{ik}c_{kj}\f$
     * @param a the right hand side rank-2 tensor
     */
    inline Rank2Tensor3d operator*(const Rank2Tensor3d &a) const{
        // return A*B(still rank-2 tensor)
        Rank2Tensor3d temp(0.0);
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                temp(i,j)=(*this)(i,0)*a(0,j)+(*this)(i,1)*a(1,j)+(*this)(i,2)*a(2,j);
            }
        }
        return temp;
    }
    /**
     * double dot(:, the double contruction) between two rank-2 tensor, the result is \f$\sum a_{ij}b_{ij}\f$, which is a scalar
     * @param a the right hand side rank-2 tensor
     */
    inline double doubledot(const Rank2Tensor3d &a) const{
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
     * double dot(:, the double contruction) between rank-2 and rank-4 tensor, 
     * the result is \f$\sum a_{ij}c_{ijkl}=b_{kl}\f$, which is a scalar
     * @param a the right hand side rank-2 tensor
     */
    Rank2Tensor3d doubledot(const Rank4Tensor3d &a) const;
    //*******************************
    //*** for *=  operator
    //*******************************
    /**
     * '*=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor3d& operator*=(const double &a) {
        for(int i=0;i<N2;i++) m_vals[i]*=a;
        return *this;
    }
    /**
     * '*=' for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor3d& operator*=(const Rank2Tensor3d &a){
        Rank2Tensor3d temp(0.0);
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
    inline Rank2Tensor3d operator/(const double &a) const{
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for / operator in rank-2 tensor");
            MessagePrinter::exitcfem();
        }
        Rank2Tensor3d temp(0.0);
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
    inline Rank2Tensor3d& operator/=(const double &a){
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
    friend Rank2Tensor3d operator*(const double &lhs,const Rank2Tensor3d &a);
    //*** for left hand vector times rank-2 tensor
    /**
     * '*' operator for left hand side vector3d
     * @param lhs the left hand side vector3d
     * @param a the right hand side rank-2 tensor
     */
    friend Vector3d operator*(const Vector3d &lhs,const Rank2Tensor3d &a);
    //*******************************************************************
    //*** for advanced math operators
    //*******************************************************************
    /**
     * get the exponetial formula of current rank-2 tensor
    */
    Rank2Tensor3d exp()const{
        Rank2Tensor3d I;
        I.setToIdentity();
        return I
              +(*this)
              +(*this)*(*this)*(1.0/(1.0*2.0))
              +(*this)*(*this)*(*this)*(1.0/(1.0*2.0*3.0))
              +(*this)*(*this)*(*this)*(*this)*(1.0/(1.0*2.0*3.0*4.0))
              +(*this)*(*this)*(*this)*(*this)*(*this)*(1.0/(1.0*2.0*3.0*4.0*5.0))
              +(*this)*(*this)*(*this)*(*this)*(*this)*(*this)*(1.0/(1.0*2.0*3.0*4.0*5.0*6.0));
    }
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
    inline void setFromVectorDyad(const Vector3d &a,const Vector3d &b){
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                (*this)(i,j)=a(i)*b(j);
            }
        }
    }
    //**********************************************************************
    //*** for some common mathematic manipulations
    //**********************************************************************
    /**
     * return the trace of a rank-2 tensor, result is \f$\sum a_{ii}\f$
     */
    inline double trace() const{
        return (*this)(0,0)+(*this)(1,1)+(*this)(2,2);
    }
    /**
     * return the determinant of the current rank-2 tensor
     */
    inline double det() const{
        // taken from http://mathworld.wolfram.com/Determinant.html
        // Eq.10
        return (*this)(0,0)*(*this)(1,1)*(*this)(2,2)
              -(*this)(0,0)*(*this)(1,2)*(*this)(2,1)
              -(*this)(0,1)*(*this)(1,0)*(*this)(2,2)
              +(*this)(0,1)*(*this)(1,2)*(*this)(2,0)
              +(*this)(0,2)*(*this)(1,0)*(*this)(2,1)
              -(*this)(0,2)*(*this)(1,1)*(*this)(2,0);
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
     * return the \f$L_{2}\f$ norm^2 of current rank-2 tensor, result is \f$\sqrt{\sum a_{ij}^{2}}\f$
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
        double trAA=((*this)*(*this)).trace();
        double trA=trace();
        return 0.5*(trA*trA-trAA);
    }
    /**
     * return the third invariant of current tensor, namely, \f$I_{3}\f$.
     */
    inline double thirdInvariant() const{
        return det();
    }
    //*** for inverse
    /**
     * return the inverse tensor of current one, namely, \f$\mathbf{A}^{-1}\f$.
     */
    inline Rank2Tensor3d inverse() const{
        double J=det();
        if(abs(J)<1.0e-16){
            MessagePrinter::printErrorTxt("inverse operation failed for a singular rank-2 tensor !");
            MessagePrinter::exitcfem();
        }
        Rank2Tensor3d inv(0.0);
        // taken from wiki:
        //   https://en.wikipedia.org/wiki/Invertible_matrix
        double A= (*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1);
        double D=-(*this)(0,1)*(*this)(2,2)+(*this)(0,2)*(*this)(2,1);
        double G= (*this)(0,1)*(*this)(1,2)-(*this)(0,2)*(*this)(1,1);
        inv(0,0)=A/J;inv(0,1)=D/J;inv(0,2)=G/J;

        double B=-(*this)(1,0)*(*this)(2,2)+(*this)(1,2)*(*this)(2,0);
        double E= (*this)(0,0)*(*this)(2,2)-(*this)(0,2)*(*this)(2,0);
        double H=-(*this)(0,0)*(*this)(1,2)+(*this)(0,2)*(*this)(1,0);
        inv(1,0)=B/J;inv(1,1)=E/J;inv(1,2)=H/J;

        double C= (*this)(1,0)*(*this)(2,1)-(*this)(1,1)*(*this)(2,0);
        double F=-(*this)(0,0)*(*this)(2,1)+(*this)(0,1)*(*this)(2,0);
        double I= (*this)(0,0)*(*this)(1,1)-(*this)(0,1)*(*this)(1,0);
        inv(2,0)=C/J;inv(2,1)=F/J;inv(2,2)=I/J;
        return inv;
    }
    /**
     * return the transpose tensor of current one, namely, \f$\mathbf{A}^{T}\f$.
     */
    inline Rank2Tensor3d transpose() const{
        Rank2Tensor3d temp(0.0);
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
        Rank2Tensor3d temp(0.0);
        temp=(*this).transpose();
        (*this)=temp;
    }
    /**
     * return the deviatoric part of current rank-2 tensor
     */
    inline Rank2Tensor3d dev()const{
        Rank2Tensor3d temp(0.0),I(0.0);
        I.setToIdentity();
        temp=(*this)-I*(this->trace()/3.0);
        return temp;
    }
    //**********************************************************************
    //*** for decomposition
    //**********************************************************************
    /**
     * calculate the eigen value and eigen vector for current rank-2 tensor
     * @param eigval the double array, which stores the eigen value
     * @param eigvec the rank-2 tensor, where each column store the related eigen vector
     */
    void calcEigenValueAndEigenVectors(double (&eigval)[3],Rank2Tensor3d &eigvec) const;
    //**********************************************************************
    //*** for decomposition
    //**********************************************************************
    /**
     * print all the elements to the terminal
     */
    inline void print() const{
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(0,0),(*this)(0,1),(*this)(0,2));
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(1,0),(*this)(1,1),(*this)(1,2));
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(2,0),(*this)(2,1),(*this)(2,2));
    }



private:
    const int N=3;/**< the dimension of current rank-2 tensor */
    const int N2=9;/**< the total length of current rank-2 tensor */
    vector<double> m_vals;/**< vector for the elements of rank-2 tensor */
};