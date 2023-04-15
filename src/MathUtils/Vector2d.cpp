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
//+++ Date   : 2022.05.13
//+++ Purpose: defines the vector with only 3-components, this will
//+++          be frequently used in shape function calculation.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/Vector2d.h"

Vector2d::Vector2d(){
    m_vals[0]=0.0;m_vals[1]=0.0;
}
Vector2d::Vector2d(const double val){
    m_vals[0]=val;m_vals[1]=val;
}
Vector2d::Vector2d(const Vector2d &a){
    m_vals[0]=a.m_vals[0];m_vals[1]=a.m_vals[1];
}
Vector2d::Vector2d(const InitMethod initmethod){
    switch (initmethod)
    {
    case Vector2d::InitMethod::ZERO:
        Vector2d(0.0);
        break;
    case Vector2d::InitMethod::IDENTITY:
        Vector2d(1.0);
        break;
    case Vector2d::InitMethod::RANDOM:
        m_vals[0]=static_cast<double>(1.0*rand()/RAND_MAX);
        m_vals[1]=static_cast<double>(1.0*rand()/RAND_MAX);
        break;
    default:
        break;
    }
}
Vector3d Vector2d::toVector3d()const{
    Vector3d tmp;
    tmp(1)=(*this)(1);
    tmp(2)=(*this)(2);
    tmp(3)=0.0;
    return tmp;
}
Vector2d operator*(const double val,const Vector2d &a){
    Vector2d temp(0.0);
    temp(1)=val*a(1);temp(2)=val*a(2);
    return temp;
}