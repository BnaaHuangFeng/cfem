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

#include "MathUtils/Vector3d.h"

Vector3d::Vector3d(){
    m_vals[0]=0.0;m_vals[1]=0.0;m_vals[2]=0.0;
}
Vector3d::Vector3d(const Vector3d &a){
    m_vals[0]=a.m_vals[0];m_vals[1]=a.m_vals[1];m_vals[2]=a.m_vals[2];
}
Vector3d::Vector3d(const InitMethod initmethod){
    switch (initmethod)
    {
    case Vector3d::InitMethod::ZERO:
        m_vals[0]=0.0;  m_vals[1]=0.0;  m_vals[2]=0.0;
        break;
    case Vector3d::InitMethod::IDENTITY:
        m_vals[0]=1.0;  m_vals[1]=1.0;  m_vals[2]=1.0;
        break;
    case Vector3d::InitMethod::RANDOM:
        m_vals[0]=static_cast<double>(1.0*rand()/RAND_MAX);
        m_vals[1]=static_cast<double>(1.0*rand()/RAND_MAX);
        m_vals[2]=static_cast<double>(1.0*rand()/RAND_MAX);
        break;
    default:
        break;
    }    
}
Vector3d operator*(const double &val,const Vector3d &a){
    Vector3d temp(Vector3d::InitMethod::ZERO);
    temp(0)=val*a(0);temp(1)=val*a(1);temp(2)=val*a(2);
    return temp;
}
