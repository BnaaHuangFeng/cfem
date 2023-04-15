#include "MathUtils/ViogtRank2Tensor2D.h"
#include "MathUtils/Rank2Tensor2d.h"
#include <cmath>
ViogtRank2Tensor2D::ViogtRank2Tensor2D(InitMethod initmethod):Vector3d(0.0){
    switch (initmethod)
    {
    case InitMethod::ZERO:
        break;
    case InitMethod::IDENTITY:
        (*this)(1)=1.0;
        (*this)(2)=1.0;
        (*this)(3)=0.0;
        break;
    case InitMethod::RANDOM:
        for(int iViogt=1;iViogt<=NViogt;iViogt++)
            (*this)(iViogt)=static_cast<double>(1.0*rand()/RAND_MAX);
        break;
    default:
        break;
    }
}

double& ViogtRank2Tensor2D::operator()(const int i,const int j){
    if(i<1||i>2||j<1||j>2){
    MessagePrinter::printErrorTxt(to_string(i)+" is out of range for Vector2");
    MessagePrinter::exitAsFem();
}
    if(i==1&&j==1)return Vector3d::operator()(1);
    else if (i!=j)return Vector3d::operator()(3);
    else return Vector3d::operator()(2);
};

double ViogtRank2Tensor2D::operator()(const int i,const int j)const{
    if(i<1||i>2||j<1||j>2){
    MessagePrinter::printErrorTxt(to_string(i)+" is out of range for Vector2");
    MessagePrinter::exitAsFem();
}
    if(i==1&&j==1)return Vector3d::operator()(1);
    else if (i!=j)return Vector3d::operator()(3);
    else return Vector3d::operator()(2);
};

Vector2d ViogtRank2Tensor2D::operator*(const Vector2d &a)const{
    Vector2d tmp;
    tmp(1)=Vector3d::operator()(1)*a(1)+Vector3d::operator()(3)*a(2);
    tmp(2)=Vector3d::operator()(3)*a(1)+Vector3d::operator()(2)*a(2);
    return tmp;
}
Rank2Tensor2d ViogtRank2Tensor2D::operator*(const Rank2Tensor2d &R){
    Rank2Tensor2d tmp(0.0);
    for(int i=1;i<=this->dim;i++){
        for(int j=1;j<=this->dim;j++){
            for(int k=1;k<=this->dim;k++){
                tmp(i,j)+=(*this)(i,k)*R(k,j);
            }
        }
    }
    return tmp;
}
Rank2Tensor2d ViogtRank2Tensor2D::toRank2Tensor2d()const{
    Rank2Tensor2d tmp(0.0);
    for(int i=1;i<=dim;i++){
        for(int j=1;j<=dim;j++){
            tmp(i,j)=(*this)(i,j);
        }
    }
    return tmp;
}
Rank2Tensor ViogtRank2Tensor2D::toRank2Tensor3d()const{
    Rank2Tensor tmp(0.0);
    for(int i=1;i<=dim;i++){
        for(int j=1;j<=dim;j++){
            tmp(i,j)=(*this)(i,j);
        }
    }
    return tmp;
}