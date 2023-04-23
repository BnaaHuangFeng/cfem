#pragma once
class Rank2Tensor
{
public:
    Rank2Tensor(){}
    virtual ~Rank2Tensor(){}
    double & operator()(const int i,const int j){}
};