#pragma once
class Vector
{
private:
    /* data */
public:
    Vector(){};
    virtual ~Vector(){};
    /**
     * get vector's i-th entry
    */
    virtual double& operator()(const int i) = 0;
    virtual double operator()(const int i)const = 0;
};