# pragma once
# include "petsc.h"
class LoadSystem
{
private:
    /* data */
public:
    LoadSystem(/* args */);
    ~LoadSystem();
public:
    DM *dmPtr;          /**< ptr to dm load Vec relied on*/
    Vec *loadVecPtr;    /**< global load Vec ptr*/
};