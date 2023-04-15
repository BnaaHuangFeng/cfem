#pragma once
#include "MathUtils/Vector2d.h"
class Vector2d;
/**
 * This class implement the calculation and information storage about the 2D-shapfun (it's a abstract class) 
 */
class Shpfun2D{
public:
    /**
     * Construction function
    */
    Shpfun2D():m_r(0.0){};
    /**
     * @param t_r the natural coords of the point
    */
    inline Shpfun2D(Vector2d t_r):m_r(t_r){};
    /**
     * destructor 
    */
    virtual ~Shpfun2D(){};
    /**
     * init the 2D shpfun, includeing set its natural coords and coordinates of a element's nodes in the reference configuration.
     * @param t_r > the natural coords of the point
     * @param t_x0 > the element's nodes' coords in the ref config (each item is for a node's coord0)
    */
    virtual void init(Vector2d t_r,Vector2d m_r[])=0;
    /**
     * set the element's nodes' coords in the ref config
     * @param t_x0 the element's nodes' coords in the ref config (each item is for a node's coord0)
    */
    virtual void setRefCoords(Vector2d t_x0[])=0;
    /**
     * set the natural coords of the point which to cal the shape function and it's derivates
     * @param t_r the natural coords of the point
    */
    inline void setNatCoords(Vector2d t_r){m_r=t_r;};
    /**
     * get the shpfun value
     * @param t_N > pointer to receive the shp fun val 
    */
    virtual void getShpfunVal(double t_N[])=0;
    /**
     * get coordinates of a element's nodes in the reference configuration, [node id sf 0](dof id sf 1)
     * @param t_x0 > coordinates of a element's nodes in the reference configuration, [node id sf 0](dof id sf 1)
    */
    virtual void getCoords0(Vector2d t_x0[])=0;
    /**
     * get the derivates of shape function to natural coords, [node id sf 0](dof id sf 1)
     * @param t_dNdr > derivates of shape function to natural coords
    */
    virtual void getDer2Nat(Vector2d t_dNdr[])=0;
    /**
     * get the derivates of shape function to reference coords, [node id sf 0](dof id sf 1)
     * @param t_dNdx0 > the derivates of shape function to reference coords, [node id sf 0](dof id sf 1)
    */
    virtual void getDer2Ref(Vector2d t_dNdx0)=0;
public:
    static const int m_dim=2;/**< dimension of shape function, i.e., 1d, 2d, and 3d. */
    Vector2d m_r;/**< the natural coords of the point which to cal the shape function and it's derivates*/
};