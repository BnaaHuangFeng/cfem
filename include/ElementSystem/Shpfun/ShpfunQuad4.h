#pragma once
#include "ElementSystem/Shpfun/Shpfun2D.h"
#include <vector>
#include "InputSystem/EnumDataType.h"
#include "MathUtils/Vector2d.h"
#include "MathUtils/Rank2Tensor2d.h"
/**
 * This class implement the calculation and data storage about the shapfun of  the planar 4-node element with reduced integration
 */
/**
 * COMPUTES SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES FOR
 * ELEMENT 'QUAD_4':
 *                         4         3
 *                          o-------o
 *                          |       |     STANDARD ISOPARAMETRIC
 *                          |       |     BI-LINEAR 4-NODE QUADRILATERAL 
 *                          |       |
 *                          o-------o
 *                         1         2
 *
 * REFERENCE: Expression (4.42)
*/
class ShpfunQuad4:Shpfun2D{
public:
    /**
     * constructor
     */
    ShpfunQuad4();
    /**
     * construction
     * @param t_r > coordinates of a element's nodes in the reference configuration, [node id sf 0](dof id sf 1)
    */
    ShpfunQuad4(Vector2d t_r):Shpfun2D(t_r){};
    /**
     * construction
     * @param t_r > the natural coords of the point
     * @param t_x0 > coordinates of a element's nodes in the reference configuration, [node id sf 0](dof id sf 1)
    */
    ShpfunQuad4(Vector2d t_r, Vector2d t_x0[]);
    /**
     * init the 2D shpfun, includeing set its natural coords and coordinates of a element's nodes in the reference configuration.
     * @param t_r > the natural coords of the point
     * @param t_x0 > the element's nodes' coords in the ref config (each item is for a node's coord0)
    */
    virtual void init(Vector2d t_r,Vector2d t_x0[]);
    /**
     * set the element's nodes' coords in the ref config
     * @param t_x0 the element's nodes' coords in the ref config (each item is for a node's coord0)
    */
    virtual void setRefCoords(Vector2d t_x0[]);
    /**
     * get the shpfun value
     * @param t_N > pointer to receive the shp fun val 
    */
    virtual void getShpfunVal(double t_N[]);
    /**
     * get coordinates of a element's nodes in the reference configuration, [node id sf 0](dof id sf 1)
     * @param t_x0 > coordinates of a element's nodes in the reference configuration, [node id sf 0](dof id sf 1)
    */
    virtual void getCoords0(Vector2d t_x0[]);
    /**
     * get the derivates of shape function to natural coords, [node id sf 0](dof id sf 1)
     * @param t_dNdr > derivates of shape function to natural coords
    */
    virtual void getDer2Nat(Vector2d t_dNdr[]);
    /**
     * get the derivates of shape function to reference coords, [node id sf 0](dof id sf 1)
     * @param t_dNdx0 > the derivates of shape function to reference coords, [node id sf 0](dof id sf 1)
    */
    virtual void getDer2Ref(Vector2d t_dNdx0[]);
public:
    static const MeshType m_mesh_type=MeshType::QUAD4;/**< the type of mesh */
    static const int m_funs=4;/**< number of shape functions */
    Vector2d m_x0[m_funs];/**< coordinates of a element's nodes in the reference configuration, [node id sf 0](dof id sf 1)*/
    Vector2d m_dNdr[m_funs];/**< the derivates of shape function to natural coords, [node id sf 0](dof id sf 1)*/
    Vector2d m_dNdx0[m_funs];/**< the derivates of shape function to reference coords, [node id sf 0](dof id sf 1)*/
};