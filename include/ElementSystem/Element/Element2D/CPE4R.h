# pragma once
# include "petsc.h"
# include "ElementSystem/Shpfun/ShpfunQuad4.h"
# include "MaterialSystem/Material2D.h"
# include "ElementSystem/Element/element.h"
/**************************************************
 *                         4         3
 *                          o-------o
 *                          |       |     STANDARD ISOPARAMETRIC
 *                          |   X   |     BI-LINEAR 4-NODE QUADRILATERAL 
 *                          |  I    |
 *                          o-------o
 *                         1         2
 *
 **************************************************
 * 4 node plane strain reduced quadrature elemnt***
 *************************************************/
class CPE4R:element{
    public:
    CPE4R():element(false){}
    CPE4R(bool nLarge):element(nLarge){}
    public:/**< inherent virtual func need to be implemented*/
    /**
     * init element
     * @param matPtr > material ptr
     * @param nLarge > large strain flag
     * @param elmtParamPtr > ptr to the elmt's params
    */
    virtual PetscErrorCode initElement(PetscInt t_elmt_rId, bool nLarge,PetscScalar *elmtParamPtr);
    /**
     * get the elmt's inner force
     * @param t_elmtCoord2 > ptr to the elmt's last converged coords (Vector2d *, Vector3d *)
     * @param t_elmtDofInc > ptr to the elmt's incremental dof values of this inrement until now (Vector2d *, Vector3d *)
     * @param t_elmtInnerForce < ptr to receive the elmt's inner force (need preallocation)
     * @param t_converged < if material update converged
    */
    virtual PetscErrorCode getElmtInnerForce(void *t_elmtCoord2, void *t_elmtDofInc, VectorXd *t_elmtInnerForce, bool *t_converged);
    /**
     * get the elmt's stiffness matrix
     * @param t_elmtCoord2 > ptr to the elmt's last converged coords (Vector2d *, Vector3d *)
     * @param t_elmtDofInc > ptr to the elmt's incremental dof values of this inrement until now (Vector2d *, Vector3d *)
     * @param t_stfMatrix < ptr to receive the elmt's stiffness matrix (need preallocation)
    */
    virtual PetscErrorCode getElmtStfMatrix(void *t_elmtCoord2, void *t_elmtDofInc, MatrixXd *t_stfMatrix);
    public:/**< static member (all elements of this kind share them)*/
    static const int m_dim;                     /**< element dimension*/
    static const int m_mDof_node;               /**< dof num per node*/
    static const int m_mNode;                   /**< a element's nodes number*/
    static const int m_mQPoint;                 /**< num of quadrature points of a elmt*/
    static const int m_QPW;                     /**< weightness of quadrature points of a elmt*/
    static ShpfunQuad4 m_shpfun;                /**< shape function relative computer*/
};