# pragma once
# include "petsc.h"
# include "ElementSystem/Shpfun/ShpfunQuad4.h"
# include "MaterialSystem/Material2D.h"
# include "ElementSystem/Element/element.h"
# include "MeshSystem/MeshSystem.h"
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
class CPE4R:public element{
    public:
    CPE4R():element(false),m_det_dx0dr(0.0){}
    CPE4R(bool nLarge):element(nLarge),m_det_dx0dr(0.0){}
    public:/**< inherent virtual func need to be implemented*/
    /**
     * init element
     * @param matPtr > material ptr
     * @param nLarge > large strain flag
     * @param elmtParamPtr > ptr to the elmt's params
    */
    virtual PetscErrorCode initElement(PetscInt t_elmt_rId, bool nLarge,MeshSystem *t_meshSysPtr,PetscScalar *elmtParamPtr);
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
    /**
     * get weighted volume quadrature of specific vector (wightness is shape function value in quadrature point)
     * @param t_valQPPtr > (qpoint id in a elmt, vector component id) -> vector value
     * @param t_valNodePtr < (node id in a elmt, vector component is) -> vector weighted quadrature value
     * @param t_mCpnt > components num of the vector
    */
    virtual PetscErrorCode getElmtWeightedVolumeInt(PetscScalar **t_valQPPtr,PetscScalar **t_valNodePtr, int t_mCpnt);
    /**
     * Get material variable of elmtVarType in PetscScalar array form
     * tensor of rank 2's vector order: 11 22 33 12 13 23
     * @param elmtVarType > required elemnt variable's type
     * @param elmtVarPtr < ptr to store the elemnt variable (1st ind is qpoint id in a elmt, 2nd ind is component) (need to preallocate)
    */
    virtual void getElmtVariableArray(ElementVariableType elmtVarType,PetscScalar **elmtVarPtr);
    virtual int getDofNum(){return m_mDof_node*m_mNode;}
    virtual int getDim(){return m_dim;}
    virtual int getNodeNum(){return m_mNode;}
    virtual int getDofPerNode(){return m_mDof_node;}
    virtual int getQpointNum(){return 1;};
    /**
     * get det of dx0/dr of i-th qpoint
     * @param i > qpoint id in a elmt 
    */
    virtual double getDetdx0dr(int i){
        if(i!=0){
            MessagePrinter::printRankError("Element of CPER4's quadrature point id must be 0");
            MessagePrinter::exitcfem();
        }
        return m_det_dx0dr;
    };
    public:
    double m_det_dx0dr;                         /**< det of dx0dr*/
    public:/**< static member (all elements of this kind share them)*/
    static const int m_dim;                     /**< element dimension*/
    static const int m_mDof_node;               /**< dof num per node*/
    static const int m_mNode;                   /**< a element's nodes number*/
    static const int m_mQPoint;                 /**< num of quadrature points of a elmt*/
    static const int m_QPW;                     /**< weightness of quadrature points of a elmt*/
    static ShpfunQuad4 m_shpfun;                /**< shape function relative computer*/
};