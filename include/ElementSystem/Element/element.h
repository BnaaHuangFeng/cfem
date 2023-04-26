# pragma once
# include "petsc.h"
# include "MaterialSystem/Material.h"
# include "MathUtils/Vector2d.h"
class element{
protected:
    /**
     * the discrete symmetric gradient operation, B-matrix
     * @param t_dNdxPtr > ptr to derivate of shpfun to current coord, one item is for one node (Vector2d *, Vector3d *)
     * @param mNodePElmt > node num per elmt
     * @param mDofPNode > dof num per node
     * @param BMatrixPtr < ptr to recieve B-Matrix,need preallocation
    */ 
    bool getBMatrix(void *t_dNdxPtr, int mNodePElmt, int mDofPNode, MatrixXd *BMatrixPtr);
    /**
     * evaluates the discrete (full) gradient operation "G" component ordering (11,21,12,22)
     * @param t_dNdxPtr > ptr to derivate of shpfun to current coord, one item is for one node (Vector2d *, Vector3d *)
     * @param mNodePElmt > node num per elmt
     * @param mDofPNode > dof num per node
     * @param GMatrixPtr < ptr to recieve G Matrix,need preallocation
    */
    bool getGMatrix(void *t_dNdxPtr, int mNodePElmt, int mDofPNode, MatrixXd *GMatrixPtr);
public:
    element(bool nLarge):m_nLarge(nLarge),m_matPtr(nullptr){}
    element():m_elmt_rId(0),m_nLarge(false),m_matPtr(nullptr){};
    virtual ~element(){if(m_matPtr)delete m_matPtr;};
    /**
     * init element
     * @param matPtr > material ptr
     * @param nLarge > large strain flag
     * @param elmtParamPtr > ptr to the elmt's params
    */
    virtual PetscErrorCode initElement(PetscInt t_elmt_rId, bool nLarge, PetscScalar *elmtParamPtr)=0;
    /**
     * get the elmt's inner force
     * @param t_elmtCoord2 > ptr to the elmt's last converged coords (Vector2d *, Vector3d *)
     * @param t_elmtDofInc > ptr to the elmt's incremental dof values of this inrement until now (Vector2d *, Vector3d *)
     * @param t_elmtInnerForce < ptr to receive the elmt's inner force (need preallocation)
     * @param t_converged < if material update converged
    */
    virtual PetscErrorCode getElmtInnerForce(void *t_elmtCoord2, void *t_elmtDofInc, VectorXd *t_elmtInnerForce, bool *t_converged)=0;
    /**
     * get the elmt's stiffness matrix
     * @param t_elmtCoord2 > ptr to the elmt's last converged coords (Vector2d *, Vector3d *)
     * @param t_elmtDofInc > ptr to the elmt's incremental dof values of this inrement until now (Vector2d *, Vector3d *)
     * @param t_stfMatrix < ptr to receive the elmt's stiffness matrix (need preallocation)
    */
    virtual PetscErrorCode getElmtStfMatrix(void *t_elmtCoord2, void *t_elmtDofInc, MatrixXd *t_stfMatrix)=0;
    virtual int getDofNum()=0;
    virtual int getDim()=0;
    virtual int getNodeNum()=0;
    virtual int getDofPerNode()=0;
public:
    PetscInt    m_elmt_rId;     /**< elmt's id in rank*/
    bool        m_nLarge;       /**< large strain flag*/
    Material    *m_matPtr;      /**< material ptr*/
    
};