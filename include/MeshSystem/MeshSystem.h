#pragma once
#include "petsc.h"
#include <vector>
#include "InputSystem/DescriptionInfo.h"
#include "Utils/Timer.h"
#include "MathUtils/MatrixXd.h"
#include "MeshSystem/SetManager.h"
/**
 * this class store the topnology structure of the mesh,including node's ID. coords, element's connectivity
 */
using namespace std;
class MeshDescription;
class MeshSystem{
public:
    /**
     * construction of MeshSystem
    */
    MeshSystem();
    MeshSystem(Timer *timerPtr);
    /**
     * destruction of MeshSystem
    */
    virtual ~MeshSystem();
//**********************************************************************************************
//** interface to creat mesh structure *********************************************************
//**********************************************************************************************
    /**
     * init the mesh systems,including data preallocation and nodes' coords set
    */
    virtual PetscErrorCode MeshSystemInit(MeshDescription *t_meshDesPtr)=0;
/**********************************************************************************************/


//**********************************************************************************************
//** interface to output mesh to outer file ****************************************************
//**********************************************************************************************
    /**
     * output the mesh data
    */
    virtual PetscErrorCode outputMeshFile()=0;
/**********************************************************************************************/


//**********************************************************************************************
//** interface to Vec access control ***********************************************************
//**********************************************************************************************
    /**
     * open access to the Vec of coords
     * @param vType > node variable type
     * @param variableVecPtr > ptr to corresponding golbal node variable Vec
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param mode > Vec access mode
    */
    virtual PetscErrorCode openNodeVariableVec(NodeVariableType vType, Vec *variableVecPtr, int state, VecAccessMode mode)=0;
    /**
     * close access to the Vec of coords
     * @param vType > node variable type
     * @param variableVecPtr > ptr to corresponding golbal node variable Vec
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param mode > Vec access mode
    */
    virtual PetscErrorCode closeNodeVariableVec(NodeVariableType vType, Vec *variableVecPtr, int state, VecAccessMode mode)=0;
/**********************************************************************************************/


//**********************************************************************************************
//** interface to writting of data in mesh node ************************************************
//**********************************************************************************************
    /**
     * update mesh configuration set converged m_nodes_coord2 and  converged incremental u m_nodes_uInc2
     * @param snesPtr > ptr to SNES
    */
    virtual PetscErrorCode updateConfig(SNES *sensPtr)=0;
    /**
     * Add a element's Jacobian (stiffness) matrix to global one by elmt id in rank (need to do MatAssembly after
     * elmts in this rank have called this func)
     * @param rId > elmt's id in this rank
     * @param matrixPtr >ptr to the elmt matrix to add
     * @param APtr > ptr to the global jacobian matrix
    */
    virtual PetscErrorCode addElmtAMatrix(PetscInt rid,MatrixXd *matrixPtr,Mat *APtr)=0;
    /**
     * Add a element's residual (unbalanced forces (f^int-f^ext) ) Vector to global one by elmt id in rank (need to do 
     * closeNodeVariableVec after elmts in this rank have called this func in order to complete assembly)
     * @param rId > elmt's id in this rank
     * @param residualPtr >ptr to the elmt matrix to add (2 ind is node id in a elmt & dof id in a node)
     * @param fPtr > ptr to the global residual matrix
    */
    virtual PetscErrorCode addElmtResidual(PetscInt rid,Vector *residualPtr, Vec *fPtr)=0;
/**********************************************************************************************/

//**********************************************************************************************
//** interface to reading of data in mesh node *************************************************
//**********************************************************************************************
    /**
     * get coords of the nodes in a element by element's id in rank
     * @param elmtRId > elment's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （nodes id in elmt, dof id in node）-> coords 
     * @param nodeNum < ptr to store the node number of this element
    */
    virtual PetscErrorCode getElmtNodeCoord(PetscInt elmtRId,int state,Vector *coordsPtr,PetscInt *nodeNum=nullptr)=0;
    /**
     * get coords of a node by its id in rank
     * @param nodeRId > node's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （dof id in node）-> coords
    */
    virtual PetscErrorCode getNodeCoord(PetscInt nodeRId,int state,Vector *coordsPtr)=0;
    /**
     * get UInc of the nodes in a element by element's id in rank
     * @param elmtRId > elment's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param uIncPtr < ptr to store the node's UInc （nodes id in elmt, dof id in node）-> UInc 
     * @param nodeNum < ptr to store the node number of this element
    */
    virtual PetscErrorCode getElmtNodeUInc(PetscInt elmtRId,int state,Vector *uIncPtr,PetscInt *nodeNum=nullptr)=0;
    /**
     * get residuals of the nodes in a element by element's id in rank
     * @param elmtRId > elment's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param residualPtr < ptr to store the node's UInc （nodes id in elmt, dof id in node）-> UInc 
     * @param nodeNum < ptr to store the node number of this element
    */
    virtual PetscErrorCode getElmtNodeResidual(PetscInt elmtRId,int state,Vector *residualPtr,PetscInt *nodeNum=nullptr)=0;
//**********************************************************************************************
//** for general utility                       *************************************************
//**********************************************************************************************

//**********************************************************************************************
//** for debug                                 *************************************************
//**********************************************************************************************
    virtual PetscErrorCode printVaribale(NodeVariableType vType, Vec *variableVecPtr, int state, int comp)=0;
public:
//**********************************************************************************************
//** Mesh system's commom implement ************************************************************
//**********************************************************************************************
    inline int dofRId2GId(int rId){return m_dof_gId[rId];}
    inline int nodeRId2GId(int rId){return m_node_gId[rId];}
    /**
     * get element's global id by id in rank
     * @param rId > node's id in this rank
    */
    inline int elmtRId2GId(int rId){return m_elmt_gId[rId];}
    /**
     * get node's id in rank id by global id
     * @param gId > node's global id
    */    
    virtual int nodeGId2RId(int gId)=0;
    /**
     * get elmt's id in rank id by global id
     * @param gId > elmt's global id
    */    
    virtual int elmtGId2RId(int gId)=0;    
    /**
     * get element's connectivity by elmt id in this rank.
     * @param rId > elmt's id in rank
     * @param elmtCnn < ptr 2 elmt's connectivity
    */
    bool getElmtCnnByRId(int rId, const vector<int> *elmtCnn);
    /**
     * check if node id in rank is in range, if not then print the error message,
     * and exit cfem
     * @param rId > node id in rank
    */   
    void checkNodeRId(int rId);
    /**
     * check if elmt id in rank is in range, if not then print the error message,
     * and exit cfem
     * @param rId > element id in rank
    */
    void checkElmtRId(int rId);
    /**
     * return the global Vec ptr corresponding to the NodeVariableType
    */
    Vec *globalVecPtr(NodeVariableType vType,int state);
    /**
     * 
    */
/**********************************************************************************************/
public:
    bool m_ifSaveMesh;                      /**< if output the mesh to outer file*/
    std::string m_outputMeshFile_Name;
    std::string m_inputMeshFile_Name;
    PetscMPIInt m_rankNum;                  /**< processor num*/
    PetscMPIInt m_rank;                     /**< current rank*/
    MeshMode m_meshMode;                    /**< mesh mode*/
    PetscInt m_dim;                         /**< mesh dimension*/
    PetscInt m_mDof_node;                   /**< dof num per node*/
    PetscInt m_mNodes;                      /**< the sum of nodes in all processors*/
    PetscInt m_mNodes_p;                    /**< nodes num in current rank*/
    PetscInt m_mElmts;                      /**< the sum of elements in all processors*/
    PetscInt m_mElmts_p;                    /**< elements num in the current rank*/
    vector<PetscInt> m_dof_gId;             /**> dof's ID in rank -> dof's global id*/
    vector<PetscInt> m_node_gId;            /**< node's ID in rank -> node's global id*/
    vector<PetscInt> m_elmt_gId;            /**< element's ID in rank -> element's global ID*/
    vector<vector<PetscInt>> m_elmt_cnn;    /**< element's ID in rank -> element's global connectivity*/
    Timer *m_timerPtr;                      /**< clock ptr*/
    SetManager m_setManager;                /**< FEM set manager*/
    DM  m_dm;                               /**< Petsc DM*/
    Vec m_nodes_coord0;                     /**< global nodes' coords in ref config*/
    Vec m_nodes_coord2;                     /**< global nodes' coords of last converged config*/
    Vec m_nodes_u2;                         /**< global nodes' coords of last converged config*/
    /**< m_nodes_uInc1 is managered by SNES*/
    Vec m_nodes_uInc2;                      /**< global nodes' incremental displacement in current config*/
    /**< m_AMatrix1 is managered by SNES*/
    Mat m_AMatrix2;                         /**< nolinear function's jacobian matrix, also tangent stiffness matrix*/
    /**< m_node_residual2 is managered by SNES*/
    Vec m_node_residual2;                   /**< nolinear function's residual Vec, also unbalanced forces (f^int-f^ext)*/
    Vec m_node_load;                        /**< node outer load Vec*/
};