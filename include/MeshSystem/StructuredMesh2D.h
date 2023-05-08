#pragma once
#include "MeshSystem/MeshSystem.h"
#include "petsc.h"
#include <string>
#include<fstream>
/**
 * this class store the topnology structure of the mesh,including node's ID. coords, element's connectivity
 * implement 2D structured mesh's init, supply query about mesh connectivity.
 */
using namespace std;
class StructuredMesh2D:public MeshSystem{
public:
    StructuredMesh2D();
    StructuredMesh2D(Timer *timerPtr);
//**********************************************************************************************
//** interface to creat mesh structure *********************************************************
//**********************************************************************************************/
    /**
     * init the mesh systems,including data preallocation and nodes' coords set
    */
    virtual PetscErrorCode MeshSystemInit(MeshDescription *t_meshDesPtr);
/**********************************************************************************************/

//**********************************************************************************************
//** interface to output mesh to outer file ****************************************************
//**********************************************************************************************
    /**
     * write mesh data to ouput mesh file
    */
    virtual PetscErrorCode outputMeshFile();
//**********************************************************************************************

//**********************************************************************************************
//** interface to Vec access control ***********************************************************
//**********************************************************************************************
    
    /**
     * open access to the Vec of coords (including create local Vec with values copy from global Vec)
     * @param vType > node variable type
     * @param variableVecPtr > ptr to corresponding golbal node variable Vec
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param mode > Vec access mode
    */
    virtual PetscErrorCode openNodeVariableVec(NodeVariableType vType, Vec *variableVecPtr, int state, VecAccessMode mode);
    /**
     * close access to the Vec of coords (including add local Vec's value to global Vec)
     * @param vType > node variable type
     * @param variableVecPtr > ptr to corresponding golbal node variable Vec
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param mode > Vec access mode
    */
    virtual PetscErrorCode closeNodeVariableVec(NodeVariableType vType, Vec *variableVecPtr, int state, VecAccessMode mode);
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
    virtual PetscErrorCode getElmtNodeCoord(PetscInt elmtRId,int state,Vector *coordsPtr,PetscInt *nodeNum=nullptr);
    /**
     * get coords of a node by its id in rank
     * @param nodeRId > node's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （dof id in node）-> coords
    */
    virtual PetscErrorCode getNodeCoord(PetscInt nodeRId,int state,Vector *coordsPtr);
    /**
     * get UInc of the nodes in a element by element's id in rank
     * @param elmtRId > elment's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param uIncPtr < ptr to store the node's UInc （nodes id in elmt, dof id in node）-> UInc 
     * @param nodeNum < ptr to store the node number of this element
    */
    virtual PetscErrorCode getElmtNodeUInc(PetscInt elmtRId,int state,Vector *uIncPtr,PetscInt *nodeNum=nullptr);
    /**
     * get residuals of the nodes in a element by element's id in rank 
     * @param elmtRId > elment's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param residualPtr < ptr to store the node's UInc （nodes id in elmt, dof id in node）-> UInc 
     * @param nodeNum < ptr to store the node number of this element
    */
    virtual PetscErrorCode getElmtNodeResidual(PetscInt elmtRId,int state,Vector *residualPtr,PetscInt *nodeNum=nullptr);
/**********************************************************************************************/

//**********************************************************************************************
//** interface to get rank local id by global id ***********************************************
//**********************************************************************************************   
    /**
     * get node's id in rank id by global id
     * @param gId > node's global id
    */    
    virtual int nodeGId2RId(int gId);
    /**
     * get elmt's id in rank id by global id
     * @param gId > elmt's global id
    */    
    virtual int elmtGId2RId(int gId); 

//**********************************************************************************************
//** interface to writting of data in mesh node ************************************************
//**********************************************************************************************
    /**
     * update mesh configuration set converged m_nodes_coord2 and  converged incremental u m_nodes_uInc2
     * @param snesPtr > ptr to SNES
    */
    virtual PetscErrorCode updateConfig(SNES *sensPtr);
    /**
     * Add a element's Jacobian (stiffness) matrix to global one (need to do MatAssembly after
     * elmts in this rank have called this func)
     * @param rId > elmt's id in this rank
     * @param matrixPtr >ptr to the elmt matrix to add
    */
    virtual PetscErrorCode addElmtAMatrix(PetscInt rid,MatrixXd *matrixPtr,Mat *APtr);
    /**
     * Add a element's residual (unbalanced forces (f^int-f^ext) ) Vector to global one (need to do 
     * closeNodeVariableVec after elmts in this rank have called this func in order to complete assembly)
     * @param rId > elmt's id in this rank
     * @param residualPtr >ptr to the elmt matrix to add (2 ind is node id in a elmt & dof id in a node)
    */
    virtual PetscErrorCode addElmtResidual(PetscInt rid,Vector *residualPtr, Vec *fPtr);
    /**
     * get the ref of the array to access node varible
     * @param vType > node variable type
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
    */
    PetscScalar *** & getNodeVariablePtrRef(NodeVariableType vType, int state);
/**********************************************************************************************/
//**********************************************************************************************
//** for general utility                       *************************************************
//**********************************************************************************************
    /**
     * for debug print
    */
    virtual PetscErrorCode printVaribale(NodeVariableType vType, Vec *variableVecPtr, int state, int comp);
private:
    /**
     * 
    */
    /**
     * get coords of the nodes in a element by elmt's x,y global index in dmda
     * @param xI > DMDA x index
     * @param yI > DMDA y index
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param variablePtr < ptr to store the node's variable （nodes id in elmt, dof id in node）-> coords 
     * @param nodeNum < ptr to store the node number of this element
    */
    PetscErrorCode getElmtNodeVariableByDmdaInd(NodeVariableType vType ,PetscInt xI,PetscInt yI,int state,Vector *variablePtr,PetscInt *nodeNum=nullptr);
    /**
     * get coords of a node in a element by node's x,y global index in dmda
     * @param vType > node variable type
     * @param xI > DMDA x index
     * @param yI > DMDA y index
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param variablePtr < ptr to store the node's variable （dof id in node）-> coords 
     * @param nodeNum < ptr to store the node number of this element
    */
    PetscErrorCode getNodeVariableByDmdaInd(NodeVariableType vType, PetscInt xI,PetscInt yI,int state,Vector *variablePtr);
    /**
     * Add a element's Jacobian (stiffness) matrix to global one by DMDA index (need to do MatAssembly after
     * elmts in this rank have called this func)
     * @param vType > node variable type
     * @param xI > DMDA x index
     * @param yI > DMDA y index
     * @param matrixPtr >ptr to the elmt matrix to add
     * @param APtr > ptr to global Jacobian Mat 
    */
    PetscErrorCode addElmtAMatrixByDmdaInd(PetscInt xI,PetscInt yI,MatrixXd *matrixPtr,Mat *APtr);  
    /**
     * Add a element's residual (unbalanced forces (f^int-f^ext) ) Vector to global one by DMDA index (need to do 
     * closeNodeVariableVec after elmts in this rank have called this func in order to complete assembly)
     * @param xI > DMDA x index
     * @param yI > DMDA y index
     * @param residualPtr >ptr to the elmt matrix to add (2 ind is node id in a elmt & dof id in a node)
     * @param fPtr > ptr to global residual Vec
    */
    PetscErrorCode addElmtResidualByDmdaInd(PetscInt xI,PetscInt yI,Vector *residualPtr,Vec *fPtr);  
    /**
     * get the ref of local Vec of node varible
     * @param vType > node variable type
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
    */
    Vec & getNodeLocalVariableVecRef(NodeVariableType vType, int state);
    /**
     * init a structured mesh
    */
    PetscErrorCode initStructuredMesh(MeshShape t_meshShape,MeshDescription *t_meshDesPtr);
    /**
     * Get the coords of rectangular mesh in ref config by DMDA index
     * @param xI > global x DMDA index
     * @param yI > global y DMDA index
     * @param aCoord < ptr to coord0 of a node
    */

    void getRectCoord0ByDmdaInd(PetscInt xI,PetscInt yI,Vector *aCoord);
    /**
     * Get the coords of sin mesh in ref config by DMDA index
     * @param xI > global x DMDA index
     * @param yI > global y DMDA index
     * @param aCoord < ptr to coord0 of a node
    */
    void getSinCoord0ByDmdaInd(PetscInt xI,PetscInt yI, Vector *aCoord);
    /**
     * Get the coords of half sin mesh in ref config by DMDA index
     * @param xI > global x DMDA index
     * @param yI > global y DMDA index
     * @param aCoord < ptr to coord0 of a node
    */
    void getHalfSinCoord0ByDmdaInd(PetscInt xI,PetscInt yI, Vector *aCoord);
    /**
     * Get the coords of cos mesh in ref config by DMDA index
     * @param xI > global x DMDA index
     * @param yI > global y DMDA index
     * @param aCoord < ptr to coord0 of a node
    */
    void getCosCoord0ByDmdaInd(PetscInt xI,PetscInt yI, Vector *aCoord);
    /**
     * Get the coords of half cos mesh in ref config by DMDA index
     * @param xI > global x DMDA index
     * @param yI > global y DMDA index
     * @param aCoord < ptr to coord0 of a node
    */
    void getHalfCosCoord0ByDmdaInd(PetscInt xI,PetscInt yI, Vector *aCoord);
    /**
     * Get element's global DMDA id via its id in rank
     * @param rId > elmt's id in rank
     * @param xIPtr > ptr to dmda x index
     * @param yIPtr > ptr to dmda y index
    */
    void getElmtDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr);
    /**
     * Get node's global DMDA id via its id in rank
     * @param rId > elmt's id in rank
     * @param xIPtr > ptr to dmda x index
     * @param yIPtr > ptr to dmda y index
    */
    void getNodeDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr);
    void openMeshOutputFile(ofstream *of,ios_base::openmode mode);
public:
    static const int vtkType=9;             /**< vtk cell type*/
    static const int m_mNode_elmt=4;        /**< node num per elmt*/
    DMDALocalInfo m_daInfo;                 /**< DMDA local info*/
    PetscScalar m_geoParam[3];              /**< geometry parameter to describe the regular mesh domain （for rectangular domain, is x length,y length in order; for sin or half sin domain, is span, amplitude, width in order)*/
    Vec m_nodes_coord0_local;               /**< local nodes' coords in ref config*/
    Vec m_nodes_coord2_local;               /**< local nodes' coords of last converged config*/
    Vec m_nodes_u2_local;                   /**< local nodes' u of last converged config*/
    Vec m_nodes_uInc1_local;                /**< local nodes' incremental displacement in current config*/
    Vec m_nodes_uInc2_local;                /**< local nodes' incremental displacement of last converged config*/
    Vec m_node_residual1_local;             /**< local nolinear function's residual Vec in current iteration, also unbalanced forces (f^int-f^ext)*/
    Vec m_node_residual2_local;             /**< local nolinear function's residual Vec in last converged increment, also unbalanced forces (f^int-f^ext)*/
    Vec m_node_load_local;                  /**< local outer load Vec*/
    PetscScalar ***m_array_nodes_coord0;    /**< ptr for access m_nodes_coord0_local*/
    PetscScalar ***m_array_nodes_coord2;    /**< ptr for access m_nodes_coord2_local*/
    PetscScalar ***m_array_nodes_uInc1;     /**< ptr for access m_nodes_uInc1_local*/
    PetscScalar ***m_array_nodes_uInc2;     /**< ptr for access m_nodes_uInc2_local*/
    PetscScalar ***m_array_nodes_u2;        /**< ptr for access m_nodes_uInc2_local*/
    PetscScalar ***m_array_nodes_residual1; /**< ptr for access m_node_residual1_local*/
    PetscScalar ***m_array_nodes_residual2; /**< ptr for access m_node_residual1_local*/
    PetscScalar ***m_array_nodes_load;      /**< ptr for access m_node_load_local*/
};