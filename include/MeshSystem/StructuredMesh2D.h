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
    inline StructuredMesh2D():m_array_nodes_coord0(nullptr),m_array_nodes_coord1(nullptr),m_array_nodes_coord2(nullptr){};
    inline StructuredMesh2D(Timer *timerPtr):MeshSystem(timerPtr),m_array_nodes_coord0(nullptr),m_array_nodes_coord1(nullptr),m_array_nodes_coord2(nullptr){};
    virtual PetscErrorCode MeshSystemInit(MeshDescription *t_meshDesPtr);
    /**
     * get coords of the nodes in a element by element's id in rank
     * @param elmtRId > elment's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （nodes id in elmt, dof id in node）-> coords 
     * @param nodeNum < ptr to store the node number of this element
    */
    virtual PetscErrorCode getElmtNodeCoord(PetscInt elmtRId,int state,PetscScalar **coordsPtr,PetscInt *nodeNum=nullptr);
    /**
     * get coords of the nodes in a element by elmt's x,y global index in dmda
     * @param xI > DMDA x index
     * @param yI > DMDA y index
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （nodes id in elmt, dof id in node）-> coords 
     * @param nodeNum < ptr to store the node number of this element
    */
    PetscErrorCode getElmtNodeCoordByDmdaInd(PetscInt xI,PetscInt yI,int state,PetscScalar **coordsPtr,PetscInt *nodeNum=nullptr);
    /**
     * get coords of a node by node's id in rank
     * @param nodeRId > node's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （nodes id in elmt, dof id in node）-> coords 
    */
    PetscErrorCode getNodeCoord(PetscInt nodeRId,int state,PetscScalar *coordsPtr);
    /**
     * get coords of a node in a element by node's x,y global index in dmda
     * @param xI > DMDA x index
     * @param yI > DMDA y index
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （dof id in node）-> coords 
     * @param nodeNum < ptr to store the node number of this element
    */
    PetscErrorCode getNodeCoordByDmdaInd(PetscInt xI,PetscInt yI,int state,PetscScalar *coordsPtr);
    /**
     * write mesh data to ouput mesh file
    */
    virtual PetscErrorCode outputMeshFile();
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
    void getRectCoord0ByDmdaInd(PetscInt xI,PetscInt yI,PetscScalar *aCoord);
    /**
     * Get the coords of sin mesh in ref config by DMDA index
     * @param xI > global x DMDA index
     * @param yI > global y DMDA index
     * @param aCoord < ptr to coord0 of a node
    */
    void getSinCoord0ByDmdaInd(PetscInt xI,PetscInt yI, PetscScalar *aCoord);
    /**
     * Get the coords of half sin mesh in ref config by DMDA index
     * @param xI > global x DMDA index
     * @param yI > global y DMDA index
     * @param aCoord < ptr to coord0 of a node
    */
    void getHalfSinCoord0ByDmdaInd(PetscInt xI,PetscInt yI, PetscScalar *aCoord);
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
    /**
     * set coord array ptr to access coord Vec
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
    */
    PetscErrorCode coordVecGetArray(int state);
    /**
     * restore coord array
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
    */
    PetscErrorCode coordVecRestoreArray(int state);
public:
    static const int vtkType=9;             /**< vtk cell type*/
    DMDALocalInfo m_daInfo;                 /**< DMDA local info*/
    PetscScalar m_geoParam[3];              /**< geometry parameter to describe the regular mesh domain （for rectangular domain, is x length,y length in order; for sin or half sin domain, is span, amplitude, width in order)*/
    PetscScalar ***m_array_nodes_coord0;    /**< ptr for access m_nodes_coord0*/
    PetscScalar ***m_array_nodes_coord1;    /**< ptr for access m_nodes_coord1*/
    PetscScalar ***m_array_nodes_coord2;    /**< ptr for access m_nodes_coord2*/
};