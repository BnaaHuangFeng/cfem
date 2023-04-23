#pragma once
#include "petsc.h"
#include <vector>
#include "InputSystem/DescriptionInfo.h"
#include "Utils/Timer.h"
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
    /**
     * output the mesh data
    */
    virtual PetscErrorCode outputMeshFile()=0;
    /**
     * init the mesh systems,including data preallocation and nodes' coords set
    */
    virtual PetscErrorCode MeshSystemInit(MeshDescription *t_meshDesPtr)=0;
    /**
     * get node's global id by id in rank
     * @param rId > node's id in this rank
    */
    /**
     * get coords of the nodes in a element by element's id in rank
     * @param elmtRId > elment's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （nodes id in elmt, dof id in node）-> coords 
     * @param nodeNum < ptr to store the node number of this element
    */
    virtual PetscErrorCode getElmtNodeCoord(PetscInt elmtRId,int state,PetscScalar **coordsPtr,PetscInt *nodeNum=nullptr)=0;
    /**
     * get coords of a node by its id in rank
     * @param nodeRId > node's id in rank
     * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
     * @param coordsPtr < ptr to store the node's coords （dof id in node）-> coords
    */
    virtual PetscErrorCode getNodeCoord(PetscInt nodeRId,int state,PetscScalar *coordsPtr)=0;
    inline int nodeRId2GId(int rId){return m_node_gId[rId];}
    // /**
    //  * get node's id in this rank by global id
    //  * @param gId > node's global id
    // */
    // virtual int nodeGId2RId(int gId)=0;
    /**
     * get element's global id by id in rank
     * @param rId > node's id in this rank
    */
    inline int elmtRId2GId(int rId){return m_elmt_gId[rId];}
    // /**
    //  * get element's id in this rank by global id
    //  * @param gId > node's global id
    // */
    // virtual int elmtGId2RId(int gId)=0;
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
public:
    bool m_ifSaveMesh;                      /**< if output the mesh to outer file*/
    std::string m_outputMeshFile_Name;
    std::string m_inputMeshFile_Name;
    PetscMPIInt m_rankNum;                  /**< processor num*/
    PetscMPIInt m_rank;                     /**< current rank*/
    PetscInt m_dim;                         /**< mesh dimension*/
    PetscInt m_mDof_node;                   /**< dof num per node*/
    PetscInt m_mNodes;                      /**< the sum of nodes in all processors*/
    PetscInt m_mNodes_p;                    /**< nodes num in current rank*/
    PetscInt m_mElmts;                      /**< the sum of elements in all processors*/
    PetscInt m_mElmts_p;                    /**< elements num in the current rank*/
    vector<PetscInt> m_node_gId;            /**< node's ID in rank -> node's global id*/
    vector<PetscInt> m_elmt_gId;            /**< element's ID in rank -> element's global ID*/
    vector<vector<PetscInt>> m_elmt_cnn;    /**< element's ID in rank -> element's global connectivity*/
    DM  m_dm;                               /**< Petsc DM*/
    Vec m_nodes_coord0;                     /**< global nodes' coords in ref config*/
    Vec m_nodes_coord1;                     /**< global nodes' coords in current config*/
    Vec m_nodes_coord2;                     /**< global nodes' coords of last converged config*/
    Timer *m_timerPtr;                      /**< clock ptr*/
};