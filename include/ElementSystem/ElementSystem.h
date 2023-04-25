#pragma once
#include <vector>
#include "InputSystem/DescriptionInfo.h"
#include "ElementSystem/Element/element.h"
#include "MeshSystem/MeshSystem.h"
#include "petsc.h"
using namespace std;
class ElementSystem
{
public:
    ElementSystem(){};
    ElementSystem(ElementDescription *elmtDesPtr,MaterialDescription *matDesPtr);
    ~ElementSystem();
    PetscErrorCode init(ElementDescription *elmtDesPtr,MaterialDescription *matDesPtr);
    PetscErrorCode assignElmtType();
    PetscErrorCode assignMatType();
    PetscErrorCode assembleAMatrix(DM *dmPtr,Mat *AMatrixPtr);
private:
/***************************************************************************************************
 *  func ptr to get the nodes' variable of elmt calculation need                                 ***
***************************************************************************************************/
/**
 * get coords of the nodes in a element by element's id in rank
 * @param elmtRId > elment's id in rank
 * @param state > configuration for node's coords to get (0: ref config; 1: current config; 2: last converged)
 * @param coordsPtr < ptr to store the node's coords （nodes id in elmt, dof id in node）-> coords 
 * @param nodeNum < ptr to store the node number of this element
*/
PetscErrorCode (*m_getElmtNodeCoord_Prt)(PetscInt elmtRId,int state,PetscScalar **coordsPtr,PetscInt *nodeNum);

public:
/***************************************************************************************************
 *  petsc rank id                                                                                ***
***************************************************************************************************/
    PetscMPIInt m_rankNum;                  /**< processor num*/
    PetscMPIInt m_rank;                     /**< current rank*/    
/***************************************************************************************************
 *  elmt type & assigment description                                                            ***
***************************************************************************************************/
    vector<string> m_elmtTypeNames;                 /**< element name (used defined in input file)*/
    vector<ElementType> m_elmtTypes;                /**< element type*/
    vector<string> m_elmtAssignSetNames;            /**< names of set ralative element assign to*/
/***************************************************************************************************
 *  material type & assigment description                                                        ***
***************************************************************************************************/
    vector<string> m_matTypeNames;                  /**< material name (used defined in input file)*/
    vector<MaterialType> m_matTypes;                /**< material type*/
    vector<nlohmann::json> m_properties;            /**< material properties json type*/
    vector<string> m_materialAssignSetNames;        /**< names of set ralative material assign to*/
/***************************************************************************************************
 *  every elmt in this rank                                                                      ***
***************************************************************************************************/
    vector<element> m_elmts;                        /**< every elmt item in this rank*/
};
