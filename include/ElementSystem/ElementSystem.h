#pragma once
#include <vector>
#include "InputSystem/DescriptionInfo.h"
#include "ElementSystem/Element/element.h"
#include "MeshSystem/MeshSystem.h"
#include "petsc.h"
using namespace std;
class ElementSystem
{
protected:
    /**
     * read elmt description
    */
    void readElmtDes(ElementDescription *elmtDesPtr);
    /**
     * read mat description
    */
    void readMatDes(MaterialDescription *matDesPtr);
    /**
     * create elmt item and put their ptr in m_elmts member
     * @param meshPtr > ptr to the based mesh system
    */
    PetscErrorCode assignElmtType(MeshSystem *meshPtr);
    /**
     * create mat item for every elmt, and direct every elmt's m_matPtr to the mat item
     * @param meshPtr > ptr to the based mesh system
    */
    PetscErrorCode assignMatType(MeshSystem *meshPtr);
public:
    ElementSystem();
    ElementSystem(Timer* timerPtr,ElementDescription *elmtDesPtr,MaterialDescription *matDesPtr);
    ~ElementSystem();
    /**
     * init ElmentSystem, including read elmt and material description, preallocate elmt item ptr's vector,
     * new elmt item and its material item
     * @param meshPtr > ptr to the based mesh system
    */
    PetscErrorCode init(ElementDescription *elmtDesPtr,MaterialDescription *matDesPtr,MeshSystem *meshPtr);
    PetscErrorCode init(MeshSystem *meshPtr);
    /**
     * check if every elmts in this rank has specify elmt type and material type.
    */
    bool checkElmtsAssigment(MeshSystem *meshPtr);
    /**
     * assemble global Jacobian matrix for Petsc's solver
     * @param t_meshPtr > ptr to the based mesh system
     * @param t_AMatrixPtr > ptr to global Jacobian matrix to assemble
    */
    PetscErrorCode assembleAMatrix(MeshSystem *t_meshPtr,Mat *t_AMatrixPtr);
    /**
     * assemble global residual Vec for Petsc's solver
     * @param t_meshPtr > ptr to the based mesh system
     * @param t_RVecPtr > ptr to global residual Vec to assemble
    */
    PetscErrorCode assemblRVec(MeshSystem *t_meshPtr, Vec *t_RVecPtr);
protected:
    Timer *m_timerPtr;                      /**< clock ptr*/
/***************************************************************************************************
 *  petsc rank id                                                                                ***
***************************************************************************************************/
    PetscMPIInt m_rankNum;                  /**< processor num*/
    PetscMPIInt m_rank;                     /**< current rank*/    
public:
/***************************************************************************************************
 *  elmt type & assigment description                                                            ***
***************************************************************************************************/
    bool m_hasElmtDes;                              /**< if has read the elmt description?*/
    vector<string> m_elmtTypeNames;                 /**< element name (used defined in input file)*/
    vector<ElementType> m_elmtTypes;                /**< element type*/
    vector<string> m_elmtAssignSetNames;            /**< names of set ralative element assign to*/
    bool m_nLarge;                                  /**< if strain large*/
/***************************************************************************************************
 *  material type & assigment description                                                        ***
***************************************************************************************************/
    bool m_hasMatDes;                               /**< if has read the mat description*/
    vector<string> m_matTypeNames;                  /**< material name (used defined in input file)*/
    vector<MaterialType> m_matTypes;                /**< material type*/
    vector<nlohmann::json> m_properties;            /**< material properties json type*/
    vector<string> m_materialAssignSetNames;        /**< names of set ralative material assign to*/
/***************************************************************************************************
 *  every elmt in this rank                                                                      ***
***************************************************************************************************/
    vector<element *> m_elmtPtrs;                   /**< every elmt item in this rank*/
};