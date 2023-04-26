#pragma once
#include <map>
#include <vector>
#include <string>
#include "petsc.h"
#include "InputSystem/EnumDataType.h"
using namespace std;
/**
 * this class manager FEM's set, a set is a vector<int> data to
 * store all FEM item's rid (id in this rank) in this set, if a set only has a
 * element and the element equal -1, then that mean this set contain
 * all item of corresponding FEM data type
*/
class SetManager{
public:
    SetManager();
    ~SetManager();
    /** create a new FEM items set*/
    bool createSet(string setName, SetType setType, vector<PetscInt> *setPtr);
    /** append new FEM items to existing set*/
    bool appendItems2Set(string setName, SetType setType, vector<PetscInt> *setPtr);
    /** push back a single FEM item to existing set*/
    bool pushItem2Set(string setName, SetType setType, PetscInt rId);
    bool renameSet(string oldName, SetType setType, string newName);
    bool deleteSet(string setName, SetType setType);
    /**
     * get a set by set name
     * @param setName > set name
     * @param setType > item in the set's FEM type
     * @param return < set ref
    */
    vector<PetscInt> & getSet(string setName, SetType setType);
public:
    map<string,vector<PetscInt>> m_nodeSets;    /**< node set name -> node set*/
    map<string,vector<PetscInt>> m_elmtSets;    /**< elmt set name -> elmt set*/
    map<string,vector<PetscInt>> m_bElmtSets;   /**< bounday elmt set name -> boundary elmt set*/
};