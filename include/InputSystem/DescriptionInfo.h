/**
 * This file contain all kinds of data structure to describe the FEM problem
 * including mesh, elements, step, nlsolver, output, bcs
*/
#pragma once
#include"InputSystem/EnumDataType.h"
#include<vector>
#include<string>
#include"petsc.h"
#include"Utils/JsonUtils.h"
#include"nlohmann/json.hpp"
# define ARCLENGTH_CYLENDER "arclength_cylender"
using namespace std;
struct MeshDescription{
    MeshMode s_mode;    /**< structured, unstructured,...*/
    MeshType s_type;    /**< quad4, quad8, ...*/
    Dimension s_dim;    /**< 1D, 2D, 3D*/
    MeshShape s_shape;  /**< the complete mesh outer shape*/
    int s_nx,s_ny,s_nz; /**< element's number in x/y/z diretion for structure mesh*/
    nlohmann::json s_size_json;/**< the complete mesh domain's outer boundary shape size*/
    string s_outputMeshFile_Name;
    string s_inputMeshFile_Name;
    bool s_ifSaveMesh;  /**< if output the mesh data*/
};
struct ElementDescription
{
    bool s_nLarge;                      /**< if strain large?*/    
    vector<string> s_names;             /**< element name (used defined in input file)*/
    vector<ElementType> s_elmtTypes;    /**< element type*/
    vector<string> s_setNames;          /**< names of set ralative element assign to*/
};
struct MaterialDescription{
    vector<string> s_names;                 /**< material name (used defined in input file)*/
    vector<MaterialType> s_matType;         /**< material type*/
    vector<nlohmann::json> s_properties;    /**< material properties json type*/
    vector<string> s_setName;               /**< names of set ralative material assign to*/
};
struct StepDescriptiom
{
    SNESType s_SNESType;        /**< Petsc nolinear solver SNES type*/
    KSPType s_KSPType;          /**< Petsc linear solver KSP type*/
    PCType s_PCType;            /**< Petsc precondition type*/
    int s_maxIterNum;           /**< max iteration number*/
    double s_t;                 /**< total time*/
    double s_dt0;               /**< first time increment*/
    double s_dtmax;             /**< max time increment*/
    double s_dtmin;             /**< min time increment*/
    double s_growFactor;        /**< factor growth ratio*/
    double s_cutbackFactor;     /**< factor cutback ratio*/
    double s_absTol; /**< absolute tolerance*/
    double s_relTol; /**< realative tolerance*/
    double s_duTol; /**< delta U tolerance*/
};
struct FieldOutputDescription{
    FiledOutputFormat s_format;
    int s_interval;
    vector<FieldVariableType> s_varTypes;
};
struct HistoryOutputDescription
{
    string s_setName;
    HistoryOutputFormat s_format;
    int s_interval;
    vector<HistoryVariableType> s_varTypes;
};
struct OutputDescription
{
    FieldOutputDescription s_FD; /**< field Output Description*/
    vector<HistoryOutputDescription> s_HD; /**< History Output Description*/
};
/**
 * Single Boundary condition description
*/
struct SingleBCDes{
    string s_BCName;
    vector<int> s_presetDofIds; /**< each item is a dof kind's vec of a BC assignment*/
    double s_bcVals; /**< each item is for a BC assignment's dof preset val*/
    string s_setName; /**< set name BCs assign to*/    
};
/**
 * Boundary condition description
*/
typedef vector<SingleBCDes> BCDescription;