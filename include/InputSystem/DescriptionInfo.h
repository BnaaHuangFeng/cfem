/**
 * This file include all kinds of data structure to describe the FEM problem
 * including mesh, elements, step, nlsolver, output, bcs
*/
#pragma once
#include"InputSystem/EnumDataType.h"
#include<vector>
#include<string>
#include"petsc.h"
using namespace std;
struct MeshDescription{
    MeshMode s_mode;    /**< structured, unstructured,...*/
    MeshType s_type;    /**< quad4, quad8, ...*/
    Dimension s_dim;    /**< 1D, 2D, 3D*/
    int s_nx,s_ny,s_nz; /**< element's number in x/y/z diretion for structure mesh*/
    double s_lx,s_ly,l_lz;/**< structure's length in x/y/z direction for rectangle domain*/
    bool s_ifSaveMesh;  /**< if output the mesh data*/
};
struct ElementDescription
{
    vector<string> s_Names; /**< element name (used difined in input file)*/
    vector<ElementType> s_elmtTypes; /** element type*/
    vector<MaterialType> s_MatType; /**< material type*/
    vector<vector<double>> s_properties; /**< material properties*/
};
struct StepDescriptiom
{
    AlgorithmType s_algType; /**< FEM iteration algorithm type*/
    SNESType s_SNESType; /**< Petsc nolinear solver SNES type*/
    PCType s_PCType; /**< Petsc precondition type*/
    int s_maxIterNum; /**< max iteration number*/
    double s_dto; /**< first time increment*/
    double s_dtmax; /**< max time increment*/
    double s_dtmin; /**< min time increment*/
    double s_absTol; /**< absolute tolerance*/
    double s_relTol; /**< realative tolerance*/
    double s_duTol; /**< delta U tolerance*/
};
struct FieldOutputDescription{
    FiledOutputFormat s_format;
    int interval;
    vector<FieldVariableType> s_varTypes;
};
struct HistoryOutputDescription
{
    string setName;
    HistoryOutputFormat s_format;
    int interval;
    vector<HistoryVariableType> s_varType;
};
struct OutputDescription
{
    FieldOutputDescription s_FD; /**< field Output Description*/
    HistoryOutputDescription s_HD; /**< History Output Description*/
};
/**
 * Boundary condition description
*/
struct BCDescription
{
    vector<string> s_BCName;
    vector<vector<int>> s_presetDofIds; /**< each item is a dof kind's vec of a BC assignment*/
    vector<double> s_bcVals; /**< each item is for a BC assignment's dof preset val*/
    vector<string> s_setName; /**< set name BCs assign to*/
};
