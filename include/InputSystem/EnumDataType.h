#pragma once
enum Dimension{
    ONE,
    TWO,
    THREE
};
enum MeshMode{
    STRUCTURED,
    UNSTRUCTURED
};
enum MeshFileFormat{
    STRUCTURE,/**< structure mesh data generated by thus program*/
    MSH2,/**< GMSH mesh file of version 2.0 (.msh)*/
    MSH4,/**< GMSH mesh file of version 4.0 (.msh)*/
};
/**
 * for mesh type
*/
enum MeshType{
    QUAD4,/**< 4 node quad mesh*/
    QUAD8,
    HYBRID,
};
/**
 * for element type
*/
enum ElementType{
    CPE4R
};
enum MaterialType{
    LINEARELASTIC,
    NEOHOOKEAN,
    VONMISES
};
enum AlgorithmType{
    TANGENTNEWTON,
    ARCLENGTH_CYLENDER
};
enum FiledOutputFormat{
    VTU
};
enum HistoryOutputFormat{
    CSV
};
enum SetType{
    NODE,
    ELEMENT,
    BCELEMENT
};

enum FieldVariableType{
    VONMISES,
    STRESS,
    LOGSTRAIN,
    U,
};

enum HistoryVariableType{
    U1,U2,U3,
    RF1,RF2,RF3
};