#pragma once
enum class Dimension{
    ONE,
    TWO,
    THREE
};
enum class MeshShape{
    RECTANGULAR,
    SIN,
    HALFSIN,
    COMPLEX
};
enum class MeshMode{
    STRUCTURED,
    UNSTRUCTURED
};
enum class MeshFileFormat{
    STRUCTURE,/**< structure mesh data generated by thus program*/
    MSH2,/**< GMSH mesh file of version 2.0 (.msh)*/
    MSH4,/**< GMSH mesh file of version 4.0 (.msh)*/
};
/**
 * for mesh type
*/
enum class MeshType{
    QUAD4,/**< 4 node quad mesh*/
    QUAD8,
    HYBRID,
};
/**
 * for element type
*/
enum class ElementType{
    CPE4R
};
enum class MaterialType{
    LINEARELASTIC,
    NEOHOOKEAN,
    VONMISESPLAS
};
// enum class AlgorithmType{
//     TANGENTNEWTON,
//     ARCLENGTH_CYLENDER
// };
enum class FiledOutputFormat{
    VTU
};
enum class HistoryOutputFormat{
    CSV
};
enum class SetType{
    NODE,
    ELEMENT,
    BCELEMENT
};

enum class FieldVariableType{
    VONMISES,
    STRESS,
    LOGSTRAIN,
    U,
};

enum class HistoryVariableType{
    U1,U2,U3,
    RF1,RF2,RF3
};
enum class ElementVariableType{
    VONMISES,
    CAUCHYSTRESS,
    LOGSTRAIN,
    PRESSURE,
    KIRCHOFFSTRESS,
    JACOBIAN
};
enum class NodeVariableType{
    COORD,
    U,
    UINC,
    RF,
    RESIDUAL
};
enum class VecAccessMode{
    WRITE,
    READ
};