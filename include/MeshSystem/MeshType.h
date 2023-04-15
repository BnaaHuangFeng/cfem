#pragma once

/**
 * the enum class for the different type of mesh supported in AsFem
 */
enum class MeshType{
    NULLTYPE,
    // for 0d case
    POINT,
    // for 1d case
    EDGE2,
    EDGE3,
    EDGE4,
    EDGE5,
    // for 2d case
    TRI3,
    TRI6,
    QUAD4,
    QUAD8,
    QUAD9,
    // for 3d case
    TET4,
    TET10,
    HEX8,
    HEX20,
    HEX27
};