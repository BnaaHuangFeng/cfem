#include "Init/SystemInit.h"
#include "MeshSystem/StructuredMesh2D.h"
PetscErrorCode MeshSystemInit(Timer *timerPtr, MeshDescription *MeshDesPtr,MeshSystem *meshSysPtr){
    MeshMode meshMode=MeshDesPtr->s_mode;
    Dimension meshDim=MeshDesPtr->s_dim;
    switch(meshDim){
        case Dimension::ONE:
            MessagePrinter::printErrorTxt("1D FEM is not developed now");
            MessagePrinter::exitcfem();
            break;
        case Dimension::TWO:
            switch (meshMode)
            {
            case MeshMode::STRUCTURED:
                meshSysPtr = new StructuredMesh2D(timerPtr);
                meshSysPtr->MeshSystemInit(MeshDesPtr);
                break;
            case MeshMode::UNSTRUCTURED:
                MessagePrinter::printErrorTxt("unstructured FEM is not developed now");
                MessagePrinter::exitcfem();
                break;
            default:
                break;
            }
            break;
        case Dimension::THREE:
            MessagePrinter::printErrorTxt("3D FEM is not developed now");
            MessagePrinter::exitcfem();
            break;
    }
    if(meshSysPtr->m_ifSaveMesh)meshSysPtr->outputMeshFile();
    return 0;
}