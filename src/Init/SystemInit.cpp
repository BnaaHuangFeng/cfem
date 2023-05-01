#include "Init/SystemInit.h"
#include "MeshSystem/StructuredMesh2D.h"
#include "BCsSystem/BCsSysStructured2d.h"
PetscErrorCode MeshSystemInit(Timer *timerPtr, MeshDescription *meshDesPtr,MeshSystem **meshSysPtrAdr){
    MeshMode meshMode=meshDesPtr->s_mode;
    Dimension meshDim=meshDesPtr->s_dim;
    switch(meshDim){
        case Dimension::ONE:
            MessagePrinter::printErrorTxt("1D FEM is not developed now");
            MessagePrinter::exitcfem();
            break;
        case Dimension::TWO:
            switch (meshMode)
            {
            case MeshMode::STRUCTURED:
                *meshSysPtrAdr = new StructuredMesh2D(timerPtr);
                (*meshSysPtrAdr)->MeshSystemInit(meshDesPtr);
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
    if((*meshSysPtrAdr)->m_ifSaveMesh)(*meshSysPtrAdr)->outputMeshFile();
    MessagePrinter::printNormalTxt("Mesh system inition is done");
    return 0;
}
PetscErrorCode ElmtSystemInit(Timer *timerPtr, ElementSystem **elmtSysPtrAdr, ElementDescription *elmtDesPtr,MaterialDescription *matDesPtr, MeshSystem *meshSysPtr){
    *elmtSysPtrAdr=new ElementSystem(timerPtr,elmtDesPtr,matDesPtr);
    (*elmtSysPtrAdr)->init(meshSysPtr);
    (*elmtSysPtrAdr)->checkInit();
    MessagePrinter::printNormalTxt("Element system inition is done");
    return 0;
}
PetscErrorCode BCsSystemInit(BCsSystem **BCsSysPtrAdr, BCDescription *BCDesPtr, MeshSystem *meshSysPtr){
    int dim=meshSysPtr->m_dim;
    MeshMode meshMode=meshSysPtr->m_meshMode;
    if(dim==2){
        switch (meshMode)
        {
        case MeshMode::STRUCTURED:
            *BCsSysPtrAdr=new BCsSysStructured2d(BCDesPtr,meshSysPtr);
            (*BCsSysPtrAdr)->init();
            (*BCsSysPtrAdr)->checkInit();
            break;
        case MeshMode::UNSTRUCTURED:
            MessagePrinter::printErrorTxt("unstructured element system is not developed now");
            MessagePrinter::exitcfem();
            break;              
        default:
            break;
        }
    }
    else{
        MessagePrinter::printErrorTxt(to_string(dim)+" dimensional element system is not developed now");
        MessagePrinter::exitcfem();        
    }
    MessagePrinter::printNormalTxt("Boundary condition system inition is done");
    return 0;
}
PetscErrorCode LoadCtrolInit(LoadController **loadCtrlPtrAdr, StepDescriptiom *stepDesPtr, MeshSystem *meshSysPtr){
    *loadCtrlPtrAdr=new LoadController(stepDesPtr,meshSysPtr);
    (*loadCtrlPtrAdr)->checkInit();
    MessagePrinter::printNormalTxt("Load controller inition is done");
    return 0;
}
PetscErrorCode SolutionSysInit(SolutionSystem **solutionSysPtrAdr, StepDescriptiom *stepDesPtr, MeshSystem *meshSysPtr,
                            ElementSystem *elmtSysPtr, BCsSystem *BCsSysPtr, LoadController *loadCtrlPtr){
    *solutionSysPtrAdr=new SolutionSystem(stepDesPtr,meshSysPtr);
    (*solutionSysPtrAdr)->init(stepDesPtr,meshSysPtr,elmtSysPtr,BCsSysPtr,loadCtrlPtr);
    (*solutionSysPtrAdr)->checkInit();
    MessagePrinter::printNormalTxt("Solution system inition is done");
    return 0;
}