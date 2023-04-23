#include "petsc.h"
// #include "InputSystem/InputSystem.h"
// #include "MeshSystem/StructuredMesh2D.h"
// #include "Utils/MessagePrinter.h"
// #include "Utils/Timer.h"
int main(int args,char *argv[]){
    PetscCall(PetscInitialize(&args,&argv,NULL,NULL));
    // Timer timer;
    // InputSystem inputSystem(&timer);
    // /******************************************************/
    // /** read the input file                             ***/
    // /******************************************************/
    // inputSystem.init(args,argv);
    // inputSystem.readFile();
    // /******************************************************/
    // /** create all kinds of system's ptr                ***/
    // /******************************************************/
    // MeshSystem *meshPtr=nullptr;
    // /******************************************************/
    // /** init all system                                 ***/
    // /******************************************************/
    // // init mesh system
    // MeshMode meshMode=inputSystem.m_meshDes.s_mode;
    // Dimension meshDim=inputSystem.m_meshDes.s_dim;
    // switch(meshDim){
    //     case Dimension::ONE:
    //         MessagePrinter::printErrorTxt("1D FEM is not developed now");
    //         MessagePrinter::exitcfem();
    //         break;
    //     case Dimension::TWO:
    //         switch (meshMode)
    //         {
    //         case MeshMode::STRUCTURED:
    //             meshPtr = new StructuredMesh2D(&timer);
    //             meshPtr->MeshSystemInit(&inputSystem.m_meshDes);
    //             break;
    //         case MeshMode::UNSTRUCTURED:
    //             MessagePrinter::printErrorTxt("unstructured FEM is not developed now");
    //             MessagePrinter::exitcfem();
    //             break;
    //         default:
    //             break;
    //         }
    //         break;
    //     case Dimension::THREE:
    //         MessagePrinter::printErrorTxt("3D FEM is not developed now");
    //         MessagePrinter::exitcfem();
    //         break;
    // }
    // if(meshPtr->m_ifSaveMesh)meshPtr->outputMeshFile();
    // /******************************************************/
    // /** delete the class vreated by new                 ***/
    // /******************************************************/
    // delete meshPtr;
    PetscCall(PetscFinalize());
    return 0;
}
