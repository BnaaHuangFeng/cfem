#include "petsc.h"
#include "InputSystem/InputSystem.h"
#include "MeshSystem/StructuredMesh2D.h"
#include "Utils/MessagePrinter.h"
#include "Utils/Timer.h"
#include "Init/SystemInit.h"
#include "unistd.h"
int main(int args,char *argv[]){
    int *FLAG=new int;
    FLAG[0]=1;
    while(FLAG[0])sleep(2);
    PetscCall(PetscInitialize(&args,&argv,NULL,NULL));
    Timer timer;
    InputSystem inputSystem(&timer);
    /******************************************************/
    /** read the input file                             ***/
    /******************************************************/
    inputSystem.init(args,argv);
    inputSystem.readFile();
    /******************************************************/
    /** create all kinds of system's ptr                ***/
    /******************************************************/
    MeshSystem *meshSysPtr=nullptr;
    /******************************************************/
    /** init all system                                 ***/
    /******************************************************/
    // init mesh system
    MeshSystemInit(&timer,&inputSystem.m_meshDes,meshSysPtr);
    /******************************************************/
    /** delete the class created by new                 ***/
    /******************************************************/
    delete meshSysPtr;
    PetscCall(PetscFinalize());
    return 0;
}
