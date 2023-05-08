#include "petsc.h"
#include "InputSystem/InputSystem.h"
#include "MeshSystem/StructuredMesh2D.h"
#include "Utils/MessagePrinter.h"
#include "Utils/Timer.h"
#include "Init/SystemInit.h"
#include "unistd.h"
int main(int args,char *argv[]){
    int *FLAG=new int;
    // FLAG[0]=1;
    // while(FLAG[0])sleep(2);
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
    timer.startTimer();
    MessagePrinter::printDashLine(MessageColor::BLUE);
    MessagePrinter::printNormalTxt("Start to dynamically create and init every system and controller",MessageColor::BLUE);    
    MeshSystem *meshSysPtr=nullptr;
    ElementSystem *elmtSysPtr=nullptr;
    BCsSystem *BCsSysPtr=nullptr;
    LoadController *loadCtrlPtr=nullptr;
    SolutionSystem *solSysPtr=nullptr;
    PostProcessSystem *postSysPtr=nullptr;
    /******************************************************/
    /** init all system                                 ***/
    /******************************************************/
    MeshSystemInit(&timer,&inputSystem.m_meshDes,&meshSysPtr);
    ElmtSystemInit(&timer,&elmtSysPtr,&inputSystem.m_ElDes,&inputSystem.m_MatDes,meshSysPtr);
    BCsSystemInit(&BCsSysPtr,&inputSystem.m_bcDes,meshSysPtr);
    LoadCtrolInit(&loadCtrlPtr,&inputSystem.m_stepDes,meshSysPtr);
    SolutionSysInit(&solSysPtr,&inputSystem.m_stepDes,meshSysPtr,elmtSysPtr,BCsSysPtr,loadCtrlPtr);
    PostSysInit(&postSysPtr,&inputSystem.m_outDes,meshSysPtr,elmtSysPtr,loadCtrlPtr);
    timer.endTimer();
    timer.printElapseTime("system and controller inition is done",false);
    MessagePrinter::printNormalTxt("All system and controller inition completed!",MessageColor::BLUE);  
    MessagePrinter::printDashLine(MessageColor::BLUE); 
    //output intial state
    postSysPtr->output(solSysPtr->m_increI,loadCtrlPtr->m_factor2);
    ++(solSysPtr->m_increI);
    bool ifConverged=true, ifCompleted=false;
    while(!ifCompleted){
        solSysPtr->run(ifConverged,&ifConverged,&ifCompleted);
        if(ifConverged&&!ifCompleted){
            meshSysPtr->updateConfig(&solSysPtr->m_snes);
            elmtSysPtr->updateConvergence();
            postSysPtr->output(solSysPtr->m_increI,loadCtrlPtr->m_factor1);
            ++(solSysPtr->m_increI);
        }
    }
    // StructuredMesh2D *meshSysPtr2=(StructuredMesh2D *)meshSysPtr;
    // meshSysPtr2->printVaribale(NodeVariableType::U,&meshSysPtr2->m_nodes_u2,2,0);
    // meshSysPtr2->printVaribale(NodeVariableType::U,&meshSysPtr2->m_nodes_u2,2,1);
    /******************************************************/
    /** delete the class created by new                 ***/
    /******************************************************/
    if(meshSysPtr) delete meshSysPtr;
    if(elmtSysPtr) delete elmtSysPtr;
    if(BCsSysPtr) delete BCsSysPtr;
    if(loadCtrlPtr) delete loadCtrlPtr;
    if(solSysPtr) delete solSysPtr;
    postSysPtr->clear();
    // if(postSysPtr) delete postSysPtr;
    PetscCall(PetscFinalize());
    delete FLAG;
    return 0;
}
