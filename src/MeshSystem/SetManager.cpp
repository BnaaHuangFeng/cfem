#include "MeshSystem/SetManager.h"
#include "Utils/MessagePrinter.h"
SetManager::SetManager(){
}
SetManager::~SetManager(){
}
bool SetManager::createSet(string setName, SetType setType, vector<PetscInt> *setPtr){
    switch (setType)
    {
    case SetType::NODE:
        m_nodeSets.insert(pair<string,vector<PetscInt>>(setName,*setPtr));
        break;
    case SetType::ELEMENT:
        m_elmtSets.insert(pair<string,vector<PetscInt>>(setName,*setPtr));
        break;
    case SetType::BELEMENT:
        m_bElmtSets.insert(pair<string,vector<PetscInt>>(setName,*setPtr));
        break;
    default:
        MessagePrinter::printErrorTxt("createSet error: unsupported set type.");
        MessagePrinter::exitcfem();
        break;
    }
    return true;
}
bool SetManager::appendItems2Set(string setName, SetType setType, vector<PetscInt> *setPtr){
    map<string,vector<PetscInt>>::iterator it;
    map<string,vector<PetscInt>> *mapPtr;
    switch (setType)
    {
    case SetType::NODE:
        mapPtr=&m_nodeSets;
        break;
    case SetType::ELEMENT:
        mapPtr=&m_elmtSets;
        break;
    case SetType::BELEMENT:
        mapPtr=&m_bElmtSets;
        break;
    default:
        MessagePrinter::printErrorTxt("appendItems2Set error: unsupported set type.");
        MessagePrinter::exitcfem();        
        break;
    }
    it=mapPtr->find(setName);
    if(it==mapPtr->end()){
        MessagePrinter::printErrorTxt("appendItems2Set error: can not append item to a unexisting set");
        MessagePrinter::exitcfem();
    }
    else{
        it->second.insert(it->second.end(),setPtr->begin(),setPtr->end());
    }
    return true;
}
bool SetManager::pushItem2Set(string setName, SetType setType, PetscInt rId){
    map<string,vector<PetscInt>>::iterator it;
    map<string,vector<PetscInt>> *mapPtr;
    switch (setType)
    {
    case SetType::NODE:
        mapPtr=&m_nodeSets;
        break;
    case SetType::ELEMENT:
        mapPtr=&m_elmtSets;
        break;
    case SetType::BELEMENT:
        mapPtr=&m_bElmtSets;
        break;
    default:
        MessagePrinter::printErrorTxt("pushItem2Set error: unsupported set type.");
        MessagePrinter::exitcfem();        
        break;
    }
    it=mapPtr->find(setName);
    if(it==mapPtr->end()){
        MessagePrinter::printErrorTxt("pushItem2Set error: can not push item to a unexisting set");
        MessagePrinter::exitcfem();
    }
    else{
        it->second.push_back(rId);
    }
    return true;    
}
bool SetManager::renameSet(string oldName, SetType setType, string newName){
    map<string,vector<PetscInt>>::iterator it;
    map<string,vector<PetscInt>> *mapPtr;
    switch (setType)
    {
    case SetType::NODE:
        mapPtr=&m_nodeSets;
        break;
    case SetType::ELEMENT:
        mapPtr=&m_elmtSets;
        break;
    case SetType::BELEMENT:
        mapPtr=&m_bElmtSets;
        break;
    default:
        MessagePrinter::printErrorTxt("renameSet error: unsupported set type.");
        MessagePrinter::exitcfem();        
        break;
    }
    it=mapPtr->find(newName);
    if(it==mapPtr->end()){
        map<string,vector<PetscInt>>::iterator it2;
        it2=mapPtr->find(oldName);
        vector<PetscInt> set;
        if(it2==mapPtr->end()){
            MessagePrinter::printErrorTxt("renameSet error: set "+oldName+" unexists.");
            MessagePrinter::exitcfem();                    
        }
        else{
            set=it2->second;
            createSet(newName,setType,&set);
            deleteSet(oldName,setType);
        }
    }
    else{
        MessagePrinter::printErrorTxt("renameSet error: set "+newName+" is existing.");
        MessagePrinter::exitcfem();            
    }
    return true;        
}
bool SetManager::deleteSet(string setName, SetType setType){
    map<string,vector<PetscInt>>::iterator it;
    map<string,vector<PetscInt>> *mapPtr;
    switch (setType)
    {
    case SetType::NODE:
        mapPtr=&m_nodeSets;
        break;
    case SetType::ELEMENT:
        mapPtr=&m_elmtSets;
        break;
    case SetType::BELEMENT:
        mapPtr=&m_bElmtSets;
        break;
    default:
        MessagePrinter::printErrorTxt("deleteSet error: unsupported set type.");
        MessagePrinter::exitcfem();        
        break;
    }
    it=mapPtr->find(setName);
    if(it==mapPtr->end()){
        MessagePrinter::printErrorTxt("deleteSet error: can not delete a unexisting set");
        MessagePrinter::exitcfem();
    }
    else{
        mapPtr->erase(it);
    }
    return true;    
}

vector<PetscInt> & SetManager::getSet(string setName, SetType setType){
    map<string,vector<PetscInt>>::iterator it;
    map<string,vector<PetscInt>> *mapPtr;
    switch (setType)
    {
    case SetType::NODE:
        mapPtr=&m_nodeSets;
        break;
    case SetType::ELEMENT:
        mapPtr=&m_elmtSets;
        break;
    case SetType::BELEMENT:
        mapPtr=&m_bElmtSets;
        break;
    default:
        MessagePrinter::printErrorTxt("getSet error: unsupported set type.");
        MessagePrinter::exitcfem();        
        break;
    }
    it=mapPtr->find(setName);
    if(it==mapPtr->end()){
        MessagePrinter::printErrorTxt("getSet error: set does not exist.");
        MessagePrinter::exitcfem();
    }
    else{
        return it->second;
    }
    return mapPtr->begin()->second;
}