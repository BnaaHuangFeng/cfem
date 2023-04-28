#include "BCsSysStructured2d.h"
PetscErrorCode BCsSysStructured2d::init(){
    if(m_ifHasReadBCDes){
        for(vector<SingleBCDes>::iterator itBCBlock=m_bcDesPtr->begin();
        itBCBlock<m_bcDesPtr->end();++itBCBlock){// loop over every single bcs description block
            
        }
    }
    return 0;
}