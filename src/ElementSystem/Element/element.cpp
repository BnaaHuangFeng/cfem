#include "ElementSystem/Element/element.h"
#include "Utils/MessagePrinter.h"
bool element::getBMatrix(void *t_dNdxPtr, int mNodePElmt, int mDofPNode, MatrixXd *BMatrixPtr){
    Vector2d *dNdxPtr=(Vector2d *)t_dNdxPtr;
    int yDofI=1,xDofI=0;
    if(mDofPNode){}
    if(!BMatrixPtr){
        MessagePrinter::printErrorTxt("B-Matrix need preallocation before calculation");
        MessagePrinter::exitcfem();
    }
    for(int nodeI=0;nodeI<mNodePElmt;nodeI++){//loop over all node 
        (*BMatrixPtr)(0,xDofI)=dNdxPtr[nodeI](0);
        (*BMatrixPtr)(0,yDofI)=0.0;
        (*BMatrixPtr)(1,xDofI)=0.0;
        (*BMatrixPtr)(1,yDofI)=dNdxPtr[nodeI](1);
        (*BMatrixPtr)(2,xDofI)=dNdxPtr[nodeI](1);
        (*BMatrixPtr)(2,yDofI)=dNdxPtr[nodeI](0);
        xDofI=yDofI+1;
        yDofI=xDofI+1;
    }
    return true;
}
bool element::getGMatrix(void *t_dNdxPtr, int mNodePElmt, int mDofPNode, MatrixXd *GMatrixPtr){
    int xDofI=0,yDofI=1;
    Vector2d *dNdxPtr=(Vector2d *)t_dNdxPtr;
    if(mDofPNode){}
    if(!GMatrixPtr){
        MessagePrinter::printErrorTxt("G Matrix need preallocation before calculation");
        MessagePrinter::exitcfem();
    }
    for(int nodeI=0;nodeI<mNodePElmt;nodeI++){
        (*GMatrixPtr)(0,xDofI)=dNdxPtr[nodeI](0);
        (*GMatrixPtr)(0,yDofI)=0.0;
        (*GMatrixPtr)(1,xDofI)=0.0;
        (*GMatrixPtr)(1,yDofI)=dNdxPtr[nodeI](0);
        (*GMatrixPtr)(2,xDofI)=dNdxPtr[nodeI](1);
        (*GMatrixPtr)(2,yDofI)=0.0;
        (*GMatrixPtr)(3,xDofI)=0.0;
        (*GMatrixPtr)(3,yDofI)=dNdxPtr[nodeI](1);   
        xDofI=yDofI+1;
        yDofI=xDofI+1;    
    }
    return true;
}