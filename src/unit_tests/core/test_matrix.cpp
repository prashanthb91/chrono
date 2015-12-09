#include "core/ChMatrix.h"
#include "core/ChTransform.h"
#include "core/ChMatrix.h"
#include "core/ChLog.h"
#include "core/ChVector.h"
#include "core/ChQuadrature.h"
#include <time.h>

#define SIZE 1000
using namespace chrono;

void test_matrix_scale(){
     chrono::ChMatrixDynamic<double> MatA(SIZE,SIZE);
     chrono::ChMatrixDynamic<double> Scale(SIZE,SIZE);

     MatA.FillRandom(-10,10);
     Scale.FillRandom(-1,1);
     
     MatA.MatrScale(Scale);
}

void test_matrix_increment(){
    chrono::ChMatrixDynamic<double> double_MatA(SIZE,SIZE); 
    chrono::ChMatrixDynamic<double> double_Incr(SIZE,SIZE);
}
void test_matrix_sub() {
     std::cout<< "Compiler, you are a big pain in the ass";
}                                    
void test_matrix_add() {
     
    chrono::ChMatrixDynamic<double> MatA(SIZE,SIZE);
    
    chrono::ChMatrixDynamic<double> MatB(SIZE,SIZE);

    
    chrono::ChMatrixDynamic<double> MatC(SIZE,SIZE);
    chrono::ChMatrixDynamic<double> MatD(SIZE,SIZE);

    MatA.FillRandom(-10,10);
    MatB.FillRandom(-5,5);

    int i=0;
    int count=0;
    for(i=0;i<10;i++) {
        MatC = MatA + MatB;
        MatD.MatrAdd(MatA,MatB); 
        count += MatD(0,i);;
    }

}

int main(int argc, char* argv[]) {


    chrono::ChMatrixDynamic<double> MatA(SIZE,SIZE);
    
    chrono::ChMatrixDynamic<double> MatB(SIZE,SIZE);

    
    chrono::ChMatrixDynamic<double> MatC(SIZE,SIZE);
    chrono::ChMatrixDynamic<double> MatD(SIZE,SIZE);

    MatA.FillRandom(-10, 10);
    MatB.FillRandom(-5,5);
    int i=0;
    int count=0;
    if(argc == 1) {
        std::cout<<"Need to know which test cases to run.";
        return 1;
    }

    
    int test_case = atoi(argv[1]);
    switch(test_case) {


    case 1: {
                test_matrix_add();
                break;
            }
    case 2: {
                test_matrix_sub();
                break;
            }
    case 3: {
                test_matrix_scale();
                break;
            }
    
    }

}
