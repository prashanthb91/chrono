#include "core/ChMatrix.h"
#include "core/ChTransform.h"
#include "core/ChMatrix.h"
#include "core/ChLog.h"
#include "core/ChVector.h"
#include "core/ChQuadrature.h"

using namespace chrono;

int main(int argc, char* argv[]) {


    chrono::ChMatrixDynamic<double> MatA(1000, 1000);
    
    chrono::ChMatrixDynamic<double> MatB(1000, 1000);

    
    chrono::ChMatrixDynamic<double> MatC(1000, 1000);
    chrono::ChMatrixDynamic<double> MatD(1000, 1000);

    MatA.FillRandom(-10, 10);
    MatB.FillRandom(-5,5);
    int i=0;
    int count=0;
    for(i=0; i<10; i++) {
        MatC = MatA+MatB;
        MatD = MatC-MatB;
        
        count +=MatD(0,i);
    }
    std::cout<<"Value of count is "<<count<<"\n";
}
