#include "core/ChMatrix.h"
#include "core/ChTransform.h"
#include "core/ChMatrix.h"
#include "core/ChLog.h"
#include "core/ChVector.h"
#include "core/ChQuadrature.h"
#include <time.h>

#define SIZE 1000
#define TIME 1
#define T_START 1
#define T_STOP 0
#define COMMENT 1 
#define PRINT if(COMMENT)


using namespace chrono;
timespec m_start,m_stop;

typedef double REAL;

void measure_time(bool flag) {
        if(!TIME) {
		return;
	}
	if(flag) {
	    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &m_start);
	}
	else {
	    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &m_stop);
	    float time_diff;
            time_diff = (m_stop.tv_sec - m_start.tv_sec)*1000000000 + (m_stop.tv_nsec - m_start.tv_nsec);
	    time_diff /= 1000000;
	    std::cout<<"TIME TAKEN : "<<time_diff<<" ms\n";
	} 
	    

}
void test_matrix_scale(){
	PRINT std::cout<<"----Executing matrix scale operation----\n";
     	chrono::ChMatrixDynamic<REAL> MatA(SIZE,SIZE);
     	chrono::ChMatrixDynamic<REAL> Scale(SIZE,SIZE);

     	MatA.FillRandom(-10,10);
     	Scale.FillRandom(-1,1);
     	measure_time(T_START);
     	MatA.MatrScale(Scale);
	measure_time(T_STOP);     
}

void test_matrix_increment(){
   	PRINT std::cout<<" -------- Executing Matrix Increment --------- \n";
    	chrono::ChMatrixDynamic<REAL> MatA(SIZE,SIZE); 
    	chrono::ChMatrixDynamic<REAL> Incr(SIZE,SIZE);
	measure_time(T_START);
    	MatA.MatrInc(Incr);
	measure_time(T_STOP);
}
void test_matrix_decrement(){
    PRINT std::cout<< " --------- Executing Matrix Decrement --------- \n";

    chrono::ChMatrixDynamic<REAL> MatA(SIZE,SIZE); 
    chrono::ChMatrixDynamic<REAL> Incr(SIZE,SIZE);
	measure_time(T_START);
    MatA.MatrDec(Incr);
	measure_time(T_STOP);
}
void test_matrix_dotproduct() {

	PRINT std::cout<<"--------- Executing Matrix DotProduct ---------\n";

	ChMatrixDynamic<REAL> MatA(SIZE,SIZE);
	ChMatrixDynamic<REAL> MatB(SIZE,SIZE);
	
	MatA.FillRandom(-10, 10);
	MatB.FillRandom(-5, 5);

	REAL dot_product;

	measure_time(T_START);
	dot_product = ChMatrix<>::MatrDot(&MatA, &MatB);
	measure_time(T_STOP);

        //PRINT	std::cout << "Dot product is :"<<dot_product<< "\n";


}

void test_matrix_sub() {
	PRINT std::cout<< " ------ Executing Matrix Subtraction ------\n";
	
	chrono::ChMatrixDynamic<REAL> MatA(SIZE,SIZE);
	chrono::ChMatrixDynamic<REAL> MatB(SIZE,SIZE);

	MatA.FillRandom(-10,10);
	MatB.FillRandom(-8,8);

	chrono::ChMatrixDynamic<REAL> MatC(SIZE, SIZE);

	measure_time(T_START);
	MatC.MatrSub(MatA, MatB);
	measure_time(T_STOP);
}   
                                 
void test_matrix_add() {
    PRINT std::cout<< " ---- Executing Matrix Addition ----\n";
     
    chrono::ChMatrixDynamic<REAL> MatA(SIZE,SIZE);
    
    chrono::ChMatrixDynamic<REAL> MatB(SIZE,SIZE);

    
    chrono::ChMatrixDynamic<REAL> MatC(SIZE,SIZE);
    chrono::ChMatrixDynamic<REAL> MatD(SIZE,SIZE);

    MatA.FillRandom(-10,10);
    MatB.FillRandom(-5,5);

    int i=0;
    int count=0;
    measure_time(T_START);
    MatC.MatrAdd(MatA, MatB);
    measure_time(T_STOP);

}

int main(int argc, char* argv[]) {


    chrono::ChMatrixDynamic<REAL> MatA(SIZE,SIZE);
    
    chrono::ChMatrixDynamic<REAL> MatB(SIZE,SIZE);

    
    chrono::ChMatrixDynamic<REAL> MatC(SIZE,SIZE);
    chrono::ChMatrixDynamic<REAL> MatD(SIZE,SIZE);

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

    case 4: {
		test_matrix_dotproduct();
		break;
	    }
    case 5: {
		test_matrix_increment();
		break;
	    }
    case 99:{
		test_matrix_add();
		test_matrix_sub();
		test_matrix_scale();
		test_matrix_dotproduct();
		test_matrix_increment();
        test_matrix_decrement();
	        break;
	    }
    
    }

}
