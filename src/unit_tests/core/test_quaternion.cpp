#include "core/ChMatrix.h"
#include "core/ChTransform.h"
#include "core/ChQuaternion.h"
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
typedef float Real;


using namespace chrono;
timespec m_start,m_stop;

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
void test_quaternion_scale(){
    PRINT std::cout<<"----Executing quaternion scale operation----\n";
    Real a = 3542.542352;
    Real scale1 = 3.434;
    Real scale2 = 2.235;
    Real b = (Real)rand();
    srand((unsigned)time(NULL));
    Real c = (Real)rand();
    Real d = (Real)rand();
    chrono::ChQuaternion<Real> A(a,b,c,d);
    chrono::ChQuaternion<Real> B(b,a,d,c);
    chrono::ChQuaternion<Real> C;
     
    int i=0;
    int count=0;
    measure_time(T_START);
    A.Scale(scale1); 
    for(i=0; i < 1000; i++)
    {
      if(i%3 ==0)
         B.Scale(scale2);
      else
         A.Scale(1);
    }
    B.Scale(scale2);
    measure_time(T_STOP);
}



void test_quaternion_sub() {
    PRINT std::cout<< " ---- Executing Quaternion Subtraction ----\n";
    Real a = 3542.542352;
    Real b = (Real)rand();
    srand((unsigned)time(NULL));
    Real c = (Real)rand();
    Real d = (Real)rand();
    chrono::ChQuaternion<Real> A(a,b,c,d);
    chrono::ChQuaternion<Real> B(b,a,d,c);
    chrono::ChQuaternion<Real> C;
     
    int i=0;
    int count=0;
    measure_time(T_START);
    C.Sub(A,B); 
    for(i=0; i < 1000; i++)
    {
      if(i%3 ==0)
         A.Sub(C,B);
      else
         B.Sub(C,A);
    }
    measure_time(T_STOP);
}   
                                 
void test_quaternion_add() {
    PRINT std::cout<< " ---- Executing Quaternion Addition ----\n";
    Real a = 3542.542352;
    Real b = (Real)rand();
    srand((unsigned)time(NULL));
    Real c = (Real)rand();
    Real d = (Real)rand();
    chrono::ChQuaternion<Real> A(a,b,c,d);
    chrono::ChQuaternion<Real> B(b,a,d,c);
    chrono::ChQuaternion<Real> C;
     
    int i=0;
    int count=0;
    measure_time(T_START);
    C.Add(A,B); 
    for(i=0; i < 1000; i++)
    {
      if(i%3 ==0)
         A.Add(C,B);
      else
         B.Add(C,A);
    }
    C.Add(A,B); 
    measure_time(T_STOP);

}

int main(int argc, char* argv[]) {
    
    Real a = 3542.542352;
    Real b = (Real)rand();
    srand((unsigned)time(NULL));
    Real c = (Real)rand();
    Real d = (Real)rand();
    chrono::ChQuaternion<Real> A(a,b,c,d);
    chrono::ChQuaternion<Real> B(b,a,d,c);

    int i=0;
    int count=0;
    if(argc == 1) {
        std::cout<<"Need to know which test cases to run.";
        return 1;
    }

    
    int test_case = atoi(argv[1]);
    switch(test_case) {


    case 1: {
                test_quaternion_add();
                break;
            }
    case 2: {
                test_quaternion_sub();
                break;
            }
    case 3: {
                test_quaternion_scale();
                break;
            }

    case 99:{
		test_quaternion_add();
		test_quaternion_sub();
		test_quaternion_scale();
	        break;
	    }
    
    }

}
