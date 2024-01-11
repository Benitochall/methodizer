/**
 *Welcome to the y'=y tester 
 * In this file We test the effecientcy of 
 * 1. Forward Euler 
 * 2. Backward Euler 
 * 3. Runge-Kutta 4
 * 
 * I first want to output how far off the soltion is at 2000 intervals 0 to 2
 *
 * Latest UPDATE: 
 * Saturday April 22 22:47
 * Need to check to see if RK4 gives the result is should, work by hand 
 *
 * Next task: Adams bashford using rk4 as the first couple steps
 *TODO check to see if rk4 is correct
 *
 *
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// the global file we will write to 
FILE *output_file; 
double exactArray[101];
double yexact(double t){

	//double ans  = pow(2.71821828459045235360287471352662497757247093699959574966967627724076630353,t); 
	double ans = exp(t); 
return ans;  

}

/*
 * For the equation y' =y the rk4 method doesn't depend on time at all
 * so we do not need to define any fxs
 */


void testForwardEuler(FILE *output_file){
int t_init = 0;
double y_init =1; // set the first value of y =0; 

//testing the time of forward Euler 
clock_t start_t, end_t; // setting clock values 
double total_t;
start_t = clock();

for(int i =1; i<=5000; i++){

	y_init = (double)y_init + 0.001*(double)y_init; 
	

	}
end_t = clock();
total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC *10E6;
printf("The total time taken by forward Euler in microseconds was: %f\n", total_t);




// file writing 
//
//printf("I think this is where the seg fault is\n"); 

printf(output_file, "%s %s %s %s %s %s\n", "The" ,"data" , "for" , "forward" , "Euler" , "is"); 
 

t_init = 0; // resets intial values 
y_init =1;
// this uses 2000 time steps between time = 0 and time = 2
fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", 0.000000, y_init, 1.00000, fabs(1.0-y_init));

for(int i =1; i<=5000; i++){

        y_init = (double)y_init + 0.001*(double)y_init;

        double y_exact = yexact((double)i/1000);

        if (i%50== 0){ // prints 100 intervals 
                double h =((double)i)/1000;
                fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", h, y_init, y_exact,fabs(y_exact-y_init));


                }
        }
}

void testBackwardEuler(FILE * output_file){

	double t_init = 0; 
	double y_init = 1; // remeber this value depends on the fucntion
	
// testing the runtime of backwards Euler	
clock_t start_t, end_t; 
double total_t; 
start_t = clock(); 

for(int i =1; i<=5000; i++){

	double stepsize = 0.001; 

	y_init = ((double)y_init)/(1-stepsize);  
		
	}
end_t = clock();
total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC *10E6;
printf("The total time taken by Backward Euler in microseconds was: %f\n", total_t);

fprintf(output_file, "The Data for backward Euler\n"); 

t_init = 0; // resets intial values 
y_init =1;
fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", 0.0, y_init, 1.000, fabs(1.0-y_init)); 

// this uses 5000 time steps between time = 0 and time = 2

for(int i =1; i<=5000; i++){

        y_init = (double)y_init + 0.001*(double)y_init;

        double y_exact = yexact((double)i/1000);

        if (i%50 == 0){ // prints 100 intervals 
                double h =((double)i)/1000;
	        fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", h, y_init, y_exact,fabs(y_exact-y_init));


                }
        }
}

void testRK4(FILE * output_file){

	double t_init = 0; 
	double y_init = 1; // remeber this value depends on the fucntion
	
// testing the runtime of RK4	
clock_t start_t, end_t; 
double total_t; 
start_t = clock(); 

double stepsize = 0.001;

for(int i =1; i<=5000; i++){

	double k1 = y_init;
	double k2 = y_init + 0.5*stepsize*k1;  
	double k3 = y_init + 0.5 * stepsize *k2; 
	double k4 = y_init + 1.0* stepsize *k3; 

	y_init =y_init + 1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4); 

	//double y_exact = yexact((double)i/1000); 
		
	}

end_t = clock();
total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC *10E6;
printf("The total time taken by RK4 in microseconds was: %f\n", total_t);

//writing to file for output data
fprintf(output_file, "The Data for RK4 Euler\n"); 

t_init = 0; // resets intial values 
y_init =1;
// this uses 2000 time steps between time = 0 and time = 2 
fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", 0.000000, y_init, 1.0000, fabs(1.0-y_init));
for(int i =1; i<=5000; i++){
	stepsize = 0.001; 
	double k1 = y_init; 
	double k2 = y_init + (0.5*stepsize*k1);
	double k3 = y_init + (0.5 * stepsize *k2);
	double k4 = y_init +( stepsize *k3); 
	y_init =y_init + (double)1/(double)6*stepsize*(k1+2*k2 +2*k3 +k4); 
	double y_exact = yexact((double)i/1000); 
	 if (i%50== 0){ // prints 100 intervals 
                double h =((double)i)/1000;
               fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", h, y_init, y_exact,fabs(y_exact-y_init));

      	 }
        }
}

/*
 *This is an implict 4 step method that uses RK4 
 * to calculate the frist 3 values 
 *
 *
 *
 */
void testAdamsBashford(FILE *output_file){
int t_init = 0; 
double y_init =1; // set the first value of y =0; 

//testing the time of forward Euler 
clock_t start_t, end_t; // setting clock values 
double total_t;

double y_0= y_init;
double y_1, y_2, y_3; 
double stepsize = 0.001; 

start_t = clock(); 

//  calculate the first 3 steps 
for(int i =1; i<=3; i++){

	double k1 = y_init;
	double k2 = y_init + 0.5*stepsize*k1;  
	double k3 = y_init + 0.5 * stepsize *k2; 
	double k4 = y_init + 1.0* stepsize *k3; 

	y_init =y_init + (double)1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4); 
	
	if (i==1){
		y_1 = y_init; }
	
	if (i==2){
		y_2 = y_init; 
	}	
	if (i==3){
		y_3 = y_init;
	}
}

for (int i=4; i<=5000; i++){

	y_init= y_2 + (double) 1 /24.0 * stepsize * (55 * y_3 + -59 * y_2 + 37 *y_1 -9 *y_0);

	// update values
		y_0 = y_1; 
		y_1= y_2; 
		y_2= y_3; 
		y_3 = y_init; 

}
	
end_t = clock();
total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC *10E6;
printf("The total time taken by Adams–Bashforth  was: %f\n", total_t); 

fprintf(output_file, "%s %s %s %s %s %s\n", "The" ,"data" , "for" , "Adams" , "Bashforth" , "is"); 

t_init = 0; // resets intial values 
y_init =1;
y_0 = y_init; 
double y_exact; 

	y_exact= yexact((double)0); 
	fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", 0.000000, y_init, y_exact, fabs(y_exact-y_init));

//  calculate the first 3 steps 
for(int i =1; i<4;i++){
	double k1 = y_init;
	double k2 = y_init + 0.5*stepsize*k1;  
	double k3 = y_init + 0.5 * stepsize *k2; 
	double k4 = y_init + 1.0* stepsize *k3; 

	y_init =y_init + (double)1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4); 
	y_exact = yexact((double)i/1000); 

	if (i==1){
		y_1 = y_init; 

	}
	
	if (i==2){
		y_2 = y_init; 
	}
	if (i==3){
		y_3 = y_init;
			
	}

}

for (int i=4; i<=5000; i++){

	//yn+4 = yn+2 + 1/24h( 55 fn+3, -59 fn+2 +37 fn+1 -9 fn)

	y_init= y_3+ (double) 1 /(double)24.0 * stepsize * (55.0 * y_3 + -59.0 * y_2 + 37.0 *y_1 -9.0 *y_0);

	// update values
		y_0 = y_1; 
		y_1= y_2; 
		y_2= y_3; 
		y_3 = y_init; 

        double y_exact = yexact((double)i/1000);

        if (i%50== 0){ // prints 100 intervals 
                double h =((double)i)/1000;
                fprintf(output_file, "%0.3f %0.15f %0.15f %0.15f\n", h, y_init, y_exact,fabs(y_exact-y_init));


                }
        }
}

void testRKF8(FILE * output_file){

	double t_init = 0; 
	double y_init = 1; // remeber this value depends on the fucntion
	
// testing the runtime of RKf8	
clock_t start_t, end_t; 
double total_t; 
start_t = clock(); 

double stepsize = 0.001;

for(int i =1; i<=5000; i++){

	double k1 = y_init;
	double k2 = y_init + stepsize *2*(k1/27.0); 
	double k3 = y_init  +  stepsize /((k1/36.0 + k2/12.0)); 

	double k4 = y_init + stepsize * (k1/24.0 + k3/8.0); 
	double k5 = y_init + stepsize *(5.0/12.0*k1 + -25.0/16.0*k3+  25.0/16.0*k4);
        double k6 = y_init  + stepsize * (k1/20.0 + k4/4.0 + k5/5.0); 
        double k7 = y_init + stepsize *(-25.0/108.0*k1 + 125.0/108.0*k4 - 65.0/27.0*k5 + 125.0/54.0*k6);

        double k8 = y_init + stepsize * (31.0/300.0*k1 + 61.0/225.0*k5 - 2.0/9.0*k6 + 13.0/900.0*k7);
        double k9 = y_init  + stepsize *(2.0*k1 + -53.0/6.0*k4 + 704.0/45.0*k5 + -107.0/9.0*k6 + 67.0/90.0*k7 + 3.0*k8);
        double k10 = y_init + stepsize *(-91/108*k1 + 23/108*k4 + -976/135*k5 + 311/54.0*k6 + -19.0/60.0*k7 + 17.0/6.0*k8 - k9/12.0);

        double k11 = y_init + stepsize*(2383.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -301.0/82.0*k6 + 2133.0/4100.0*k7 + 45.0/82.0*k8 + 45.0/164.0*k9 + 18.0/41.0*k10); 
        double k12 = y_init  + stepsize * (3.0/205.0*k1 + -6.0/41.0*k6 + -3.0/205.0*k7 + -3.0/41.0*k8 + 3.0/41.0*k9 + 6.0/41.0*k10);
        double k13 = y_init + stepsize *(-1777.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -289.0/82.0*k6 + 2193.0/4100.0*k7 + 51.0/82.0*k8 + 33.0/164.0*k9 + 12.0/41.0*k10 + k12); 

	y_init =y_init + stepsize*(34.0/105.0*k6 + 9.0/35.0*k7 + 9.0/35.0*k8 + 9.0/280.0*k9 + 9.0/280.0*k10 + 41.0/840.0*k12 + 41.0/840.0*k13);

	}

end_t = clock();
total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC *10E6;
printf("The total time taken by RKF78 in microseconds was: %f\n", total_t);

//writing to file for output data
//
fprintf(output_file, "The Data for RKF8 Euler\n"); 


t_init = 0; // resets intial values 
y_init =1;
// this uses 2000 time steps between time = 0 and time = 2 
fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", 0.000000, y_init, 1.0000, fabs(1.0-y_init));
for(int i =1; i<=5000; i++){

	stepsize = 0.001; 
	double k1 = y_init;
	double k2 = y_init +stepsize *2*k1/27.0; 
	double k3 = y_init  +(( stepsize) /(((k1/36.0) + (k2/12.0)))); 

	double k4 = y_init + stepsize * (k1/24.0 + k3/8.0); 
	double k5 = y_init + stepsize *(5.0/12.0*k1  -25.0/16.0*k3+  25.0/16.0*k4);
        double k6 = y_init  + stepsize * (k1/20.0 + k4/4.0 + k5/5.0); 
        double k7 = y_init + stepsize *(-25.0/108.0*k1 + 125.0/108.0*k4 - 65.0/27.0*k5 + 125.0/54.0*k6);

        double k8 = y_init + stepsize * (31.0/300.0*k1 + 61.0/225.0*k5 - 2.0/9.0*k6 + 13.0/900.0*k7);
        double k9 = y_init  + stepsize *(2.0*k1 + -53.0/6.0*k4 + 704.0/45.0*k5 + -107.0/9.0*k6 + 67.0/90.0*k7 + 3.0*k8);
        double k10 = y_init + stepsize *(-91.0/108.0*k1 + 23.0/108.0*k4 + -976.0/135.0*k5 + 311.0/54.0*k6 + -19.0/60.0*k7 + 17.0/6.0*k8 - k9/12.0);

        double k11 = y_init + stepsize*(2383.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -301.0/82.0*k6 + 2133.0/4100.0*k7 + 45.0/82.0*k8 + 45.0/164.0*k9 + 18.0/41.0*k10); 
        double k12 = y_init  + stepsize * (3.0/205.0*k1 + -6.0/41.0*k6 + -3.0/205.0*k7 + -3.0/41.0*k8 + 3.0/41.0*k9 + 6.0/41.0*k10);
        double k13 = y_init + stepsize *(-1777.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -289.0/82.0*k6 + 2193.0/4100.0*k7 + 51.0/82.0*k8 + 33.0/164.0*k9 + 12.0/41.0*k10 + k12); 

	y_init =y_init + stepsize*(34.0/105.0*k6 + 9.0/35.0*k7 + 9.0/35.0*k8 + 9.0/280.0*k9 + 9.0/280.0*k10 + 41.0/840.0*k12 + 41.0/840.0*k13);



        
	 if (i%50== 0){ // prints 100 intervals 
                double time =((double)i)/1000;
		double y_exact = yexact(time); 
		
               fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", time, y_init, y_exact,fabs(y_exact-y_init));


      	 }
        }
}

int main(){

FILE *output_file = fopen("y_y_out.txt", "w"); // write only
// test if file is NULL

if (output_file == NULL){
	printf("Could not open file"); 
	exit(-1); 
} 

testForwardEuler(output_file); 

testBackwardEuler(output_file);

testRK4(output_file); 

testAdamsBashford(output_file); 

testRKF8(output_file); 

fclose(output_file); 
return 0; 
}
