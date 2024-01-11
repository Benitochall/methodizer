#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"
#include "lsoda.h"

/*
 * This file runs varius iterative methods on
 * a complicated polynomial
 *
 *
 *
 *
 */

// the global file we will write to 
FILE *output_file; 

// the global array of exact soultions
extern double exactArray[];
double exactArray[101]; 
/*Sets up funtion to integrate
*/
int fex(double t, double *y, double *ydot, void *data)
{
        ydot[0] = -0.3*pow(t,4)+pow(t,3)-2 +y[0]*pow(t,2)+ -0.586*t*y[0]; 
        return(0);
}

//Method to calculate the exact solutions
double poly(double t, double y){

        return -0.3*pow(t,4)+pow(t,3)-2 +y*pow(t,2)+ -0.586*t*y; 
}

int getExactSolution(void){
        // atol rtol need to be size of how many dimentions 
        // t = starting value of time 
        // tout = the first place you want output 
        //
        double          atol[1], rtol[1], t, tout, y[1];
        int             neq = 1;
        int             i;

        y[0] = 0.02; // the starting values of each paramter
                      // ie the initial coditions
        t = 0.0E0; // the starting value of t
        tout = 0.001E0; // the first point where we want an answer


        struct lsoda_opt_t opt = {0}; // a strct that i dont know the purpose of yet
        opt.ixpr = 0;//setting ixpr value 
        opt.rtol = rtol;// setting rtol value
        opt.atol = atol;
        opt.itask = 1; // set to 1 for normal compuations 

        rtol[0] = 1.0E-10; // realative tolerance param
        atol[0] = 1.0E-10;// absolute tollerance param


        struct lsoda_context_t ctx = {
                .function = fex, // the function to be passed in
                .neq = neq, // number of first orcer eqns
                .data = NULL,
                .state = 1,
        };

        lsoda_prepare(&ctx, &opt);// binds together ctx and opt

	// build the array of exact solutions
	exactArray[0] = y[0];
      	int j=1; 	
        for (i= 1; i <= 5000; i++) {

                lsoda(&ctx, y, &t, tout);
		 
                if (i%50 == 0){
	
                        //double exact = abs(y[0] - exactexp(t));
			exactArray[j] =y[0];
                        
			j= j+1; 
		       // sets the exact value 	
                }

                tout = tout + 0.001; // increments time step by 0.001
        }
	
        lsoda_free(&ctx); // frees the struct 
        return(0);
}

void testForwardEuler(FILE *output_file){
double t_init = 0;
double y_init =0.02; 
double total_fxn_evals =0; 

for(int i =1; i<=5000; i++){

        double tim = (double) i/1000;

        y_init = (double)y_init + 0.001*poly(tim,y_init); 
	total_fxn_evals+=1; 

        }
printf("The total fucntion evaluations by Forward Euler was: %f\n", total_fxn_evals);




// file writing  

fprintf(output_file, "%s %s %s %s %s %s\n", "The" ,"data" , "for" , "forward" , "Euler" , "is");
 

t_init = 0; // resets intial values 
y_init =0.02;
// this uses 2000 time steps between time = 0 and time = 2
fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", t_init, y_init, exactArray[0], (y_init-exactArray[0])); 
int j =1; 
for(int i =1; i<=5000; i++){

        double tim = (double) i/1000;

        y_init = (double)y_init + 0.001*poly(tim,y_init); 

        if (i%50== 0){ // prints 100 intervals 
                double h =((double)i)/1000;
		double y_exact = exactArray[j];
		j = j+1; 
                fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", h, y_init,y_exact, y_init-y_exact); 


                }
        }
}
void testRK4(FILE * output_file){
 
        double y_init = 0.02; // remeber this value depends on the fucntion
        
// testing fxn evals
double stepsize = 0.001;
double fxn_evals =0; 
for(int i =1; i<=5000; i++){
	double tim = (double) i/1000;

        
	fxn_evals +=4; 
        double k1 = poly(tim,y_init);
        double k2 = poly(tim+ stepsize/2,y_init + 0.5*stepsize*k1);  
        double k3 = poly(tim+ stepsize/2,y_init + 0.5 * stepsize *k2); 
        double k4 = poly(tim,y_init + 1.0* stepsize *k3); 

        y_init =y_init + 1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4); 

        }

printf("The total fxn evals for RK4 was: %f\n", fxn_evals);

//writing to file for output data
fprintf(output_file, "The Data for RK4 Euler\n"); 

double t_init = 0; // resets intial values 
y_init =0.02;
// this uses 2000 time steps between time = 0 and time = 2 
fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", t_init, y_init, exactArray[0], (y_init-exactArray[0]));
int j=1; 
for(int i =1; i<=5000; i++){
        stepsize = 0.001; 
        double tim = (double) i/1000;

        
	fxn_evals +=4; 
        double k1 = poly(tim,y_init);
        double k2 = poly(tim+ stepsize/2,y_init + 0.5*stepsize*k1);  
        double k3 = poly(tim+ stepsize/2,y_init + 0.5 * stepsize *k2); 
        double k4 = poly(tim,y_init + 1.0* stepsize *k3); 

        y_init =y_init + 1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4); 
        if (i%50== 0){ // prints 100 intervals 
                double h =((double)i)/1000;
                double y_exact = exactArray[j];
                j = j+1; 
                fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", h, y_init,y_exact, y_init-y_exact); 
         }
        }
}

void testAdamsBashford(FILE *output_file){
	double t_init = 0; 
	double y_init =0.02; // set the first value of y =0; 

//testing the fxn evaluations 
double total_fxn_evals=0; 
double y_0= y_init;
double y_1, y_2, y_3; 
double stepsize = 0.001; 
 

//  calculate the first 3 steps 
for(int i =1; i<=3; i++){

        double tim = (double) i/1000;

        double k1 = poly(tim,y_init);
        double k2 = poly(tim+ stepsize/2,y_init + 0.5*stepsize*k1);  
        double k3 = poly(tim+ stepsize/2,y_init + 0.5 * stepsize *k2); 
        double k4 = poly(tim,y_init + 1.0* stepsize *k3); 

        y_init =y_init + 1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4); 
	total_fxn_evals +=4; 
        
        if (i==1){
                y_1 = y_init; }
        
        if (i==2){
                y_2 = y_init; 
        }       
        if (i==3){
                y_3 = y_init;
        }
}
// Running Adams-Bashforth for the rest of the intervals
for (int i=4; i<=5000; i++){
        y_3 = y_init; 
        double tim = (double) i/1000;

        y_init= y_3 + (double) 1 /24.0 * stepsize * (55 *poly(tim,y_3) + -59 *poly(tim, y_2) + 37 *poly(tim,y_1) -9 *poly(tim, y_0));
	total_fxn_evals +=4; 
        // update values
                y_0 = y_1; 
                y_1= y_2; 
                y_2= y_3; 
                y_3 = y_init; 

}

printf("The total Function Evaluations of Adams–Bashforth  was: %f\n", total_fxn_evals); 

fprintf(output_file, "%s %s %s %s %s %s\n", "The" ,"data" , "for" , "Adams" , "Bashforth" , "is"); 
t_init = 0; // resets intial values 
y_init =0.02;
y_0 = y_init; 
 
        fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", t_init, y_init, exactArray[0], (y_init-exactArray[0]));

//  calculate the first 3 steps 
for(int i =1; i<4;i++){
        double tim = (double) i/1000;

        double k1 = poly(tim,y_init);
        double k2 = poly(tim+ stepsize/2,y_init + 0.5*stepsize*k1);  
        double k3 = poly(tim+ stepsize/2,y_init + 0.5 * stepsize *k2); 
        double k4 = poly(tim,y_init + 1.0* stepsize *k3); 

        y_init =y_init + 1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4); 
	total_fxn_evals +=4; 
	
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
int j=1; 

for (int i=4; i<=5000; i++){

        y_3 = y_init; 
        double tim = (double) i/1000;

        y_init= y_3 + (double) 1 /24.0 * stepsize * (55 *poly(tim,y_3) + -59 *poly(tim, y_2) + 37 *poly(tim,y_1) -9 *poly(tim, y_0));
	total_fxn_evals +=4; 
        // update values
                y_0 = y_1; 
                y_1= y_2; 
                y_2= y_3; 
                y_3 = y_init; 
        if (i%50== 0){ // prints 100 intervals 
               double y_exact = exactArray[j];
		double h =((double)i)/1000;
		j=j+1; 
		fprintf(output_file, "%0.3f %0.15f %0.15f %0.15f\n", h, y_init, y_exact, y_init -y_exact);
		}
	}	
}
void testRKF8(FILE * output_file){
// Determining the amount of function evaluations
double y_init = 0.02; // remeber this value depends on the fucntion
double total_fxn_evals =0; 

double stepsize = 0.001;

for(int i =1; i<=5000; i++){
        //RUnning RFK8
        double tim = (double) i/1000;



        double k1 =poly(tim, y_init);
        double k2 = poly(tim,y_init + stepsize *2*(k1/27.0)); 
        double k3 = poly(tim, y_init  +  stepsize /((k1/36.0 + k2/12.0))); 

        double k4 =poly(tim, y_init + stepsize * (k1/24.0 + k3/8.0)); 
        double k5 = poly(tim,y_init + stepsize *(5.0/12.0*k1 + -25.0/16.0*k3+  25.0/16.0*k4));
        double k6 = poly(tim,y_init  + stepsize * (k1/20.0 + k4/4.0 + k5/5.0)); 
        double k7 = poly(tim,y_init + stepsize *(-25.0/108.0*k1 + 125.0/108.0*k4 - 65.0/27.0*k5 + 125.0/54.0*k6));

        double k8 = poly(tim,y_init + stepsize * (31.0/300.0*k1 + 61.0/225.0*k5 - 2.0/9.0*k6 + 13.0/900.0*k7));
        double k9 = poly(tim,y_init  + stepsize *(2.0*k1 + -53.0/6.0*k4 + 704.0/45.0*k5 + -107.0/9.0*k6 + 67.0/90.0*k7 + 3.0*k8));
        double k10 = poly(tim, y_init + stepsize *(-91/108*k1 + 23/108*k4 + -976/135*k5 + 311/54.0*k6 + -19.0/60.0*k7 + 17.0/6.0*k8 - k9/12.0));

 //       double k11 = y_init + stepsize*(2383.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -301.0/82.0*k6 + 2133.0/4100.0*k7 + 45.0/82.0*k8 + 45.0/164.0*k9 + 18.0/41.0*k10); 
        double k12 = poly(tim,y_init  + stepsize * (3.0/205.0*k1 + -6.0/41.0*k6 + -3.0/205.0*k7 + -3.0/41.0*k8 + 3.0/41.0*k9 + 6.0/41.0*k10));
        double k13 = poly(tim,y_init + stepsize *(-1777.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -289.0/82.0*k6 + 2193.0/4100.0*k7 + 51.0/82.0*k8 + 33.0/164.0*k9 + 12.0/41.0*k10 + k12)); 

	total_fxn_evals +=12; 

        y_init =y_init + stepsize*(34.0/105.0*k6 + 9.0/35.0*k7 + 9.0/35.0*k8 + 9.0/280.0*k9 + 9.0/280.0*k10 + 41.0/840.0*k12 + 41.0/840.0*k13);

        }

printf("The total number of function Evaluations was %f\n", total_fxn_evals);

//writing to file for output data
//
fprintf(output_file, "The Data for RKF8 Euler\n"); 


double t_init = 0; // resets intial values 
y_init =0.02;
// this uses 2000 time steps between time = 0 and time = 2 

fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", t_init, y_init, exactArray[0], y_init-exactArray[0]);
int j =1; 
for(int i =1; i<=5000; i++){

        double tim = (double) i/1000;

        stepsize = 0.001;        
       	double k1 =poly(tim,y_init);
        double k2 = poly(tim,y_init + stepsize *2*(k1/27.0)); 
        double k3 = poly(tim, y_init  +  stepsize /((k1/36.0 + k2/12.0))); 

        double k4 =poly(tim, y_init + stepsize * (k1/24.0 + k3/8.0)); 
        double k5 = poly(tim,y_init + stepsize *(5.0/12.0*k1 + -25.0/16.0*k3+  25.0/16.0*k4));
        double k6 = poly(tim,y_init  + stepsize * (k1/20.0 + k4/4.0 + k5/5.0)); 
        double k7 = poly(tim,y_init + stepsize *(-25.0/108.0*k1 + 125.0/108.0*k4 - 65.0/27.0*k5 + 125.0/54.0*k6));

        double k8 = poly(tim,y_init + stepsize * (31.0/300.0*k1 + 61.0/225.0*k5 - 2.0/9.0*k6 + 13.0/900.0*k7));
        double k9 = poly(tim,y_init  + stepsize *(2.0*k1 + -53.0/6.0*k4 + 704.0/45.0*k5 + -107.0/9.0*k6 + 67.0/90.0*k7 + 3.0*k8));
        double k10 = poly(tim, y_init + stepsize *(-91.0/108.0*k1 + 23.0/108.0*k4 + -976.0/135.0*k5 + 311.0/54.0*k6 + -19.0/60.0*k7 + 17.0/6.0*k8 - k9/12.0));

 //     double k11 = y_init + stepsize*(2383.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -301.0/82.0*k6 + 2133.0/4100.0*k7 + 45.0/82.0*k8 + 45.0/164.0*k9 + 18.0/41.0*k10); 
        double k12 = poly(tim,y_init  + stepsize * (3.0/205.0*k1 + -6.0/41.0*k6 + -3.0/205.0*k7 + -3.0/41.0*k8 + 3.0/41.0*k9 + 6.0/41.0*k10));
        double k13 = poly(tim,y_init + stepsize *(-1777.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -289.0/82.0*k6 + 2193.0/4100.0*k7 + 51.0/82.0*k8 + 33.0/164.0*k9 + 12.0/41.0*k10 + k12));
        
	y_init =y_init + stepsize*(34.0/105.0*k6 + 9.0/35.0*k7 + 9.0/35.0*k8 + 9.0/280.0*k9 + 9.0/280.0*k10 + 41.0/840.0*k12 + 41.0/840.0*k13);



        
         if (i%50== 0){ // prints 100 intervals 
                double time =((double)i)/1000;
                double y_exact = exactArray[j];
		j= j+1;	
               	fprintf(output_file, "%0.2f %0.15f %0.15f %0.15f\n", time, y_init, y_exact, y_init-y_exact);


         }
        }
}

int main(){

FILE *output_file = fopen("ysin1y_out.txt", "w"); // write only
// test if file is NULL

if (output_file == NULL){
	printf("Could not open file"); 
	exit(-1); 
} 
// creates and fills an array of exact solutions
getExactSolution();

// Test all Iterative methods used and prints their results
testForwardEuler(output_file); 

testRK4(output_file); 

testAdamsBashford(output_file); 

testRKF8(output_file); 

fclose(output_file); 
return 0; 
}
