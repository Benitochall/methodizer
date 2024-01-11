#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "lsoda.h"
extern double exactArray[];
double exactArray[5+1];
extern double ans[];
double ans[4];
int fex(double t, double *y, double *ydot, void *data){
   ydot[0] =2 * t *  y[0]  / (1+ pow(t,2))
;         
   return(0);
}
int getExactSolution(void){
double          atol[1], rtol[1], t, tout, y[1];
int             neq = 1;
int             i;
y[0] = 1.000000;
t = 0.0E0;
tout = 1.000000000000E0;
struct lsoda_opt_t opt = {0};
opt.ixpr = 0;
opt.rtol = rtol;
opt.atol = atol;
opt.itask = 1; 
rtol[0] = 1.0E-13;
atol[0] = 1.0E-13;
struct lsoda_context_t ctx = {
.function = fex, 
.neq = neq,
.data = NULL,
.state = 1,
};
lsoda_prepare(&ctx, &opt);
exactArray[0] = y[0];
for (i= 1; i <= 5; i++) {
lsoda(&ctx, y, &t, tout);
exactArray[i] =y[0];
tout = tout + 1.0000000000;
}
lsoda_free(&ctx);
return(0);
}

double func(double t, double y){
return 2 * t * y / (1+ pow(t,2))
;
}
double calcEfficacy(double * absValArray){
double errCounter =0;
double sumError=0;  
for (int i=1; i<5 +1; i++){
if (absValArray[i] > 10e-16){
sumError += log(absValArray[i]);
errCounter++;
}
else {
continue; 
}
}
double finalError = sumError/errCounter + log(500);
return -1 *finalError; 
}

double * testForwardEuler(){
double static errorCalc[5+1];
errorCalc[0] = 0; 
double y_init =1.000000;
double total_fxn_evals =0;
double tim = (double) 0;
for (int i= 1; i <= 5; i++) {
y_init = (double)y_init + 1.000000000000*func(tim,y_init);
total_fxn_evals +=1; 
tim = tim+ 1.000000000000;
errorCalc[i] = fabs(y_init-exactArray[i]); 
}
printf("-----------------------------------------------------------------\n");
printf("The total function evaluations by Forward Euler were %0.2f \n",   total_fxn_evals);
printf("The error after 5 intervals was %0.15f \n",   y_init-exactArray[5]);
printf("-----------------------------------------------------------------\n");

return errorCalc;
}
double * testRK4(){
double static errorCalc1[5+1];
errorCalc1[0] = 0; 
double y_init = 1.000000; 
double stepsize = 1.000000000000;
double total_fxn_evals =0; 
double tim = (double) 0;
for (int i= 1; i <= 5; i++) {
tim = tim+ 1.000000000000;
double k1 =func(tim, y_init); 
double k2 = func(tim + stepsize/2.0, y_init + (0.5*stepsize*k1));
double k3 = func(tim + stepsize/2.0, y_init + (0.5 * stepsize *k2));
double k4 = func(tim + stepsize, y_init +( stepsize *k3)); 
total_fxn_evals+=4; 
y_init =y_init + (double)1/(double)6*stepsize*(k1+2*k2 +2*k3 +k4);
errorCalc1[i] = fabs(y_init-exactArray[i]); 
}
sleep(1);
printf("The total function evaluations by RK4 were %0.2f \n",   total_fxn_evals);
printf("The error after 5 intervals was %0.15f \n",   y_init-exactArray[5]);
printf("-----------------------------------------------------------------\n");

return errorCalc1;
}
double * testAdamsBashforth(){ 
double static errorCalc2[5+1];
errorCalc2[0] = 0; 
double y_init =1.000000; 
double y_0 = y_init;  
double total_fxn_evals =0;  
double y_1, y_2, y_3, y_4; 
double stepsize = 1.000000000000;
double tim = (double) 0;
for(int i =1; i<5;i++){ 
tim = tim+ 1.000000000000;
double k1 =func(tim, y_init);  
double k2 = func(tim + stepsize/2.0, y_init + (0.5*stepsize*k1)); 
double k3 = func(tim + stepsize/2.0, y_init + (0.5 * stepsize *k2)); 
double k4 = func(tim + stepsize, y_init +( stepsize *k3));  
y_init =y_init + (double)1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4);  
total_fxn_evals +=4;  
errorCalc2[i] = fabs(y_init-exactArray[i]); 
if (i==1){ 
y_1 = y_init;  
} 
if (i==2){ 
y_2 = y_init;  
} 
if (i==3){ 
y_3 = y_init; 
} 
y_4 = y_init; 
} 
double tim_0 = 0.0;
double tim_1 = 1.000000000000;
double tim_2 = 2.000000000000;
double tim_3 = 3.000000000000;
double tim_4 = 4.000000000000;
for (int i= 5; i <= 5; i++) {
y_4 = y_init;  
y_init= y_4+ 1.0 /720.0 * stepsize * ( 1901.0 * func(tim_4,y_4) -2774.0 * func(tim_3,y_3) + 2616.0 *func(tim_2, y_2) - 1274.0 *func(tim_1,y_1) +251.0 *func(tim_0,y_0)); 
y_0 = y_1;  
y_1= y_2;  
y_2= y_3;  
y_3 = y_4;  
y_4 = y_init;  
tim_0 +=1.000000000000;
tim_1 +=1.000000000000;
tim_2 +=1.000000000000;
tim_3 +=1.000000000000;
tim_4 +=1.000000000000;
total_fxn_evals +=5; 
} 
sleep(1);
printf("The total function evaluations by Adams-Bashforth were %0.2f \n",   total_fxn_evals);
printf("The error after 5 intervals was %0.15f \n",   y_init-exactArray[5]);
printf("-----------------------------------------------------------------\n");

return errorCalc2;
}	 
double * testRKF8(){ 
double y_init = 1.000000; 
double static errorCalc3[5+1];
errorCalc3[0] = 0; 
double total_fxn_evals =0;  
double stepsize = 1.000000000000;
double tim = (double) 0;
for (int i= 1; i <= 5; i++) {
tim = tim+ 1.000000000000;
double k1 =func(tim, y_init);  
double k2 = func(tim +2*stepsize/27.0, y_init + stepsize *2*(k1/27.0));  
double k3 = func(tim+3*stepsize/9.0, y_init  +  stepsize /((k1/36.0 + k2/12.0)));  
double k4 =func(tim +stepsize/6.0, y_init + stepsize * (k1/24.0 + k3/8.0));  
double k5 = func(tim + 5.0/12.0*stepsize, y_init + stepsize *(5.0/12.0*k1 + -25.0/16.0*k3+  25.0/16.0*k4)); 
double k6 = func(tim + stepsize/2.0, y_init + stepsize * (k1/20.0 + k4/4.0 + k5/5.0));  
double k7 = func(tim + 5.0/6.0*stepsize, y_init + stepsize *(-25.0/108.0*k1 + 125.0/108.0*k4 - 65.0/27.0*k5 + 125.0/54.0*k6)); 
double k8 = func(tim +stepsize/6.0, y_init + stepsize * (31.0/300.0*k1 + 61.0/225.0*k5 - 2.0/9.0*k6 + 13.0/900.0*k7)); 
double k9 = func(tim + 2.0/3.0*stepsize, y_init  + stepsize *(2.0*k1 + -53.0/6.0*k4 + 704.0/45.0*k5 + -107.0/9.0*k6 + 67.0/90.0*k7 + 3.0*k8)); 
double k10 = func(tim + stepsize/3.0, y_init + stepsize *(-91.0/108.0*k1 + 23/108.0*k4 + -976.0/135.0*k5 + 311.0/54.0*k6 + -19.0/60.0*k7 + 17.0/6.0*k8 - k9/12.0)); 
double k12 = func(tim, y_init  + stepsize * (3.0/205.0*k1 + -6.0/41.0*k6 + -3.0/205.0*k7 + -3.0/41.0*k8 + 3.0/41.0*k9 + 6.0/41.0*k10)); 
double k13 = func(tim + stepsize, y_init + stepsize *(-1777.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -289.0/82.0*k6 + 2193.0/4100.0*k7 + 51.0/82.0*k8 + 33.0/164.0*k9 + 12.0/41.0*k10 + k12));  
total_fxn_evals +=12;  
y_init =y_init + stepsize*(34.0/105.0*k6 + 9.0/35.0*k7 + 9.0/35.0*k8 + 9.0/280.0*k9 + 9.0/280.0*k10 + 41.0/840.0*k12 + 41.0/840.0*k13); 
errorCalc3[i] = fabs(y_init-exactArray[i]); 
} 
sleep(1);
printf("The total function evaluations by RK8 were %0.2f \n",   total_fxn_evals);
printf("The error after 5 intervals was %0.15f \n",   y_init-exactArray[5]);
printf("-----------------------------------------------------------------\n");

return errorCalc3;
} 

int main(){
getExactSolution();
double * euler_error;
double * rk4_error;
double * AB_error;
double * rkf8_error;
euler_error = testForwardEuler(); 
double euler_efficacy = calcEfficacy(euler_error); 
printf("The Efficacy for Forward Euler was %0.8f\n",euler_efficacy); 
printf("-----------------------------------------------------------------\n");
rk4_error = testRK4(); 
double rk4_efficacy = calcEfficacy(rk4_error); 
printf("The Efficacy for RK4 was %0.8f\n",rk4_efficacy); 
printf("-----------------------------------------------------------------\n");
 AB_error = testAdamsBashforth(); 
double AB_efficacy = calcEfficacy(AB_error); 
printf("The Efficacy for Adams-Bashforth was %0.8f\n",AB_efficacy); 
printf("-----------------------------------------------------------------\n");
rkf8_error = testRKF8(); 
double rkf8_efficacy = calcEfficacy(rkf8_error); 
printf("The Efficacy for RKF8 was %0.8f\n",rkf8_efficacy); 
printf("-----------------------------------------------------------------\n");
double efficacyArray[] = {euler_efficacy, rkf8_efficacy, AB_efficacy, rk4_efficacy}; 
double bestEfficacy = -1000; 
for (int i =0; i<4; i++){
if (efficacyArray[i] > bestEfficacy){
bestEfficacy = efficacyArray[i]; 
continue;
}
}
if (bestEfficacy == efficacyArray[0]){
printf("The method you should use is Forward Euler\n");
}
if (bestEfficacy == efficacyArray[1]){
printf("The method you should use is RKF8\n");
}
if (bestEfficacy == efficacyArray[2]){
printf("The method you should use is Adams-Bashforth\n");
}
if (bestEfficacy == efficacyArray[3]){
printf("The method you should use is RK4\n");
}
return 0;
}
