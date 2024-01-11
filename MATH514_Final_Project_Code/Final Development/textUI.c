#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

/*
 *
 * Welcome to the source file behind the methodizer
 *
 * This file creates a unique ode.c file that
 * is then complied and run with the users specified
 * inputs.
 *
 /
 
/*
*This method replaces any instance of y with y[0]
* so that the equation can be passed into lsoda
*/
void replace_y_with_y0(char *str) {
    int len = strlen(str);
    char newstr[len*2+1]; 
    int newlen = 0; 
    for (int i = 0; i < len; i++) {
        if (str[i] == 'y') {
            
            newstr[newlen++] = ' ';
            newstr[newlen++] = 'y';
            newstr[newlen++] = '[';
            newstr[newlen++] = '0';
            newstr[newlen++] = ']';
            newstr[newlen++] = ' ';
        } else {
            // Copy current character to new string
            newstr[newlen++] = str[i];
        }
    }

    // Add null terminator to new string
    newstr[newlen] = '\0';

    // Replace original string with new string
    strcpy(str, newstr);
}



/*
*This method builds a unique file for every differnt funtion passed in 
*fxn: the function to integrate 
* orig: the original function before y was replaced with y[0]
* h: the specified stepsize
* ystart: the intial value of y
*/
int createFile(char *fxn, char * orig, char* h, double ystart){ 

    // Creates a file to run
    char *filename = "ode.c";
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Failed to open file %s\n", filename);
        return 1;
    }

    //Calculates all constants needed to pass in to the method 
    double int_val = 5./strtod(h,NULL); 
    double step_size = strtod(h,NULL);

    int num_intervals = round(int_val); 
    int num_intervals_10 = num_intervals/10.0;

    
    // writes the unique file later to be compiled
    fprintf(fp, "#include <math.h>\n");
    fprintf(fp, "#include <unistd.h>\n");
    fprintf(fp, "#include <stdlib.h>\n");
    fprintf(fp, "#include <stdio.h>\n"); 
    fprintf(fp, "#include \"common.h\"\n");
    fprintf(fp, "#include \"lsoda.h\"\n");


    // writing global variables 
    fprintf(fp,"extern double exactArray[];\n");
    fprintf(fp, "double exactArray[%d+1];\n", num_intervals); //TODO need to change this value 
    fprintf(fp,"extern double ans[];\n");
    fprintf(fp, "double ans[4];\n"); 
    
    // the fex function requried for lsoda
    fprintf(fp, "int fex(double t, double *y, double *ydot, void *data){\n");
    fprintf(fp, "   ydot[0] =");
    fprintf(fp, "%s",fxn);
    fprintf(fp, ";         \n");
    fprintf(fp, "   return(0);\n");
    fprintf(fp, "}\n"); 

    // using lsoda to get a value of the exact solution
    fprintf(fp, "int getExactSolution(void){\n"); 
    fprintf(fp, "double          atol[1], rtol[1], t, tout, y[1];\n");
    fprintf(fp, "int             neq = 1;\n"); 
    fprintf(fp, "int             i;\n"); 

        fprintf(fp, "y[0] = %f;\n", ystart); 
        fprintf(fp, "t = 0.0E0;\n");
        fprintf(fp, "tout = %0.12fE0;\n", step_size);
        fprintf(fp, "struct lsoda_opt_t opt = {0};\n");
        fprintf(fp, "opt.ixpr = 0;\n");
        fprintf(fp, "opt.rtol = rtol;\n");
        fprintf(fp, "opt.atol = atol;\n");
        fprintf(fp, "opt.itask = 1; \n");

        fprintf(fp, "rtol[0] = 1.0E-13;\n");
        fprintf(fp, "atol[0] = 1.0E-13;\n");


        fprintf(fp, "struct lsoda_context_t ctx = {\n");
                fprintf(fp, ".function = fex, \n");
                fprintf(fp, ".neq = neq,\n");
                fprintf(fp, ".data = NULL,\n");
                fprintf(fp, ".state = 1,\n");
        fprintf(fp, "};\n");

        fprintf(fp, "lsoda_prepare(&ctx, &opt);\n");

	
	fprintf(fp, "exactArray[0] = y[0];\n");
        // fills in exact array with lsoda values 
        fprintf(fp, "for (i= 1; i <= %d; i++) {\n", num_intervals);
                fprintf(fp, "lsoda(&ctx, y, &t, tout);\n");
			fprintf(fp, "exactArray[i] =y[0];\n");


         

                fprintf(fp, "tout = tout + %0.10f;\n", strtod(h,NULL)); // increments time step by 0.001 // TODO check this
        fprintf(fp, "}\n");
	
        fprintf(fp, "lsoda_free(&ctx);\n"); // frees the struct 
        
        fprintf(fp, "return(0);\n");
        fprintf(fp, "}\n");
        
        fprintf(fp,"\n");

        //a unique function for each ode is created
        fprintf(fp,"double func(double t, double y){\n");

        fprintf(fp,"return %s;\n", orig);
        fprintf(fp,"}\n");

    
    // method for calculating efficacy
    fprintf(fp,"double calcEfficacy(double * absValArray){\n");
    fprintf(fp,"double errCounter =0;\n");
    fprintf(fp,"double sumError=0;  \n");
    fprintf(fp,"for (int i=1; i<%d +1; i++){\n", num_intervals);
    fprintf(fp,"if (absValArray[i] > 10e-16){\n");
    fprintf(fp,"sumError += log(absValArray[i]);\n");
    fprintf(fp,"errCounter++;\n");
    fprintf(fp,"}\n");
    fprintf(fp,"else {\n");
    fprintf(fp,"continue; \n");
    fprintf(fp,"}\n");
    fprintf(fp,"}\n");
    fprintf(fp,"double finalError = sumError/errCounter + log(500);\n");
    fprintf(fp,"return -1 *finalError; \n");
    fprintf(fp,"}\n");
    fprintf(fp,"\n");


        // testing the forward euler method 
        fprintf(fp, "double * testForwardEuler(){\n");
        fprintf(fp, "double static errorCalc[%d+1];\n", num_intervals);

        fprintf(fp, "errorCalc[0] = 0; \n"); 

        fprintf(fp, "double y_init =%f;\n", ystart);
        fprintf(fp, "double total_fxn_evals =0;\n");
        fprintf(fp, "double tim = (double) 0;\n"); 


        fprintf(fp, "for (int i= 1; i <= %d; i++) {\n", num_intervals);
        
        

        fprintf(fp, "y_init = (double)y_init + %0.12f*func(tim,y_init);\n", strtod(h,NULL));
        fprintf(fp, "total_fxn_evals +=1; \n");
        fprintf(fp, "tim = tim+ %0.12f;\n", step_size);

        // need to update errorCalc

        fprintf(fp, "errorCalc[i] = fabs(y_init-exactArray[i]); \n"); 



        fprintf(fp, "}\n");
    fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");
    fprintf(fp, "printf(\"The total function evaluations by Forward Euler were %%0.2f \\n\",   total_fxn_evals);\n");

    fprintf(fp, "printf(\"The error after %d intervals was %%0.15f \\n\",   y_init-exactArray[%d]);\n", num_intervals, num_intervals);
    fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");
    fprintf(fp,"\n");
    fprintf(fp, "return errorCalc;\n");

    fprintf(fp, "}\n");

    // Construting the rk4 subroutine 
    fprintf(fp, "double * testRK4(){\n");

    fprintf(fp, "double static errorCalc1[%d+1];\n", num_intervals);

    fprintf(fp, "errorCalc1[0] = 0; \n"); 

 
    fprintf(fp, "double y_init = %f; \n", ystart);

    fprintf(fp, "double stepsize = %0.12f;\n", step_size);
    fprintf(fp, "double total_fxn_evals =0; \n");
    fprintf(fp, "double tim = (double) 0;\n"); 
    fprintf(fp, "for (int i= 1; i <= %d; i++) {\n", num_intervals);

    fprintf(fp, "tim = tim+ %0.12f;\n", step_size);
        
    fprintf(fp, "double k1 =func(tim, y_init); \n");
    fprintf(fp, "double k2 = func(tim + stepsize/2.0, y_init + (0.5*stepsize*k1));\n");
    fprintf(fp, "double k3 = func(tim + stepsize/2.0, y_init + (0.5 * stepsize *k2));\n");
    fprintf(fp, "double k4 = func(tim + stepsize, y_init +( stepsize *k3)); \n");

        fprintf(fp, "total_fxn_evals+=4; \n");

        fprintf(fp, "y_init =y_init + (double)1/(double)6*stepsize*(k1+2*k2 +2*k3 +k4);\n");
        fprintf(fp, "errorCalc1[i] = fabs(y_init-exactArray[i]); \n"); 
        fprintf(fp, "}\n");
        fprintf(fp, "sleep(1);\n");
        fprintf(fp, "printf(\"The total function evaluations by RK4 were %%0.2f \\n\",   total_fxn_evals);\n");

        
        fprintf(fp, "printf(\"The error after %d intervals was %%0.15f \\n\",   y_init-exactArray[%d]);\n", num_intervals, num_intervals);
        
        fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");
        fprintf(fp,"\n");
       fprintf(fp, "return errorCalc1;\n");

        fprintf(fp, "}\n");

    //construct adams bashforth subroutine
    fprintf(fp, "double * testAdamsBashforth(){ \n");
     fprintf(fp, "double static errorCalc2[%d+1];\n", num_intervals);

    fprintf(fp, "errorCalc2[0] = 0; \n"); 


    fprintf(fp, "double y_init =%f; \n", ystart);
    fprintf(fp, "double y_0 = y_init;  \n");
    fprintf(fp, "double total_fxn_evals =0;  \n");
    fprintf(fp,"double y_1, y_2, y_3, y_4; \n");
    fprintf(fp, "double stepsize = %0.12f;\n", step_size);
    fprintf(fp, "double tim = (double) 0;\n"); 


    fprintf(fp, "for(int i =1; i<5;i++){ \n");
    fprintf(fp, "tim = tim+ %0.12f;\n", step_size);
        
        fprintf(fp, "double k1 =func(tim, y_init);  \n");
        fprintf(fp, "double k2 = func(tim + stepsize/2.0, y_init + (0.5*stepsize*k1)); \n");
        fprintf(fp, "double k3 = func(tim + stepsize/2.0, y_init + (0.5 * stepsize *k2)); \n");
        fprintf(fp, "double k4 = func(tim + stepsize, y_init +( stepsize *k3));  \n");
        fprintf(fp, "y_init =y_init + (double)1.0/6.0*stepsize*(k1+2*k2 +2*k3 +k4);  \n");
        fprintf(fp, "total_fxn_evals +=4;  \n");

        fprintf(fp, "errorCalc2[i] = fabs(y_init-exactArray[i]); \n"); 
	
        fprintf(fp, "if (i==1){ \n");
                fprintf(fp, "y_1 = y_init;  \n");

        fprintf(fp, "} \n");
        fprintf(fp, "if (i==2){ \n");
                fprintf(fp, "y_2 = y_init;  \n");
        fprintf(fp, "} \n");
        fprintf(fp, "if (i==3){ \n");
                fprintf(fp, "y_3 = y_init; \n");
                        
        fprintf(fp, "} \n");
        printf(fp, "if (i==4){ \n");
                fprintf(fp, "y_4 = y_init; \n");
                        
        fprintf(fp, "} \n");

	fprintf(fp, "double tim_0 = 0.0;\n");
	fprintf(fp, "double tim_1 = %0.12f;\n", 1*step_size);
	fprintf(fp, "double tim_2 = %0.12f;\n", 2*step_size);
    fprintf(fp, "double tim_3 = %0.12f;\n", 3*step_size);
    fprintf(fp, "double tim_4 = %0.12f;\n", 4*step_size);

       fprintf(fp, "for (int i= 5; i <= %d; i++) {\n", num_intervals);

        fprintf(fp, "y_4 = y_init;  \n");

        fprintf(fp, "y_init= y_4+ 1.0 /720.0 * stepsize * ( 1901.0 * func(tim_4,y_4) -2774.0 * func(tim_3,y_3) + 2616.0 *func(tim_2, y_2) - 1274.0 *func(tim_1,y_1) +251.0 *func(tim_0,y_0)); \n");
        fprintf(fp, "y_0 = y_1;  \n");
        fprintf(fp, "y_1= y_2;  \n");
        fprintf(fp, "y_2= y_3;  \n");
        fprintf(fp, "y_3 = y_4;  \n");
        fprintf(fp, "y_4 = y_init;  \n");


	
	
    fprintf(fp, "tim_0 +=%0.12f;\n", step_size);
    fprintf(fp, "tim_1 +=%0.12f;\n", step_size);
    fprintf(fp, "tim_2 +=%0.12f;\n", step_size);
    fprintf(fp, "tim_3 +=%0.12f;\n", step_size);
    fprintf(fp, "tim_4 +=%0.12f;\n", step_size);

        fprintf(fp, "total_fxn_evals +=5; \n");
		
	
        fprintf(fp, "} \n");
        fprintf(fp, "sleep(1);\n");

        fprintf(fp, "printf(\"The total function evaluations by Adams-Bashforth were %%0.2f \\n\",   total_fxn_evals);\n");

        fprintf(fp, "printf(\"The error after %d intervals was %%0.15f \\n\",   y_init-exactArray[%d]);\n", num_intervals, num_intervals);
        fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");
        fprintf(fp,"\n");
       fprintf(fp, "return errorCalc2;\n");
        fprintf(fp, "}	 \n");

    //constructing RKF8 subroutine 

    fprintf(fp, "double * testRKF8(){ \n");
fprintf(fp, "double y_init = %f; \n", ystart);

fprintf(fp, "double static errorCalc3[%d+1];\n", num_intervals);

fprintf(fp, "errorCalc3[0] = 0; \n"); 

        

fprintf(fp, "double total_fxn_evals =0;  \n");

fprintf(fp, "double stepsize = %0.12f;\n", step_size);
fprintf(fp, "double tim = (double) 0;\n"); 

fprintf(fp, "for (int i= 1; i <= %d; i++) {\n", num_intervals);

fprintf(fp, "tim = tim+ %0.12f;\n", step_size);




        fprintf(fp, "double k1 =func(tim, y_init);  \n");
        fprintf(fp, "double k2 = func(tim +2*stepsize/27.0, y_init + stepsize *2*(k1/27.0));  \n");
        fprintf(fp, "double k3 = func(tim+3*stepsize/9.0, y_init  +  stepsize /((k1/36.0 + k2/12.0)));  \n");

        fprintf(fp, "double k4 =func(tim +stepsize/6.0, y_init + stepsize * (k1/24.0 + k3/8.0));  \n");
        fprintf(fp, "double k5 = func(tim + 5.0/12.0*stepsize, y_init + stepsize *(5.0/12.0*k1 + -25.0/16.0*k3+  25.0/16.0*k4)); \n");
        fprintf(fp, "double k6 = func(tim + stepsize/2.0, y_init + stepsize * (k1/20.0 + k4/4.0 + k5/5.0));  \n");
        fprintf(fp, "double k7 = func(tim + 5.0/6.0*stepsize, y_init + stepsize *(-25.0/108.0*k1 + 125.0/108.0*k4 - 65.0/27.0*k5 + 125.0/54.0*k6)); \n");

        fprintf(fp, "double k8 = func(tim +stepsize/6.0, y_init + stepsize * (31.0/300.0*k1 + 61.0/225.0*k5 - 2.0/9.0*k6 + 13.0/900.0*k7)); \n");
        fprintf(fp, "double k9 = func(tim + 2.0/3.0*stepsize, y_init  + stepsize *(2.0*k1 + -53.0/6.0*k4 + 704.0/45.0*k5 + -107.0/9.0*k6 + 67.0/90.0*k7 + 3.0*k8)); \n");
        fprintf(fp, "double k10 = func(tim + stepsize/3.0, y_init + stepsize *(-91.0/108.0*k1 + 23/108.0*k4 + -976.0/135.0*k5 + 311.0/54.0*k6 + -19.0/60.0*k7 + 17.0/6.0*k8 - k9/12.0)); \n");

        fprintf(fp, "double k12 = func(tim, y_init  + stepsize * (3.0/205.0*k1 + -6.0/41.0*k6 + -3.0/205.0*k7 + -3.0/41.0*k8 + 3.0/41.0*k9 + 6.0/41.0*k10)); \n");
        fprintf(fp, "double k13 = func(tim + stepsize, y_init + stepsize *(-1777.0/4100.0*k1 + -341.0/164.0*k4 + 4496.0/1025.0*k5 + -289.0/82.0*k6 + 2193.0/4100.0*k7 + 51.0/82.0*k8 + 33.0/164.0*k9 + 12.0/41.0*k10 + k12));  \n");

	    fprintf(fp, "total_fxn_evals +=12;  \n");

        fprintf(fp, "y_init =y_init + stepsize*(34.0/105.0*k6 + 9.0/35.0*k7 + 9.0/35.0*k8 + 9.0/280.0*k9 + 9.0/280.0*k10 + 41.0/840.0*k12 + 41.0/840.0*k13); \n");
        fprintf(fp, "errorCalc3[i] = fabs(y_init-exactArray[i]); \n"); 

        fprintf(fp, "} \n");
        fprintf(fp, "sleep(1);\n");

        fprintf(fp, "printf(\"The total function evaluations by RK8 were %%0.2f \\n\",   total_fxn_evals);\n");

        fprintf(fp, "printf(\"The error after %d intervals was %%0.15f \\n\",   y_init-exactArray[%d]);\n", num_intervals, num_intervals);
        fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");
        fprintf(fp,"\n");
         fprintf(fp, "return errorCalc3;\n");


    fprintf(fp, "} \n");

    fprintf(fp,"\n");

    //The main method of the ODE function

    fprintf(fp, "int main(){\n");

    fprintf(fp, "getExactSolution();\n"); // calculates exact solution 

    //define variables for error
    fprintf(fp, "double * euler_error;\n");
    fprintf(fp, "double * rk4_error;\n");
    fprintf(fp, "double * AB_error;\n");
    fprintf(fp, "double * rkf8_error;\n");


    // run methods and calculate their efficacy
    fprintf(fp, "euler_error = testForwardEuler(); \n");

    fprintf(fp, "double euler_efficacy = calcEfficacy(euler_error); \n");
    fprintf(fp, "printf(\"The Efficacy for Forward Euler was %%0.8f\\n\",euler_efficacy); \n");
    fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");


    fprintf(fp, "rk4_error = testRK4(); \n");

    fprintf(fp, "double rk4_efficacy = calcEfficacy(rk4_error); \n");
    fprintf(fp, "printf(\"The Efficacy for RK4 was %%0.8f\\n\",rk4_efficacy); \n");
    fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");

    fprintf(fp, " AB_error = testAdamsBashforth(); \n");

    fprintf(fp, "double AB_efficacy = calcEfficacy(AB_error); \n");
    fprintf(fp, "printf(\"The Efficacy for Adams-Bashforth was %%0.8f\\n\",AB_efficacy); \n");
    fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");

    fprintf(fp, "rkf8_error = testRKF8(); \n");

    fprintf(fp, "double rkf8_efficacy = calcEfficacy(rkf8_error); \n");
    fprintf(fp, "printf(\"The Efficacy for RKF8 was %%0.8f\\n\",rkf8_efficacy); \n");
    fprintf(fp, "printf(\"-----------------------------------------------------------------\\n\");\n");
    
    
    //An array to keep track of the efficacy values of each method 

     fprintf(fp, "double efficacyArray[] = {euler_efficacy, rkf8_efficacy, AB_efficacy, rk4_efficacy}; \n"); 


    //Figure out the best efficacy and then print that method
     fprintf(fp, "double bestEfficacy = -1000; \n"); 
     fprintf(fp, "for (int i =0; i<4; i++){\n"); 
         fprintf(fp, "if (efficacyArray[i] > bestEfficacy){\n"); 
             fprintf(fp, "bestEfficacy = efficacyArray[i]; \n"); 
             fprintf(fp, "continue;\n"); 
         fprintf(fp, "}\n"); 
     fprintf(fp, "}\n"); 
     fprintf(fp, "if (bestEfficacy == efficacyArray[0]){\n"); 
         fprintf(fp, "printf(\"The method you should use is Forward Euler\\n\");\n"); 
     fprintf(fp, "}\n"); 
    
     fprintf(fp, "if (bestEfficacy == efficacyArray[1]){\n"); 
         fprintf(fp, "printf(\"The method you should use is RKF8\\n\");\n"); 
     fprintf(fp, "}\n"); 

      fprintf(fp, "if (bestEfficacy == efficacyArray[2]){\n"); 
         fprintf(fp, "printf(\"The method you should use is Adams-Bashforth\\n\");\n"); 
     fprintf(fp, "}\n"); 

      fprintf(fp, "if (bestEfficacy == efficacyArray[3]){\n"); 
         fprintf(fp, "printf(\"The method you should use is RK4\\n\");\n"); 
     fprintf(fp, "}\n"); 



    fprintf(fp, "return 0;\n"); 

    fprintf(fp, "}\n");


    fclose(fp);


    return 0; 
}
/*
 * This method creates a makefile to then compile the newly generated ODE file
 *
 * 
*/
void createMakefile(){

    char *filename = "Makefile";
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Failed to open file %s\n", filename);
        exit(1);
    }

    fprintf(fp, "include ../Make.inc\n");
    fprintf(fp, "\n");

    fprintf(fp, "CFLAGS-add += -I../src -DNDEBUG -w\n");
    fprintf(fp, "LDFLAGS-add += -L../src -llsoda -lm\n");
    fprintf(fp, "\n");
    fprintf(fp, "first: textUI\n");
    fprintf(fp, "second: ode\n");
    fprintf(fp, "\n");

    fprintf(fp, "textUI.o: textUI.c\n");
    fprintf(fp, "\t$(CC) $(CPPFLAGS) $(CFLAGS-add) $(CFLAGS) $(fPIC) -c $< -o $@\n");
    fprintf(fp, "textUI: textUI.o\n");
    fprintf(fp, "\t$(CC) $(LDFLAGS) -o $@ $^ $(LDFLAGS-add)\n");

    fprintf(fp, "ode.o: ode.c\n");
    fprintf(fp, "\t$(CC) $(CPPFLAGS) $(CFLAGS-add) $(CFLAGS) $(fPIC) -c $< -o $@\n");
    fprintf(fp, "ode: ode.o\n");
    fprintf(fp, "\t$(CC) $(LDFLAGS) -o $@ $^ $(LDFLAGS-add)\n");

    fprintf(fp, "clean:\n");
    fprintf(fp, "\trm -f *.o ode textUI ysin1y ode\n");

    
    fclose(fp);
    char *make_command = "make second";
    system(make_command);

    //unlink(filename);


}

int main(){

    printf("////////////////////////Welcome to the Methodizer///////////////////////////");
    printf("\n"); 
    printf("Your Peronsal 1-D ODE solver\n");

   printf("Verison 1.0\n");
   sleep(1);
   printf("Here is some info about how to use this program\n");
   sleep(1); 
   

   printf("* The fucntions need to be written in the form y' = t * y etc. where there is a space bewteen every opperation\n");
  sleep(1); 
   printf("* The methods can only handle first order ODE's\n");
   sleep(1); 

   printf("* Also every expontial and trig and power funtions will have to be written in the form pow(t,3) and tan(t)\n"); 
    //printf("* Buffer every variable with whitespace\n");
    printf("* The program cannot handle divide by 0 so you cannot pass in a function that divides by t\n");

   printf("Here are some examples\n");
   sleep(1); 

   printf("tan(t) + 7 * y + 6\n");
   printf("6 * pow(t,3) + 54 * y\n");
   sleep(1); 

   printf("To exit the program hit cntrl z at any time\n");
   
   

   // assume buffer size is less than 150 
   
   char str[150];
    char h[100];
    char init[20]; 

   //while (*truthval == 'n'){
      printf("Enter your funcition y'=");
      
      
      if (fgets(str, 150, stdin) ==NULL){
        printf("Error reading string"); 
        exit(1); 
      }
      char orig[150]; 
      strcpy(orig,str); 

        replace_y_with_y0(str);

        printf("Choose an inintial value for y\n");
        printf("ie y[0] =");

        if (fgets(init, 150, stdin) ==NULL){
        printf("Error reading string\n"); 
        exit(1); 
        }


        printf("Choose a stepsize from 0 to 1:\n");

        if (fgets(h, 150, stdin) ==NULL){
        printf("Error reading string\n"); 
        exit(1); 
        }

        printf("The function you have chosen is: %s\n", orig);
         printf("with intial value y[0] = %s\n", init);

        printf("The step size you have chosen is: %s\n", h);
        
        double ystart = strtod(init, NULL); 

      createFile(str,orig,h, ystart);  // eventually add in step size 

      createMakefile(); 

    char *run_command = "./ode";
    system(run_command);

    printf("Thank you for using One Dimentional ODE solver\n"); 
   
   

   return 0;
}
