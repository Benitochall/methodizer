#include <time.h>

// A simple program to test the time it takes for a function to add 1000000
// numbers in c


int main(){
    int total = 0;
    clock_t start_t, end_t; // setting clock values 
    double total_t;
    start_t = clock(); 



for (int i=0; i< 1000000; i++){
    total+=i; 
}
end_t = clock();
total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC *10E6;
printf("The total time taken by in microseconds was: %f\n", total_t);


}
