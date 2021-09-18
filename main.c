#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double stabledevd(float alpha) {

  // we define beta = 0
  // for a symmetric stable distribution. 

  double theta, W, holder, left, right;

  theta = M_PI*((float)rand()/RAND_MAX - 0.5);
  W = -log((float)rand()/RAND_MAX); // takes natural log

  // if beta = 0 then c(alpha,beta)=1; theta_0 = 0
  // expression reduces to sin alpha.theta / (cos theta) ^1/alpha
  //  * (cos (theta - alpha theta)/W) ^(1-alpha)/alpha

  left = (sin(alpha*theta)/pow(cos(theta), 1.0/alpha));
  right = pow(cos(theta*(1.0 - alpha))/W, ((1.0-alpha)/alpha));
  holder = left*right;
  return(holder);
}

int main(void) {
    srand(time(NULL));
    int length = 10000;
    float *x = malloc(length * sizeof(float));
    if (x == NULL) {
        fprintf(stderr, "malloc failed\n");
        return -1;
    }
    
    for (int i=0; i<length; i++) {
        x[i] = stabledevd(1.6);
        printf("%.2f\n",x[i]);
    }

    free(x);
    return 0;
}
