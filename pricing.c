#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double randn(double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

double stabledevd(float alpha) {

  // we define beta = 0
  // for a symmetric stable distribution. 

  double theta, W, holder, left, right;

  theta=M_PI*((float)rand()/RAND_MAX - 0.5);
  W = -log((float)rand()/RAND_MAX); // takes natural log

  // if beta = 0 then c(alpha,beta)=1; theta_0 = 0
  // expression reduces to sin alpha.theta / (cos theta) ^1/alpha
  //  * (cos (theta - alpha theta)/W) ^(1-alpha)/alpha

  left = (sin(alpha*theta)/pow(cos(theta), 1.0/alpha));
  right= pow(cos(theta*(1.0 - alpha))/W, ((1.0-alpha)/alpha));
  holder=left*right;
  return(holder);
}

float path(float K, float S_0, float sigma_S, float lambda, int length) {
    srand(time(NULL));
    float *x = malloc(length * sizeof(float));
    float *omega = malloc((length+1) * sizeof(float));
    float *S = malloc((length+1) * sizeof(float));

    if (x == NULL || omega == NULL || S == NULL) {
        fprintf(stderr, "malloc failed\n");
        return -1;
    }
    
    omega[0] = 0;
    for (int i=0; i<length; i++) {
        x[i] = stabledevd(1.7);
        omega[i+1] = omega[i]+x[i]*0.0001;
        /*  printf("%.2f\n",omega[i]); */
    }
    free(x);

    S[0] = S_0;
    for (int j=0; j<length; j++) {
        S[j+1] = S[j]+omega[j]*(lambda-S[j])+S[j]*sigma_S*randn(0,1);
    }
    
    float payoff;

    if (S[length]-K > 0) {
       payoff = S[length]-K;
    }
    else {
       payoff = 0;
    }

    free(omega);
    free(S);

    return payoff;
}

int main(void) {

    int iterations = 1000000;
    float *sum = malloc(iterations * sizeof(float));

    if (sum == NULL) {
        fprintf(stderr, "malloc failed\n");
        return -1;
    }
    
    FILE *f = fopen("call_price.txt", "w");
    if (f == NULL) {
       printf("Error opening file!\n");
       return -1;
    }

    sum[0] = 0;
    for (int k = 1; k <= iterations; k++) {
        sum[k] = sum[k-1] + path(10.0, 10.0, 0.02, 0, 90);
        fprintf(f,"%.2f\n",sum[k]/k);
    }

    return 0;
}
