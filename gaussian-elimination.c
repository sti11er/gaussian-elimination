#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <unistd.h>
#include <sys/wait.h>
#define N 10
# define EPS 0.0000001


void 
initMatrix(double** A1, double *B1,double** A2, double *B2){
    double buf;
    for (int i=1; i<N; i++){
        for (int j=1; j<N; j++){
            buf=(double)1/(i+1);
            A1[i-1][j-1]=pow(buf,j-1);
        }
    }
    for (int i=1; i<N+2; i++){
        for (int j=1; j<N+2; j++){
            if (j==(N+1) || i==j){
                A2[i-1][j-1]=(double)1;
            }
            else if (j>i) {
                A2[i-1][j-1]=0;
            }
            else {
                A2[i-1][j-1]=j-i;
                }
            }
        }
    for (int i=1; i<N; i++){
        buf=(double)1/(i+1);
        B1[i-1]=pow(buf,N-2);
    }
    for (int i=1; i<N+2; i++){
        B2[i-1]=(double)(1.5);
    }

}

void Nevyazka(double** A, double *f, double *x, unsigned size){
  for (int i = 0; i < size; i++)
  {
      f[i] = 0-f[i];
      for (int k = 0; k < size; k++){
        f[i] += A[i][k] * x[k];
    }
  }

}

double Norma(double* f, unsigned size){
  double s=0;
  for (int i = 0; i < size; i++)
  {
    s+=f[i]*f[i];
  }
  s=pow(s, 0.5);
  return s;
}

// Функция, которая меняет местами две строки матрицы
void 
swap_rows(double **A, double *b, int i, int j) {
    double *temp_row = A[i];
    A[i] = A[j];
    A[j] = temp_row;
    double temp_b = b[i];
    b[i] = b[j];
    b[j] = temp_b;
}

void 
gaussian_elimination(double **A,  double *b, double *x, unsigned size, int step) {
    int i, j;
    double buf;

    if (size == 1) {
        x[step] = b[step] / A[step][step];
        return;
    }

    // Находим индекс строки с максимальным элементом в i-м столбце
    int i_max = step;
    for (int j = step+1; j < (size+step); j++) {
        if (fabs(A[j][step]) > fabs(A[i_max][step])) {
            i_max = j;
        }
    }
    // Меняем местами i-ю и i_max-ю строки
    swap_rows(A, b, step, i_max);
    // Проверяем, что главный элемент не нулевой
    if (A[step][step] == 0) {
        printf("Матрица вырождена или несовместна\n");
        exit(1);
    }

    b[step] /= A[step][step];
    for (i = 1; i < size; i++){
        A[step][step + i] /= A[step][step];
    }
    A[step][step] = 1;

    for (i = 1; i < size; i ++){
        b[i + step] -= b[step] * A[i + step][step];
        for (j = 1; j < size; j ++){
            A[i + step][j + step] -= A[step][j + step] * A[i + step][step];
        }
        A[i + step][step] = 0;
    }
    gaussian_elimination(A, b, x, size - 1, step + 1);


    x[step] = b[step];
    for (i = 1; i < size; i ++){
        x[step] -= x[i + step] * A[step][i + step];
    }
    x[step] /= A[step][step];
    return;
}


int 
main(int argc, char* argv[]) {
    double ** A1 = calloc((N - 1), sizeof(double *));
    for (int i = 0; i < N; i ++) {
        A1[i] = calloc((N-1), sizeof(double));
    }
    double ** A2 = calloc((N + 1), sizeof(double *));
    for (int i = 0; i <= N+1; i ++) {
        A2[i] = calloc((N+1), sizeof(double));
    }
    
    double *B1=calloc((N-1), sizeof(double));
    double *B2=calloc((N+1), sizeof(double));
    
    double *X1=calloc((N-1), sizeof(double));
    double *X2=calloc((N+1), sizeof(double));


    initMatrix(A1, B1, A2, B2);

    for (int i = 0; i < N+1; i ++) {
        for (int j = 0; j < N+1; j ++) {
            printf("%.3lf ", A2[i][j]);
        }
        printf("| %.3lf ", B2[i]);
        printf("\n");
    }
    printf("\n");

    gaussian_elimination(A2,B2,X2,N+1,0);
    //Gauss(A2,B2,X2,M+1,0);
    
    printf("Gaussian Elimination\n");
    for (int i = 0; i < N+1; i ++) {
        for (int j = 0; j < N+1; j ++) {
            printf("%.3lf ", A2[i][j]);
        }
        printf("| %.3lf ", B2[i]);
        printf("\n");
    }

    // printf("\n");
    // for (int i = 0; i < N+1; i ++) {
    //     printf("%.3lf ", X2[i]);
    // }
    // printf("\n");
    // printf("\n");

    // Nevyazka(A2, B2, X2, N+1);
    // printf("NORMA NEVYAZKI ");
    // printf("%.100lf", Norma(B2,N+1));

    // initMatrix(A1, B1, A2, B2);
    // Gauss(A2,B2,X2,N+1,0);
    // printf("\nGauss\n");
    // for (int i = 0; i < N+1; i ++) {
    //     for (int j = 0; j < N+1; j ++) {
    //         printf("%.3lf ", A2[i][j]);
    //     }
    //     printf("| %.3lf ", B2[i]);
    //     printf("\n");
    // }

    // printf("\n");
    // for (int i = 0; i < N+1; i ++) {
    //     printf("%.3lf ", X2[i]);
    // }
    // printf("\n");
    // printf("\n");

    // Nevyazka(A2, B2, X2, N+1);
    // printf("NORMA NEVYAZKI ");
    // printf("%.100lf", Norma(B2,N+1));

    return 0;
}