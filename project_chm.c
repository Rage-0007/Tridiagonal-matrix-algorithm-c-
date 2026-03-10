#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Функция решения СЛАУ с трехдиагональной матрицей методом прогонки
// A - массив диагонали под главной (размер n-1)
// B - массив главной диагонали (размер n)
// C - массив диагонали над главной (размер n-1)
// F - массив правых частей (размер n)
// X - массив решений (размер n)
// n - размер системы
int thomas_algorithm(double *A, double *B, double *C, double *F, double *X, int n) {
    if (n <= 0) return -1;
    
    // Массивы для прогоночных коэффициентов
    double *alpha = (double*)malloc((n) * sizeof(double));
    double *beta = (double*)malloc((n) * sizeof(double));
    
    if (!alpha || !beta) {
        free(alpha);
        free(beta);
        return -1;
    }
    
    // Прямой ход прогонки
    // Для первого уравнения
    double denominator = B[0];
    if (fabs(denominator) < 1e-15) {
        free(alpha);
        free(beta);
        return -1; // Вырожденная матрица
    }
    
    alpha[0] = -C[0] / denominator;
    beta[0] = F[0] / denominator;
    
    // Для остальных уравнений i = 1..n-2
    for (int i = 1; i < n - 1; i++) {
        denominator = B[i] + A[i-1] * alpha[i-1];
        if (fabs(denominator) < 1e-15) {
            free(alpha);
            free(beta);
            return -1; // Вырожденная матрица
        }
        
        alpha[i] = -C[i] / denominator;
        beta[i] = (F[i] - A[i-1] * beta[i-1]) / denominator;
    }
    
    // Для последнего уравнения
    denominator = B[n-1] + A[n-2] * alpha[n-2];
    if (fabs(denominator) < 1e-15) {
        free(alpha);
        free(beta);
        return -1; // Вырожденная матрица
    }
    
    beta[n-1] = (F[n-1] - A[n-2] * beta[n-2]) / denominator;
    
    // Обратный ход прогонки
    X[n-1] = beta[n-1];
    
    for (int i = n - 2; i >= 0; i--) {
        X[i] = alpha[i] * X[i+1] + beta[i];
    }
    
    free(alpha);
    free(beta);
    return 0;
}

// Пример использования
int main() {
    int n = 5; // Размер системы
    
    // Выделение памяти для матрицы и правых частей
    double *A = (double*)malloc((n-1) * sizeof(double)); // Нижняя диагональ
    double *B = (double*)malloc(n * sizeof(double));     // Главная диагональ
    double *C = (double*)malloc((n-1) * sizeof(double)); // Верхняя диагональ
    double *F = (double*)malloc(n * sizeof(double));     // Правая часть
    double *X = (double*)malloc(n * sizeof(double));     // Решение
    
    if (!A || !B || !C || !F || !X) {
        printf("Ошибка выделения памяти\n");
        return 1;
    }
    
    // Заполнение нижней диагоналей
    for (int i = 0; i < n-1; i++) {
        A[i] = 1.0; //  диагональ
    }

    // Заполнение верхнее диагоналей
    for (int i = 0; i < n-1; i++) {
        C[i] = 1.0;
    }

    // Главная диагональ
    for (int i = 0; i < n; i++) {
        B[i] = 4.0;
    }
    
    // Правая часть
    for (int i = 0; i < n; i++) {
        F[i] = i * 0.02;
    }

    // Вывод системы
    printf("Трехдиагональная система уравнений:\n");
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            printf("%.1f*x%d + ", A[i-1], i);
        }
        
        printf("%.1f*x%d", B[i], i+1);
        
        if (i < n-1) {
            printf(" + %.1f*x%d", C[i], i+2);
        }
        
        printf(" = %.2f\n", F[i]);
    }
    
    // Решение системы
    int result = thomas_algorithm(A, B, C, F, X, n);
    
    if (result == 0) {
        // Вывод решения
        printf("\nРешение системы:\n");
        for (int i = 0; i < n; i++) {
            printf("x%d = %.6f\n", i+1, X[i]);
        }
    } else {
        printf("Ошибка при решении системы (возможно, вырожденная матрица)\n");
    }
    
    // Освобождение памяти
    free(A);
    free(B);
    free(C);
    free(F);
    free(X);
    
    return 0;
}