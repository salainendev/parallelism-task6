#include <iostream>
#include <boost/program_options.hpp>
#include <cmath>



double linearInterpolation(double x, double x1, double y1, double x2, double y2) {
    // делаем значение y(щначение клетки)используя формулу линейной интерполяции
    return y1 + ((x - x1) * (y2 - y1) / (x2 - x1));
}

int main(int argc, char const *argv[])
{
    int N = 128;

    std::unique_ptr<double[]> matrix(new double[128*128]); // матрица температуры
    
    // заполнение матрицы нулями
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            matrix[i*N + j] = 0;
        }
        
    }
    // заполнение углов матрицы по очереди ( как я понял по часовой стрелке)
    matrix[0] = 10.0;
    matrix[N-1] = 20.0;
    matrix[(N-1)*N + (N-1)] =30.0;
    matrix[(N-1)*N] = 20.0; 
    // линейная интерполяция между углами матриц 
    
    // верхняя граница 
    for (size_t i = 1; i < N-1; i++)
    {
        matrix[0*N+i] = linearInterpolation(i,0.0,matrix[0],N-1,matrix[N-1]);
    }
    // left граница 
    for (size_t i = 1; i < N-1; i++)
    {
        matrix[i*N+0] = linearInterpolation(i,0.0,matrix[0],N-1,matrix[(N-1)*N]);
    }
    // right граница

    for (size_t i = 0; i < N-1; i++)
    {
        matrix[i*N+(N-1)] = linearInterpolation(i,0.0,matrix[N-1],N-1,matrix[(N-1)*N + (N-1)]);
    }
    
    // нижняя граница 
    for (size_t i = 1; i < N-1; i++)
    {
        matrix[(N-1)*N+i] = linearInterpolation(i,0.0,matrix[(N-1)*N],N-1,matrix[(N-1)*N + (N-1)]);
    }
    

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            std::cout << matrix[i*N+j] << " ";
        }
        std::cout << std::endl;
    }
    
    

    return 0;
}
