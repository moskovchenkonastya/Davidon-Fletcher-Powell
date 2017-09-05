//
//  main.cpp
//  main
//
//  Created by Anastasiya Moskovchenko on 01.11.16.
//  Copyright © 2016 Anastasiya Moskovchenko. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <cstdio>
#include <cmath>


#define eps 0.000001
#define x0 1.5
#define x1 1.2
#define x2 0.4
#define x3 0.3
#define x4 0.1

using namespace std;

double my_function(vector<double> &x){
    //return x[0] + x[1];
    double sum = 0;
    for (int i = 1; i < 5; ++i)
        sum += pow((x[i] - (i + 1) * x[0]), 4);
    
    return 150 * sum + pow(x[0] - 1, 2);
}

double barrier_function(vector<double> &x, double z){
    
    double sum = 0;
    for (int i = 0; i < 5; ++i)
        sum += (i + 1) * pow(x[i], 2);
    return 1/(z - sum);
    
}

double new_function(vector<double> &x, int k, double z){
    return my_function(x) + (1.0/(k * k)) * barrier_function(x, z);
}

vector<double> func_grad(vector<double> &x, int k, double z, const double delta = 0.000001){
    
    vector<double> grad(x.size());
    vector<double> curr(x.size());
    curr = x;
    for (int i = 0; i < x.size(); ++i){
        curr[i] += delta;
        grad[i]  = (new_function(curr, k, z) - new_function(x, k, z)) / delta;
        curr[i] -= delta;
    }
    
    return grad;
}

vector<double> vector_alpha(vector<double> &x, double alpha){
    vector<double> result(x.size());
    for (int i = 0; i < x.size(); ++i)
        result[i] = x[i] * alpha;
    return result;
}

vector<double> vector_minus(vector<double> &x, vector<double> &y){
    vector<double> result(x.size());
    for (int i = 0; i < x.size(); ++i)
        result[i] = x[i] - y[i];
    return result;
}

vector<double> vector_add(vector<double> &x, vector<double> &y){
    vector<double> result(x.size());
    for (int i = 0; i < x.size(); ++i)
        result[i] = x[i] + y[i];
    return result;
}


double vector_scalar(vector<double> &x, vector<double> &y){
    double pr = 0;
    for (int i = 0; i < x.size(); ++i){
        pr += x[i] * y[i];
    }
    return pr;
}

vector<vector<double> > vector_multipl(vector<double> &x, vector<double> &y){
    vector<vector<double> > M(x.size(), vector<double>(x.size()));
    for (int i = 0; i < x.size(); ++i){
        for (int j = 0; j < x.size(); ++j){
            M[i][j] = x[i] * y[j];
        }
    }
    return M;
}

vector<vector<double> > vector_multipl_T(vector<double> &x, vector<double> &y){
    vector<vector<double> > M(x.size(), vector<double>(x.size()));
    for (int i = 0; i < x.size(); ++i){
        for (int j = 0; j < x.size(); ++j){
            M[i][j] = x[i] * y[j];
        }
    }
    return M;
}



vector<vector<double> > divide_matrix(vector<vector<double> > &M, double alpha){
    vector<vector<double> > result(M.size(), vector<double>(M.size()));
    for (int i  = 0; i < M.size(); ++i){
        for (int j = 0; j < M.size(); ++j){
            result[i][j] = M[i][j] / alpha;
        }
    }
    return result;
}

vector<vector<double> > sub_matrix(vector<vector<double> > &M1, vector<vector<double> > &M2){
    vector<vector<double> > result(M1.size(), vector<double>(M1.size()));
    for (int i = 0; i < M1.size(); ++i)
        for (int j = 0; j < M1.size(); ++j)
            result[i][j] = M1[i][j] - M2[i][j];
    return result;
}

vector<vector<double> > sum_matrix(vector<vector<double> > &M1, vector<vector<double> > &M2){
    vector<vector<double> > result(M1.size(), vector<double>(M1.size()));
    for (int i = 0; i < M1.size(); ++i)
        for (int j = 0; j < M1.size(); ++j)
            result[i][j]  = M1[i][j] + M2[i][j];
    return result;
}
vector<double> matrix_vector(vector<vector<double> > &M, vector<double> &x){
    vector<double> result(x.size(),0);
    for (int i = 0; i < x.size(); ++i){
        for (int j = 0; j< x.size(); ++j){
            result[i] += M[i][j] * x[j];
        }
    }
    return result;
}


vector<vector<double> > copy_matrix(vector<vector<double> > &M){
    vector<vector<double> > result(M.size(), vector<double>(M.size()));
    
    for (int i = 0; i < M.size(); ++i){
        
        result[i] = M[i];
    }
    return result;
}

double find_min_new(vector<double> &x, vector<double> &d, double a, double b, int count, double  z){
    if(fabs(b - a) < eps)
        return (a + b) / 2.0;
    
    
    double e = 0;
    double x_1 = 0;
    double x_2 = 0;
    vector<double> V_a1;
    vector<double> V_a2;
    
    
    e = (b - a) * 0.0001; // 1E-2 * 1E-3;
    x_1 = (b + a) / 2.0 - e;
    x_2 = (b + a) / 2.0 + e;
        
    V_a1 = vector_alpha(d, x_1);
    V_a2 = vector_alpha(d, x_2);
        
    V_a1 = vector_add(x, V_a1);
    V_a2 = vector_add(x, V_a2);
        
    double f1 = new_function(V_a1, count, z);
    double f2 = new_function(V_a2, count, z);
    
    if (fabs(f1) > fabs(f2)){
        return find_min_new(x, d, (a + b) /2.0, b, count, z);
        }else{
            return find_min_new(x, d, a, (a + b) /2.0, count, z);
    }

    
    return 0.0;
    
}




int main( int argc, char *argv[]) {
    
  /*
    if (argc < 6) {
        printf("Error: Некорректное количество переменных. Введите значения eps и 5 переменных.\n");
        return 0;
    }
  
    double x0, x1, x2, x3, x4, x5, eps;

    
    if(argc == 0){
        
        eps = 0.000001;
        x0  = 0.1;
        x1  = 1.5;
        x2  = 0.4;
        x3  = 1.2;
        x4  = 1.3;
    
    }
    

    x0 = atof(argv[0]);
    x1 = atof(argv[1]);
    x2 = atof(argv[2]);
    x3 = atof(argv[3]);
    x4 = atof(argv[4]);
    x5 = atof(argv[5]);
    
    
    eps = atof(argv[6]);
   
    printf("Kоординаты начальной точки: X_0 (%s, %s, %s, %s, %s, %s)\n",argv[0], argv[1], argv[2], argv[3], argv[4]);
    */
    
    double z =224.0;
    do{
    double alpha;

    int count = 0;

    double vec_scalar = 0;
    

    vector<double> V_a1;
    vector<double> V_a2;
    
    vector<double> gradient;
    vector<vector<double> > matrix;
    //vector<vector<double> > tmp;
    const double arr0[] = {0, 0, 0, 0, 0};
    vector<double> x;
    x.assign(arr0, arr0 + sizeof arr0 / sizeof * arr0);
    
    vector<double> y;
    y.assign(arr0, arr0 + sizeof arr0 / sizeof * arr0);
    
    vector<double> d;
    d.assign(arr0, arr0 + sizeof arr0  / sizeof * arr0);
    
    double norm;
   
    
    vector<double> old_x = x;
    vector<double> old_grad = gradient;
    
    vector<double> delta_x(x.size());
    vector<double> delta_d(x.size());
    vector<double> m_v(x.size());
    
   
    // initial estimate x0
    const double arr1[] = {x0, x1, x2, x3, x4};
    vector<double> values;
    values.assign(arr1, arr1 + sizeof arr1 / sizeof * arr1);
    
    // проверка принадлежности начальной точки к множеству
    double sum = 0;
    for (int i = 0; i < 5; ++i)
        sum += pow(i * values[i], 2);
    
    if(sum > 224){
        cout<<sum;
        cout<<"Начальная точка не принадлежит множеству ограничений"<<endl;
        return 0;
    }

    // select a positive definite matrix S0
    vector<vector<double> > vec_mult(values.size(), vector<double>(values.size()));
    vector<vector<double> > S(values.size(), vector<double>(values.size()));
    vector<vector<double> > dev_S(values.size(), vector<double>(values.size()));
    vector<vector<double> > tmp(values.size(), vector<double>(values.size()));
    
    double cur_f, prev_f;

    count = 1;
    vector<double> X0;
    
    do{
        //++count;
        x = values;
        X0 = values;

        
        //
        for (int i = 0; i < values.size(); ++i){
            for (int j = 0; j < values.size(); ++j){
                if (i == j){
                    S[i][j] = 1.0;
                }else{
                    S[i][j] = 0.0;
                }
            }
        }
        
        
        do{
            
            tmp = copy_matrix(S);
     
            // compulate gradient
            gradient = func_grad(x, count, z);
            
            old_grad = gradient;
            
            // compulate d
            d = matrix_vector(S, gradient);
            d = vector_alpha(d, -1.0);
            
            // looking for the minimum of F(x - alpha * d) using the method of one-dimensional optimization
            alpha = find_min_new(x, d, 0, 1, count, z);
            
            prev_f = new_function(x, count, z);
            delta_x = vector_alpha(d, alpha);
            
            x = vector_add(x, delta_x);
            
            ++count;
            gradient = func_grad(x, count, z);
            
            //delta_x = vector_alpha(delta_x, -1.0);

            delta_d = vector_minus(gradient, old_grad);
            
            m_v = matrix_vector(tmp, delta_d);
            
            vec_mult = vector_multipl(delta_x, delta_x);
            vec_scalar = vector_scalar(delta_x, delta_d);
            dev_S = divide_matrix(vec_mult, vec_scalar);
            
            S = sum_matrix(tmp, dev_S);
            
            vec_mult = vector_multipl(m_v, m_v);
            vec_scalar = vector_scalar(m_v, delta_d);
            dev_S = divide_matrix(vec_mult, vec_scalar);
            
            S = sub_matrix(S, dev_S);
            cur_f = new_function(x, count, z);
           
            norm = sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]);
            
            
        }while(fabs(cur_f - prev_f) > eps && count < 1E6 && norm >= eps);
        
        
        values = x;
        x.clear();
        
        
    } while (fabs(my_function(X0) - my_function(values)) > eps);
        
        printf("%.2f  |  ",z);
   
    
    //cout<< "Значение точки минимума:\n";
    printf("%f  |  %f  |  %f  |  %f  |  %f  |  ", values[0], values[1], values[2], values[3], values[4]);
    
    //cout<< "Значение функции в точке минимума:\n";
    cout<<my_function(values)<<endl;
    z -= 0.01;
      
    }while(z >= 223.0);
    return 0;
}
