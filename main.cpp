#include <cmath>
//#include <vector>
#include <iostream>
#include <random>

//Функтор возращающий растояние между точкой детектора (x,y,z) и точкой реактора (r, fi, h)
class distance{
    public:
    distance(double L, double theta): L(L), theta(theta) {};

    void setL(double a){L = a;};
    void setTheta(double a){theta = a;};

    double getL(){return L;};
    double getTheta(){return theta;};

    double operator()(double x, double y, double z, double r, double fi, double h){
        return std::sqrt(std::pow(L-r*std::cos(fi)-y*std::cos(theta) + x*std::sin(theta), 2)+ std::pow(r*std::sin(fi)+x*std::cos(theta)+y*std::sin(theta), 2)+ std::pow(z-h, 2));
    };
    private:
    double L;
    double theta;
};

class histogram{
    public:
    histogram(double x_max, int N): x_max(x_max), N(N) {data = new int[N];};
    histogram(double x_max, double err): x_max(x_max) {
        N = int(x_max/err);
        data = new int[N];
    };

    void setX_max(double a){x_max = a;};
    void setN(int a){N = a;};

    double getX_max(){return x_max;};
    int getN(){return N;};
    int* gatData() {return data;};

    void push_point(double x){
        if (x < x_max){
            int n = int(x/x_max *N);
            data[n]++;
        };
    };

    double operator()(double x){
        int n = int(x/x_max *N);
        return data[n];
    };

    private:
    double x_max;
    int N;
    int *data;
};

int main(){
    //Создаём наш функтор и распределение
    distance dist = distance(100.0, 0.0);
    histogram F = histogram(150.0, 150);  

    double H = 10;
    double R = 10;
    double a = 20;
    
    //В тупую считаем распределение перебирая все точки
    double err = 0.1;   //шаг сетки
    int H_max = int(H/err);
    int R_max = int(R/err);
    int a_max = int(a/err);

    for (int i_z = 0; i_z<= a_max; i_z++){
        for (int i_y = -a_max/2; i_y<= a_max/2; i_y++){
            for (int i_x = -a_max/2; i_x <= a_max/2; i_x++){
                for (int i_r = 0; i_r<= R_max; i_r++){
                    for (int i_fi = 0; i_fi< int(2*M_PI/err); i_fi++){
                        for (int i_h = 0; i_h<= H_max; i_h++){
                            F.push_point(dist(i_x, i_y, i_z, i_r, i_fi, i_h));    
                        };
                    };
                };
            };
        };
    };

    //Метод Монте-Карло - перебираем только часть точек
    
    //Создаём ГПЧ
    std::mt19937 engine;
    std::random_device device;
    engine.seed(device());

    

    //Вывод результата
    int * data = F.gatData();
    for (int i =0 ; i<= F.getN(); i++){
        std::cout << data[i] << std::endl;
    };

    return 1;
};