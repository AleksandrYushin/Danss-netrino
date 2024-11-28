#include <cmath>
//#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <string>

//Функтор возращающий растояние между точкой детектора (x,y,z) и точкой реактора (r, phi, h)
/*Модель: реактор - цилиндр с высотой H и радиусом R (координаты точки - цилиндрические от центра основания)
          детектор - куб с длиной ребра a (координаты точки - декартовы от центра основания)
Взаимное расположение: куб и цилиндр стоят на одной плоскости. L - растояие между центрами фигур, theta - угол поворота куба вокруг своей оси (т.е. угл между прямой соед. центры фигур и нормалью грани)*/
class distance{
    public:
    distance(double L, double theta): L(L), theta(theta) {};

    void setL(double a){L = a;};
    void setTheta(double a){theta = a;};

    double getL(){return L;};
    double getTheta(){return theta;};

    double operator()(double x, double y, double z, double r, double phi, double h){
        return std::sqrt(std::pow(L-r*std::cos(phi)-y*std::cos(theta) + x*std::sin(theta), 2)+ std::pow(r*std::sin(phi)+x*std::cos(theta)+y*std::sin(theta), 2)+ std::pow(z-h, 2));
    };
    private:
    double L;
    double theta;
};


// Класс гистограммы. x_max - верхяя граница интервала значений (x_min = 0 по умолчанию), N - число бинов.
class histogram{
    public:
    histogram(double x_max, double x_min, int N): x_max(x_max), x_min(x_min), N(N) {data = new int[N];};
    histogram(double x_max, double x_min, double err): x_max(x_max), x_min(x_min) {
        N = int((x_max-x_min)/err); //err - требуемая max ошибка округления при заполнении гистограммы. Из него получаем N.
        data = new int[N];
    };
    ~histogram(){
        delete[] data;
    };

    void setX_max(double a){x_max = a;};
    void setX_min(double a){x_min = a;};
    void setN(int a){N = a;};

    double getX_max(){return x_max;};
    double getX_min(){return x_min;};
    int getN(){return N;};
    int* gatData() {return data;};

    void push_point(double x){
        if (x < x_max && x > x_min){
            int n = int((x-x_min)/(x_max-x_min) *N);
            data[n]++;
        };
    };

    double operator()(double x){
        if (x < x_max && x > x_min){
            int n = int((x-x_min)/(x_max-x_min) *N);
            return data[n];
        }
        return 0;
    };

    private:
    double x_max;
    double x_min;
    int N;
    int *data;
};

void BruteForceMethod(double delita, distance* dist, histogram* F, double H, double R, double a){
    //В тупую считаем распределение перебирая все точки
    int H_max = int(H/delita);
    int R_max = int(R/delita);
    int a_max = int(a/delita);

    for (int i_z = 0; i_z<= a_max; i_z++){
        for (int i_y = -a_max/2; i_y<= a_max/2; i_y++){
            for (int i_x = -a_max/2; i_x <= a_max/2; i_x++){
                for (int i_r = 0; i_r<= R_max*R_max; i_r++){
                    for (int i_phi = 0; i_phi< int(2*M_PI/delita); i_phi++){
                        for (int i_h = 0; i_h<= H_max; i_h++){
                            F->push_point(dist->operator()(i_x, i_y, i_z, std::sqrt(i_r), i_phi, i_h));    
                        };
                    };
                };
            };
        };
    };
};

void MonteKarloMethod(int N, distance* dist, histogram* F, double H, double R, double a){
    //Метод Монте-Карло - перебираем только часть точек
    std::uniform_real_distribution<double> unif(0,1);   //ГПСЧ
    std::default_random_engine re;

    for (int i = 0; i< N; i++){
        double r_x = unif(re)*a - a/2;
        double r_y = unif(re)*a - a/2;
        double r_z = unif(re)*a;
        double r_r = std::sqrt(unif(re)*R*R);
        double r_phi = unif(re)*2*M_PI;
        double r_h = unif(re)*H;

        F->push_point(dist->operator()(r_x, r_y, r_z, r_r, r_phi, r_h));  
    };
};

int main(){
    //Создаём наш функтор и распределение
    distance dist = distance(10.0, 0);
    histogram F = histogram(15.0, 5.0, 0.01);  

    double H = 2;
    double R = 2;
    double a = 1;

    // BruteForceMethod(0.001, &dist, &F, H, R, a); 
    MonteKarloMethod(100000000, &dist, &F, H, R, a);    


    //Вывод результата
    std::ofstream f;
    f.open("Data_end"); //std::to_string(i)
    int * data = F.gatData();
    f << "#parameters:" << std::endl;
    f << "#H = " << H << std::endl;
    f << "#R = " << R << std::endl;
    f << "#a = " << a << std::endl;
    f << "#L = " << dist.getL() << std::endl;
    f << "#theta = " << dist.getTheta() << std::endl;
    f << "#x max/min = " << F.getX_max() << " " << F.getX_min() << std::endl;
    double err = (F.getX_max() - F.getX_min() )/F.getN();
    f << "#err = " << err << std::endl;
    f << "#data:" << std::endl;
    for (int i =0 ; i< F.getN(); i++){
        f << data[i] << std::endl;
    };
    f.close();

    return 1;
};