// Markus Buchholz, 2023

// g++ attractor_halvorsen.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <vector>
#include <tuple>
#include <math.h>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//--------------------------------------------------------------------------------

float a = 1.4;

float dt = 0.01;

//--------------------------------------------------------------------------------
// dot x
float function1(float x, float y, float z)
{
    return -a * x - 4 * y - 4 * z - y * y;
}

//--------------------------------------------------------------------------------
// dot y

float function2(float x, float y, float z)
{
    return -a * y - 4 * z - 4 * x - z * z;
}

//--------------------------------------------------------------------------------
// dot z

float function3(float x, float y, float z)
{
    return -a * z - 4 * x - 4 * y - x * x;
}

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> methodRungeKutta1Diff()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;
    std::vector<float> diffEq3;

    std::vector<float> time;

    // init values
    float x1 = 0.1; // x1
    float x2 = 0.0; // x2
    float x3 = 0.0; // x3
    float t = 0.0;  // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    diffEq3.push_back(x3);
    time.push_back(t);

    for (int ii = 0; ii < 6000; ii++)
    {
        t = t + dt;
        float k11 = function1(x1, x2, x3);
        float k12 = function2(x1, x2, x3);
        float k13 = function3(x1, x2, x3);

        float k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13);
        float k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13);
        float k23 = function3(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13);

        float k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23);
        float k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23);
        float k33 = function3(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23);

        float k41 = function1(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33);
        float k42 = function2(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33);
        float k43 = function3(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);
        x3 = x3 + dt / 6.0 * (k13 + 2 * k23 + 2 * k33 + k43);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        diffEq3.push_back(x3);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, diffEq3, time);
}

//---------------------------------------------------------------------------------------------------------

void plot2D(std::vector<float> xX, std::vector<float> yY)
{
    plt::title("Halvorsen attractor ");
    plt::named_plot("solution", xX, yY);
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::legend();
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::show();
}

//---------------------------------------------------------------------------------------------------------

void plot3D(std::vector<float> xX, std::vector<float> yY, std::vector<float> zZ)
{

    plt::plot3(xX, yY, zZ);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("z");
    plt::show();
}

//---------------------------------------------------------------------------------------------------------

int main()
{
    auto halvorsen = methodRungeKutta1Diff();
    auto xX = std::get<0>(halvorsen);
    auto yY = std::get<1>(halvorsen);
    auto zZ = std::get<2>(halvorsen);
    plot3D(xX, yY, zZ);
    plot2D(xX, yY);
}