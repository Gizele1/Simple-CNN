#include <iostream>
#include <sys/time.h>
#include "header.h"
#pragma GCC optimize(3, "Ofast", "inline")
using namespace std;

int main(int argc, char *argv[])
{
    struct timeval t_start, t_end, t1_s, t1_e, t2_s, t2_e, t3_s, t3_e, t4_s, t4_e, t5_s, t5_e, t6_s, t6_e;
    cv::Mat A = cv::imread(argv[1]);
    extern conv_param conv_params[3];
    extern fc_param fc_params[1];
    matrices<float> mat0(128, 128, 3);
    Convert(A, mat0);

    gettimeofday(&t_start, NULL);

    gettimeofday(&t1_s, NULL);
    matrices<float> mat1 = ConvBNReLU(conv_params[0], mat0);
    gettimeofday(&t1_e, NULL);
    double deltaT1 = (t1_e.tv_sec - t1_s.tv_sec) * 1000000 + t1_e.tv_usec - t1_s.tv_usec;
    cout << "Fist ConvBNReLU = " << deltaT1 << "us" << endl;

    gettimeofday(&t2_s, NULL);
    matrices<float> mat2 = MaxPool(2, mat1);
    gettimeofday(&t2_e, NULL);
    double deltaT2 = (t2_e.tv_sec - t2_s.tv_sec) * 1000000 + t2_e.tv_usec - t2_s.tv_usec;
    cout << "Fist MaxPool = " << deltaT2 << "us" << endl;

    gettimeofday(&t3_s, NULL);
    matrices<float> mat3 = ConvBNReLU(conv_params[1], mat2);
    gettimeofday(&t3_e, NULL);
    double deltaT3 = (t3_e.tv_sec - t3_s.tv_sec) * 1000000 + t3_e.tv_usec - t3_s.tv_usec;
    cout << "Second ConvBNReLU = " << deltaT3 << "us" << endl;

    gettimeofday(&t4_s, NULL);
    matrices<float> mat4 = MaxPool(2, mat3);
    gettimeofday(&t4_e, NULL);
    double deltaT4 = (t4_e.tv_sec - t4_s.tv_sec) * 1000000 + t4_e.tv_usec - t4_s.tv_usec;
    cout << "Second MaxPool = " << deltaT4 << "us" << endl;

    gettimeofday(&t5_s, NULL);
    matrices<float> mat5 = ConvBNReLU(conv_params[2], mat4);
    gettimeofday(&t5_e, NULL);
    double deltaT5 = (t5_e.tv_sec - t5_s.tv_sec) * 1000000 + t5_e.tv_usec - t5_s.tv_usec;
    cout << "Third ConvBNReLU= " << deltaT5 << "us" << endl;

    gettimeofday(&t6_s, NULL);
    matrices<float> mat6 = FullConnectSoftMax(fc_params[0], mat5);
    gettimeofday(&t6_e, NULL);
    double deltaT6 = (t6_e.tv_sec - t6_s.tv_sec) * 1000000 + t6_e.tv_usec - t6_s.tv_usec;
    cout << "FullConnectSoftMax= " << deltaT6 << "us" << endl;

    gettimeofday(&t_end, NULL);
    double deltaT = (t_end.tv_sec - t_start.tv_sec) * 1000000 + t_end.tv_usec - t_start.tv_usec;
    cout << "all the time = " << deltaT << "us" << endl;
    cout << mat6 << endl;
}