#pragma once
#include "matrices.hpp"
#include <cmath>
#include <opencv2/opencv.hpp>
#include <immintrin.h>

typedef struct conv_param
{
    size_t pad;
    size_t stride;
    size_t kernel_size;
    size_t in_channels;
    size_t out_channels;
    float *p_weight;
    float *p_bias;
} conv_param;

typedef struct fc_param
{
    size_t in_features;
    size_t out_features;
    float *p_weight;
    float *p_bias;
} fc_param;

template <class T>
void Convert(cv::Mat &A, matrices<T> &mat)
{
    int size = 128 * 128 * 3;
    for (int j = 2; j >= 0; --j)
    {
        int j_index_mat = 128 * 128 * j;
        for (int i = 2 - j, i_index = 0; i < size; i += 3, i_index++)
        {
            mat.data[j_index_mat + i_index] = A.data[i];
        }
    }
    for (int i = 0; i < 3 * 128 * 128; i++)
    {
        mat.data[i] = mat.data[i] / 255.f;
    }
}

template <class T>
matrices<T> ConvBNReLU(conv_param param, matrices<T> &mat)
{
    // check the conv_param
    if (param.stride == 0)
    {
        fprintf(stderr, "conv_param_stride error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (param.kernel_size == 0)
    {
        fprintf(stderr, "conv_param_kernel_size error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (param.kernel_size > mat.getRow() || param.kernel_size > mat.getCol())
    {
        fprintf(stderr, "conv_param_kernel_size error:out of range! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (param.in_channels == 0)
    {
        fprintf(stderr, "conv_param_in_channels error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (param.out_channels == 0)
    {
        fprintf(stderr, "conv_param_out_channels error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (param.p_weight == NULL)
    {
        fprintf(stderr, "conv_param_p_weight error:NULL! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (param.p_bias == NULL)
    {
        fprintf(stderr, "conv_param_p_bias error:NULL! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    // check the param of mat
    if (mat.data == NULL)
    {
        fprintf(stderr, "mat error:mat is empty! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getChannel() == 0)
    {
        fprintf(stderr, "mat_channel error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getRow() == 0)
    {
        fprintf(stderr, "mat_row error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getCol() == 0)
    {
        fprintf(stderr, "mat_col error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    // some parameters
    int channel_in = param.in_channels;
    int channel_out = param.out_channels;
    int kernel_size = param.kernel_size;
    int pad = param.pad;
    int stride = param.stride;
    int row = mat.getRow();
    int col = mat.getCol();
    int row_out = (row - kernel_size + 2 * pad) / stride + 1;
    int col_out = (col - kernel_size + 2 * pad) / stride + 1;
    matrices<T> result(row_out, col_out, channel_out);
    // the convolution part
    // version2
    matrices<T> kernel(channel_out, channel_in * kernel_size * kernel_size, 1);
    kernel.data = param.p_weight;
    ++kernel.refcount;
    int H_data = channel_in * kernel_size * kernel_size;
    int W_data = row_out * col_out;
    matrices<T> matrix(H_data, W_data, 1);

    for (int i = -pad, i_index = 0; i <= col - kernel_size + pad; i += stride, i_index++)
    {
        for (int j = -pad, j_index = 0; j <= row - kernel_size + pad; j += stride, j_index++)
        {
            int i_index_matrix = i_index * col_out + j_index;
            for (int k = 0; k < channel_in; k++)
            {
                int k_index_matrix = k * kernel_size * kernel_size;
                int k_index_mat = k * row * col;
                for (int u = 0; u < kernel_size; u++)
                {
                    int u_index_matrix = u * kernel_size;
                    int u_index_mat = (i + u) * col;
                    for (int v = 0; v < kernel_size; v++)
                    {
                        if (i + u < 0 || i + u >= row || j + v < 0 || j + v >= row)
                        {
                        }
                        else
                        {
                            matrix.data[(k_index_matrix + u_index_matrix + v) * W_data + i_index_matrix] = mat.data[k_index_mat + u_index_mat + (j + v)];
                        }
                    }
                }
            }
        }
    }

    result = kernel * matrix;
    result.setRow(row_out);
    result.setCol(col_out);
    result.setChannel(channel_out);

    // //version1
    //     for (int k = 0; k < channel_in; k++)
    //     {
    //         for (int i = -pad, i_index = 0; i <= col - kernel_size + pad; i += stride, i_index++)
    //         {
    //             for (int j = -pad, j_index = 0; j <= row - kernel_size + pad; j += stride, j_index++)
    //             {
    //                 //产生一个（i,j,k）对应的矩阵所对应的16个值
    //                 for (int chan = 0; chan < channel_out; chan++)
    //                 {
    //                     //产生一个（i,j,k）对应的矩阵与所对应一个out_channel的值
    //                     for (int u = 0; u < kernel_size; u++)
    //                     {
    //                         for (int v = 0; v < kernel_size; v++)
    //                         {
    //                             //相对应的矩阵进行点乘
    //                             if (i + u < 0 || i + u >= row || j + v < 0 || j + v >= row)
    //                             {
    //                             }
    //                             else
    //                             {
    //                                 result.data[chan * row_out * col_out + i_index * col_out + j_index] += mat.data[k * col * row + (i + u) * col + (j + v)] * param.p_weight[chan * channel_in * kernel_size * kernel_size + k * kernel_size * kernel_size + u * kernel_size + v];
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }

    for (int k = 0; k < channel_out; k++)
    {
        float bias = param.p_bias[k];
        int k_index_bais = k * col_out * row_out;
        for (int i = 0; i < col_out * row_out; i++)
        {
            result.data[k_index_bais + i] += bias;
        }
    }
    for (int i = 0; i < channel_out * row_out * col_out; i++)
    {
        result.data[i] = std::max(0.f, result.data[i]);
    }
    return result;
}

template <class T>
matrices<T> MaxPool(int downsampling, matrices<T> &mat)
{
    // check the arguments
    if (downsampling < 0 || downsampling == 0)
    {
        fprintf(stderr, "downsampling error:negative or zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getCol() % downsampling != 0 || mat.getRow() % downsampling != 0)
    {
        fprintf(stderr, "downsampling error:size do not fit! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.data == NULL)
    {
        fprintf(stderr, "mat_data error:NULL! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getChannel() == 0)
    {
        fprintf(stderr, "mat_channel error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getRow() == 0)
    {
        fprintf(stderr, "mat_row error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getCol() == 0)
    {
        fprintf(stderr, "mat_col error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    // some parameters
    int channel = mat.getChannel();
    int row = mat.getRow();
    int col = mat.getCol();
    int row_out = row / downsampling;
    int col_out = col / downsampling;
    matrices<T> result(row_out, col_out, channel);
    matrices<T> temp(row_out, col, channel);

    for (int k = 0; k < channel; ++k)
    {
        int k_index_temp = k * row_out * col;
        int k_index_mat = k * row * col;
        for (int i = 0, i_index = 0; i < row; i += downsampling, ++i_index)
        {
            int i_index_temp = i_index * col;
            int i_index_mat = i * col;
            for (int j = 0; j < col; j += downsampling)
            {
                temp.data[k_index_temp + i_index_temp + j] = std::max(mat.data[k_index_mat + i_index_mat + j], mat.data[k_index_mat + i_index_mat + j + 1]);
            }
        }
    }

    for (int k = 0; k < channel; ++k)
    {
        int k_index_temp = k * row_out * col;
        int k_index_mat = k * row * col;
        for (int i = 1, i_index = 0; i < row; i += downsampling, ++i_index)
        {
            int i_index_temp = i_index * col;
            int i_index_mat = i * col;
            for (int j = 0; j < col; j += downsampling)
            {
                temp.data[k_index_temp + i_index_temp + j + 1] = std::max(mat.data[k_index_mat + i_index_mat + j], mat.data[k_index_mat + i_index_mat + j + 1]);
            }
        }
    }

    for (int k = 0; k < channel; ++k)
    {
        int k_index_result = k * row_out * col_out;
        int k_index_temp = k * row_out * col;
        for (int i = 0; i < row_out; ++i)
        {
            int i_index_result = i * col_out;
            int i_index_temp = i * col;
            for (int j = 0, j_index = 0; j < col; j += downsampling, ++j_index)
            {
                result.data[k_index_result + i_index_result + j_index] = std::max(temp.data[k_index_temp + i_index_temp + j], temp.data[k_index_temp + i_index_temp + j + 1]);
            }
        }
    }
    return result;
}

template <class T>
matrices<T> FullConnectSoftMax(fc_param param, matrices<T> &mat)
{
    // check the paramter!
    if (param.in_features == 0)
    {
        fprintf(stderr, "fc_param_in_feature error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (param.out_features == 0)
    {
        fprintf(stderr, "fc_param_out_feature error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.data == NULL)
    {
        fprintf(stderr, "mat_data error:NULL! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getChannel() == 0)
    {
        fprintf(stderr, "mat_channel error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getRow() == 0)
    {
        fprintf(stderr, "mat_row error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }
    if (mat.getCol() == 0)
    {
        fprintf(stderr, "mat_col error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    }


    matrices<T> result(param.out_features, 1, 1);
    float expsum = 0;
    int out_features = param.out_features;
    int in_features = param.in_features;
    for (int i = 0; i < out_features; i++)
    {
        int i_index = i * in_features;
        float temp = 0;
        for (int j = 0; j < in_features; j++)
        {
            temp += param.p_weight[i_index + j] * mat.data[j];
        }
        result.data[i] = temp + param.p_bias[i];
        expsum += exp(result.data[i]);
    }

    for (int i = 0; i < param.out_features; i++)
    {
        result.data[i] = exp(result.data[i]) / expsum;
    }
    return result;
}
