#pragma once
#include<iostream>
using namespace std;

template <class T>
class matrices
{
private:
    size_t row;
    size_t col;
    size_t channel;
public:
    T* data;
    int* refcount;
    matrices();
    matrices(size_t row, size_t col, size_t channel);
    matrices(const matrices<T> & mat);
    ~matrices();
    matrices<T> operator=(const matrices<T> &mat);
    T &operator()(size_t i, size_t j,size_t channel);
    bool operator==(const matrices<T> &mat);
    bool operator!=(const matrices<T> &mat);
    matrices<T> operator+(const matrices<T> &mat);
    matrices<T> operator-(const matrices<T> &mat);
    matrices<T> operator+=(const matrices<T> &mat);
    matrices<T> operator-=(const matrices<T> &mat);
    matrices<T> operator*(matrices<T> &mat);
    matrices<T> operator*(T scalar);
    matrices<T> operator/(const matrices<T> &mat);
    matrices<T> operator/(T scalar);
    void at(size_t channel,size_t row_index,size_t col_index,T value);
    int getRow() const;
    int getCol() const;
    int getChannel() const;  
    void setCol(size_t col);
    void setRow(size_t row);
    void setChannel(size_t channel);
    friend std::ostream & operator <<(std::ostream & os, const matrices<T> mat){
        os <<"mat infos( col: " <<mat.col << " row: " << mat.row <<" channel: " << mat.channel << " )matrix :"<< endl;
        for(int i = 0;i < mat.channel*mat.col*mat.row;i++){
            os<< mat.data[i] <<endl;
        }
        os <<endl;
        return os; 
    }
};


template <class T>
matrices<T> ::matrices():row(1),col(1),channel(1)
{
data = new T[1]{0};
refcount = new int[1];
refcount[0] = 1;
} 

template <class T>
matrices<T> ::matrices(size_t row, size_t col, size_t channel):row(row),col(col),channel(channel)
{
if(row == 0){
   fprintf(stderr, "mat_row error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
}
if(col == 0){
   fprintf(stderr, "mat_col error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
}
if(channel == 0){
   fprintf(stderr, "mat_channel error:zero! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
}
data = new T[row*col*channel]{0};
refcount = new int[1];
refcount[0] = 1;
} 

template <class T>
matrices<T> ::matrices(const matrices<T> &c)
{
row = c.row;
col = c.col;
channel = c.channel;
data = c.data;
refcount = c.refcount; 
refcount[0] = refcount[0] +1;
} 

template <class T>
matrices<T> :: ~matrices()
{
if(refcount[0] == 1){
    delete[] data;
    delete[] refcount;
}
else{
    refcount[0] = refcount[0] -1;
}
data = NULL;
refcount = NULL;
} 

template <class T>
matrices<T>  matrices<T>::operator=(const matrices<T> &c)
{
if(refcount[0] == 1){
    delete[] data;
    delete[] refcount;
}
else{
    refcount[0] = refcount[0] -1;
}
this->col =c.col;
this->row = c.row;
this->channel = c.channel;
this->data = c.data;
this->refcount = c.refcount;
refcount[0] = refcount[0]+1;
return *this;
} 

template <class T>
matrices<T> matrices<T>::operator+(const matrices<T> &mat){
if(this->col == mat.col && this->row == mat.row && this->channel == mat.channel){
matrices<T> temp1(row,col,channel);
for(int i = 0; i< row*col*channel;i++){
  temp1.data[i] = this->data[i] + mat.data[i];
  
}
return temp1;
}
else{
    matrices<T> temp2;
    fprintf(stderr, "mat error:size do not fit! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    return temp2;
}
}

template <class T>
matrices<T> matrices<T>::operator-(const matrices<T> &mat){
if(this->col == mat.col && this->row == mat.row && this->channel == mat.channel){
matrices<T> temp1(row,col,channel);
for(int i = 0; i< row*col*channel;i++){
  temp1.data[i] = this->data[i] - mat.data[i];
}
return temp1;
}
else{
    matrices<T> temp2;
    fprintf(stderr, "mat error:size do not fit! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
    return temp2;
}
}


template <class T>
matrices<T> matrices<T>::operator+=(const matrices<T> &mat){
if(this->col == mat.col && this->row == mat.row && this->channel == mat.channel){
for(int i = 0; i< row*col*channel;i++){
  this->data[i] = this->data[i] +mat.data[i];
}
}
else{
fprintf(stderr, "mat error:size do not fit! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
}
return *this;
}

template <class T>
matrices<T> matrices<T>::operator-=(const matrices<T> &mat){
if(this->col == mat.col && this->row == mat.row && this->channel == mat.channel){
for(int i = 0; i< row*col*channel;i++){
  this->data[i] = this->data[i] -mat.data[i];
}

}
else{
fprintf(stderr, "mat error:size do not fit! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
}
return *this;
}

template <class T>
bool matrices<T>::operator==(const matrices<T> &mat){
bool temp = true;
if(row != mat.row){
    temp = false;
}
if(col != mat.col){
    temp = false;
}
if(channel != mat.channel){
    temp = false;
}
for(int i = 0;i < row*col*channel;i++){
    if(data[i] != mat.data[i]){
        temp = false;
        break;
    }
}
return temp;
}

template <class T>
bool matrices<T>:: operator!=(const matrices<T> &mat){
bool temp = (*this == mat);
return !temp;
}

template <class T>
matrices<T> matrices<T>::operator*(T scalar){
matrices<T> temp(row,col,channel);
for(int i = 0;i < row*col*channel;i++){
    temp.data[i] = this->data[i] * scalar;
}
return temp;
}

template <class T>
matrices<T> matrices<T>::operator*(matrices<T> &mat){
if(col != mat.row || channel != mat.channel){
matrices<T> temp;
fprintf(stderr, "mat error:size do not fit! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
return temp;
}
else{
matrices<T> temp(row,mat.col,channel);
 int m = row;
 int n = mat.col;
 int p = col;
 int size1 = m*n;
 int size2 = m*p;
 int size3 = n*p;
        for(int chan = 0; chan < channel;chan++){
                for (int i = 0; i < m; i++)
                {     
                    for (int k = 0; k < p; k++)
                    {
                        T r = data[chan*size2 + i*p + k];
                        for (int j = 0; j < n; j++)
                        {
                            temp.data[chan*size1 + i*n + j] += r * mat.data[chan*size3 + k*n + j];
                        }
                    }
                }
        }
return temp;
}
}

template <class T>
void matrices <T>::at(size_t channel,size_t row_index,size_t col_index,T value){
this->data[(channel-1)*col*row+row_index*col+col_index] = value;
}

template <class T>
matrices<T> matrices<T>::operator/(const matrices<T> &mat){
if(row != mat.row || col != mat.col){
matrices<T> temp;
fprintf(stderr, "mat error:size do not fit! %s(%d)-%s\n", __FILE__, __LINE__, __FUNCTION__);
return temp;
}
else{
matrices<T> temp(row,col,channel);
for(int i = 0;i< row*col*channel;i++){
    temp.data[i] = this->data[i]/mat.data[i]; 
}
return temp;
}
}

template <class T>
matrices<T> matrices<T>::operator/(T scalar){
matrices<T> temp(row,col,channel);
for(int i = 0;i < col*row*channel;i++){
    temp.data[i] = this->data[i]/scalar;
}
return temp;
}

template <class T>
int matrices<T>::getRow() const{
return row;
}

template <class T>
int matrices<T>::getCol() const{
return col;
}

template <class T>
int matrices<T>::getChannel() const{
return channel;
}

template <class T>
T &matrices<T>::operator()(size_t i, size_t j, size_t channel)
{
    return data[(channel - 1) *row*col + i*col + j];
}

template <class T>
void matrices<T>::setCol(size_t col){
    this->col = col;
}

template <class T>
void matrices<T>::setRow(size_t row){
    this->row = row;
}

template <class T>
void matrices<T>::setChannel(size_t channel){
    this->channel = channel;
}
