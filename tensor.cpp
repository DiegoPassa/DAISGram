<<<<<<< HEAD
#include "tensor.h"

#include <math.h>

#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "dais_exc.h"
=======
#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"
>>>>>>> 5994379c2e9a2fef0f3b6192386a473179f024f1

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

<<<<<<< HEAD
/**
 * struct implementation usato per nascondere l'implementazione dei tensori
 * 
 */
struct Tensor::Impl {
    float* data;
    float** cols_p;
    float*** matrix_p;
};

/**
 * Class constructor
 * 
 * Parameter-less class constructor 
 */
Tensor::Tensor() {
    pimpl = nullptr;
    this->col = 0;
    this->row = 0;
    this->dep = 0;
}

/**
 * metodo privato usato per allocare il tensore così da evitare ridondanza del codice
 * l'allocazione prevede che la matrice sia flattened così da permettere una maggiore località spaziale
 * e in'oltre l'accesso avviene tramite due altri array di pointer che puntano alla matrice permettendo di accedervi tramite 
 * parentesi quadre e senza dover ogni volta utilizzare la formula di accesso a matrice flattened
 * 
 * @param row 
 * @param col 
 * @param dep 
 */
void Tensor::allocate_matrix(int row, int col, int dep) {
    if (row == 0 || col == 0 || dep == 0)
        throw(unknown_exception());

    //alloco l'array contenente i dati del tensore
    pimpl->data = new float[row * col * dep];

    //alloco i due array utilizzati per accedere al tensore
    pimpl->cols_p = new float*[row * col];
    pimpl->matrix_p = new float**[row];

    //sistemo gli array per far si che si possa dereferenziare matrix_p con le parentesi quadre
    for (int i = 0; i < row; i++) {
        pimpl->matrix_p[i] = &(pimpl->cols_p[i * col]);

        for (int j = 0; j < col; j++)
            pimpl->matrix_p[i][j] = pimpl->data + i * col * dep + j * dep;
    }
}

/**
 * Copy constructor
 * 
 * This constructor copies the data from another Tensor
 *      
 * @return the new Tensor
 */
Tensor::Tensor(const Tensor& that) {
    pimpl = new Impl;

    row = that.row;
    col = that.col;
    dep = that.dep;

    allocate_matrix(row, col, dep);

    for (int i = 0; i < row * col * dep; i++) {
        pimpl->data[i] = that.pimpl->data[i];
    }
}

/**
 * Class constructor
 * 
 * Creates a new tensor of size r*c*d initialized at value v
 * 
 * @param r
 * @param c
 * @param d
 * @param v
 * @return new Tensor
 */
Tensor::Tensor(int r, int c, int d, float v) {
    init(r, c, d, v);
}

/**
 * Class distructor
 * 
 * Cleanup the data when deallocated
 */
Tensor::~Tensor() {
    if (pimpl) {
        delete[] pimpl->data;
        delete[] pimpl->cols_p;
        delete[] pimpl->matrix_p;
        delete pimpl;
        pimpl = nullptr;
    }
}

/**
 * Constant Initialization
 * 
 * Perform the initialization of the tensor to a value v
 * 
 * @param r The number of rows
 * @param c The number of columns
 * @param d The depth
 * @param v The initialization value
 */
void Tensor::init(int r, int c, int d, float v) {
    if (r < 0 || c < 0 || d < 0)
        throw(unknown_exception());

    pimpl = new Impl;
    allocate_matrix(r, c, d);

    row = r;
    col = c;
    dep = d;

    for (int i = 0; i < row * col * dep; i++) {
        pimpl->data[i] = v;
    }
}

/** 
 * Get minimum 
 * 
 * Compute the minimum value considering a particular index in the third dimension
 * 
 * @return the minimum of data( , , k)
 */
float Tensor::getMin(int k) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    float min = FLT_MAX;

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (pimpl->matrix_p[i][j][k] < min) min = pimpl->matrix_p[i][j][k];
        }
    }

    return min;
}

/** 
 * Get maximum 
 * 
 * Compute the maximum value considering a particular index in the third dimension
 * 
 * @return the maximum of data( , , k)
 */
float Tensor::getMax(int k) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    float max = FLT_MIN;

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (pimpl->matrix_p[i][j][k] > max) max = pimpl->matrix_p[i][j][k];
        }
    }

    return max;
}

/**
 * Operator overloading ==
 * 
 * It performs the point-wise equality check between two Tensors.
 * 
 * The equality check between floating points cannot be simply performed using the 
 * operator == but it should take care on their approximation.
 * 
 * This approximation is known as rounding (do you remember "Architettura degli Elaboratori"?)
 *  
 * For example, given a=0.1232 and b=0.1233 they are 
 * - the same, if we consider a rounding with 1, 2 and 3 decimals 
 * - different when considering 4 decimal points. In this case b>a
 * 
 * So, given two floating point numbers "a" and "b", how can we check their equivalence? 
 * through this formula:
 * 
 * a ?= b if and only if |a-b|<EPSILON
 * 
 * where EPSILON is fixed constant (defined at the beginning of this header file)
 * 
 * Two tensors A and B are the same if:
 * A[i][j][k] == B[i][j][k] for all i,j,k 
 * where == is the above formula.
 * 
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 * 
 * @return returns true if all their entries are "floating" equal
 */
bool Tensor::operator==(const Tensor& rhs) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    if (row != rhs.row || col != rhs.col || dep != rhs.dep)
        throw(dimension_mismatch());

    bool equals = true;

    for (size_t i = 0; i < (size_t)row * col * dep && equals; i++) {
        if (fabs(pimpl->data[i] - rhs.pimpl->data[i]) >= EPSILON) {
            equals = false;
        }
    }

    return equals;
}

/**
 * Operator overloading ()
 * 
 * if indexes are out of bound throw index_out_of_bound() exception
 * 
 * @return the value at location [i][j][k]
 */
float Tensor::operator()(int i, int j, int k) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    if (i < 0 || i >= row || j < 0 || j >= col || k < 0 || k >= dep)
        throw(index_out_of_bound());

    return pimpl->matrix_p[i][j][k];
}

/**
 * Operator overloading ()
 * 
 * Return the pointer to the location [i][j][k] such that the operator (i,j,k) can be used to 
 * modify tensor data.
 * 
 * If indexes are out of bound throw index_out_of_bound() exception
 * 
 * @return the pointer to the location [i][j][k]
 */
float& Tensor::operator()(int i, int j, int k) {
    if (!pimpl)
        throw(tensor_not_initialized());

    if (i < 0 || i >= row || j < 0 || j >= col || k < 0 || k >= dep)
        throw(index_out_of_bound());

    return pimpl->matrix_p[i][j][k];
}

/**
 * Operator overloading <<
 * 
 * Use the overaloading of << to show the content of the tensor.
 * 
 * You are free to chose the output format, btw we suggest you to show the tensor by layer.
 * 
 * [..., ..., 0]
 * [..., ..., 1]
 * ...
 * [..., ..., k]
 */
ostream& operator<<(ostream& stream, const Tensor& obj) {
    if (!obj.pimpl)
        throw(tensor_not_initialized());

    for (int i = 0; i < obj.row; i++) {
        for (int j = 0; j < obj.col; j++) {
            stream << "[";
            for (int k = 0; k < obj.dep - 1; k++) {
                stream << obj.pimpl->matrix_p[i][j][k] << ",";
            }
            stream << obj.pimpl->matrix_p[i][j][obj.dep - 1] << "] ";
        }
        stream << "\n";
    }
    stream << "\n";

    return stream;
}

/**
 * Tensor Clamp
 * 
 * Clamp the tensor such that the lower value becomes low and the higher one become high.
 * 
 * @param low Lower value
 * @param high Higher value 
 */
void Tensor::clamp(float low, float high) {
    if (!pimpl)
        throw(tensor_not_initialized());

    if (low > high)
        throw(unknown_exception());

    for (int i = 0; i < row * col * dep; i++) {
        if (pimpl->data[i] < low)
            pimpl->data[i] = low;
        else if (pimpl->data[i] > high)
            pimpl->data[i] = high;
    }
}

/**
 * Tensor Rescaling
 * 
 * Rescale the value of the tensor following this rule:
 * 
 * newvalue(i,j,k) = ((data(i,j,k)-min(k))/(max(k)-min(k)))*new_max
 * 
 * where max(k) and min(k) are the maximum and minimum value in the k-th channel.
 * 
 * new_max is the new value for the maximum
 * 
 * @param new_max New maximum vale
 */
void Tensor::rescale(float new_max) {
    if (!pimpl)
        throw(tensor_not_initialized());

    for (int k = 0; k < dep; k++) {
        float min = getMin(k);
        float max = getMax(k);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                pimpl->matrix_p[i][j][k] = ((pimpl->matrix_p[i][j][k] - min) / (max - min)) * new_max;
            }
        }
    }
}

/**
 * Tensor padding
 * 
 * Zero pad a tensor in height and width, the new tensor will have the following dimensions:
 * 
 * (rows+2*pad_h) x (cols+2*pad_w) x (depth) 
 * 
 * @param pad_h the height padding
 * @param pad_w the width padding
 * @return the padded tensor
 */
Tensor Tensor::padding(int pad_h, int pad_w) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    Tensor new_t{row + pad_h * 2, col + pad_w * 2, dep};

    for (int k = 0; k < new_t.dep; k++) {
        for (int i = 0; i < new_t.row; i++) {
            for (int j = 0; j < new_t.col; j++) {
                if (i >= pad_h && i < new_t.row - pad_h && j >= pad_w && j < new_t.col - pad_w)
                    new_t(i, j, k) = pimpl->matrix_p[i - pad_h][j - pad_w][k];
                else
                    new_t(i, j, k) = 0;
            }
        }
    }

    return new_t;
}

/**
 * Subset a tensor
 * 
 * retuns a part of the tensor having the following indices:
 * row_start <= i < row_end  
 * col_start <= j < col_end 
 * depth_start <= k < depth_end
 * 
 * The right extrema is NOT included
 * 
 * @param row_start 
 * @param row_end 
 * @param col_start
 * @param col_end
 * @param depth_start
 * @param depth_end
 * @return the subset of the original tensor
 */
Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    if (row_start < 0 || (int)row_start > row + 1 || row_end < 0 || (int)row_end > row + 1 ||
        col_start < 0 || (int)col_start > col + 1 || col_end < 0 || (int)col_end > col + 1 ||
        depth_start < 0 || (int)depth_start > dep + 1 || depth_end < 0 || (int)depth_end > dep + 1)
        throw(index_out_of_bound());

    if (row_end <= row_start || col_end <= col_start || depth_end <= depth_start)
        throw(unknown_exception());

    Tensor new_t{(int)(row_end - row_start), (int)(col_end - col_start), (int)(depth_end - depth_start)};

    for (unsigned int k = depth_start; k < depth_end; k++) {
        for (unsigned int i = row_start; i < row_end; i++) {
            for (unsigned int j = col_start; j < col_end; j++) {
                new_t(i - row_start, j - col_start, k - depth_start) = pimpl->matrix_p[i][j][k];
            }
        }
    }

    return new_t;
}

/** 
 * Concatenate 
 * 
 * The function concatenates two tensors along a give axis
 * 
 * Example: this is of size 10x5x6 and rhs is of 25x5x6
 * 
 * if concat on axis 0 (row) the result will be a new Tensor of size 35x5x6
 * 
 * if concat on axis 1 (columns) the operation will fail because the number 
 * of rows are different (10 and 25).
 * 
 * In order to perform the concatenation is mandatory that all the dimensions 
 * different from the axis should be equal, other wise throw concat_wrong_dimension(). 
 *  
 * @param rhs The tensor to concatenate with
 * @param axis The axis along which perform the concatenation 
 * @return a new Tensor containing the result of the concatenation
 */
Tensor Tensor::concat(const Tensor& rhs, int axis) const {
    if (!pimpl || !rhs.pimpl)
        throw(tensor_not_initialized());

    Tensor new_t;

    switch (axis) {
        case 0:
            if (col != rhs.col || dep != rhs.dep)
                throw(concat_wrong_dimension());

            new_t.init(row + rhs.row, col, dep);
            for (int k = 0; k < dep; k++) {
                for (int i = 0; i < row + rhs.row; i++) {
                    for (int j = 0; j < col; j++) {
                        if (i < row)
                            new_t(i, j, k) = pimpl->matrix_p[i][j][k];
                        else
                            new_t(i, j, k) = rhs.pimpl->matrix_p[i - row][j][k];
                    }
                }
            }

            break;
        case 1:
            if (row != rhs.row || dep != rhs.dep)
                throw(concat_wrong_dimension());

            new_t.init(row, col + rhs.col, dep);
            for (int k = 0; k < dep; k++) {
                for (int i = 0; i < row; i++) {
                    for (int j = 0; j < col + rhs.col; j++) {
                        if (j < col)
                            new_t(i, j, k) = pimpl->matrix_p[i][j][k];
                        else
                            new_t(i, j, k) = rhs.pimpl->matrix_p[i][j - col][k];
                    }
                }
            }

            break;
        case 2:
            if (col != rhs.col || row != rhs.row)
                throw(concat_wrong_dimension());

            new_t.init(row, col, dep + rhs.dep);
            for (int k = 0; k < dep + rhs.dep; k++) {
                for (int i = 0; i < row; i++) {
                    for (int j = 0; j < col; j++) {
                        if (k < dep)
                            new_t(i, j, k) = pimpl->matrix_p[i][j][k];
                        else
                            new_t(i, j, k) = rhs.pimpl->matrix_p[i][j][k - dep];
                    }
                }
            }

            break;
        default:
            throw(unknown_exception());
    }

    return new_t;
}

/** 
 * Convolution 
 * 
 * This function performs the convolution of the Tensor with a filter.
 * 
 * The filter f must have odd dimensions and same depth. 
 * 
 * Remeber to apply the padding before running the convolution
 *  
 * @param f The filter
 * @return a new Tensor containing the result of the convolution
 */
Tensor Tensor::convolve(const Tensor& f) const {
    if (!pimpl || !f.pimpl)
        throw(tensor_not_initialized());

    if (f.col % 2 == 0 || f.row % 2 == 0) throw(filter_odd_dimensions());
    if (dep != f.dep) throw dimension_mismatch();

    int pad_h = (f.row - 1) / 2;
    int pad_w = (f.col - 1) / 2;

    Tensor conv = {row, col, dep};
    Tensor padded = padding(pad_h, pad_w);

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            for (int k = 0; k < dep; k++) {
                float result = 0;

                for (int i_filter = 0; i_filter < f.row; i_filter++) {
                    for (int j_filter = 0; j_filter < f.col; j_filter++) {
                        result += padded(i + i_filter, j + j_filter, k) * f(i_filter, j_filter, k);
                    }
                }

                conv(i, j, k) = result;
            }
        }
    }

    return conv;
}

/**
 * Operator overloading = (assignment) 
 * 
 * Perform the assignment between this object and another
 * 
 * @return a reference to the receiver object
 */
Tensor& Tensor::operator=(const Tensor& other) {
    if (this != &other) {
        if (pimpl) {
            delete[] pimpl->data;
            delete[] pimpl->cols_p;
            delete[] pimpl->matrix_p;
            delete pimpl;
            pimpl = nullptr;
        }

        if (!other.pimpl) {
            pimpl = nullptr;
            col = 0;
            row = 0;
            dep = 0;
        } else {
            init(other.row, other.col, other.dep);
            for (int i = 0; i < row * col * dep; i++) {
                pimpl->data[i] = other.pimpl->data[i];
            }
        }
    }

    return (*this);
}

/**
 * Operator overloading -
 * 
 * It performs the point-wise difference between two Tensors.
 * 
 * result(i,j,k)=this(i,j,k)-rhs(i,j,k)
 * 
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator-(const Tensor& rhs) const {
    if (!pimpl || !rhs.pimpl)
        throw(tensor_not_initialized());

    if (row != rhs.row || col != rhs.col || dep != rhs.dep) {
        throw(dimension_mismatch());
    }

    Tensor newTensor{*this};
    for (int i = 0; i < row * col * dep; i++) {
        newTensor.pimpl->data[i] = pimpl->data[i] - rhs.pimpl->data[i];
    }

    return newTensor;
}

/**
 * Operator overloading +
 * 
 * It performs the point-wise sum between two Tensors.
 * 
 * result(i,j,k)=this(i,j,k)+rhs(i,j,k)
 * 
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator+(const Tensor& rhs) const {
    if (!pimpl || !rhs.pimpl)
        throw(tensor_not_initialized());

    if (row != rhs.row || col != rhs.col || dep != rhs.dep) {
        throw(dimension_mismatch());
    }

    Tensor newTensor{*this};
    for (int i = 0; i < row * col * dep; i++) {
        newTensor.pimpl->data[i] = pimpl->data[i] + rhs.pimpl->data[i];
    }

    return newTensor;
}

/**
 * Operator overloading *
 * 
 * It performs the point-wise product between two Tensors.
 * 
 * result(i,j,k)=this(i,j,k)*rhs(i,j,k)
 * 
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator*(const Tensor& rhs) const {
    if (!pimpl || !rhs.pimpl)
        throw(tensor_not_initialized());

    if (row != rhs.row || col != rhs.col || dep != rhs.dep) {
        throw(dimension_mismatch());
    }

    Tensor newTensor{*this};
    for (int i = 0; i < row * col * dep; i++) {
        newTensor.pimpl->data[i] = pimpl->data[i] * rhs.pimpl->data[i];
    }

    return newTensor;
}

/**
 * Operator overloading /
 * 
 * It performs the point-wise division between two Tensors.
 * 
 * result(i,j,k)=this(i,j,k)/rhs(i,j,k)
 * 
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator/(const Tensor& rhs) const {
    if (!pimpl || !rhs.pimpl)
        throw(tensor_not_initialized());

    if (row != rhs.row || col != rhs.col || dep != rhs.dep) {
        throw(dimension_mismatch());
    }

    Tensor newTensor{*this};
    for (int i = 0; i < row * col * dep; i++) {
        newTensor.pimpl->data[i] = pimpl->data[i] / rhs.pimpl->data[i];
    }

    return newTensor;
}

/**
 * Operator overloading +
 * 
 * It performs the point-wise sum between a Tensor and a constant
 * 
 * result(i,j,k)=this(i,j,k)+rhs
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator+(const float& rhs) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    Tensor newTensor{*this};
    for (int i = 0; i < row * col * dep; i++) {
        newTensor.pimpl->data[i] += rhs;
    }

    return newTensor;
}

/**
 * Operator overloading - 
 * 
 * It performs the point-wise difference between a Tensor and a constant
 * 
 * result(i,j,k)=this(i,j,k)-rhs
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator-(const float& rhs) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    Tensor newTensor{*this};
    for (int i = 0; i < row * col * dep; i++) {
        newTensor.pimpl->data[i] -= rhs;
    }

    return newTensor;
}

/**
 * Operator overloading *
 * 
 * It performs the point-wise product between a Tensor and a constant
 * 
 * result(i,j,k)=this(i,j,k)*rhs
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator*(const float& rhs) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    Tensor newTensor{*this};
    for (int i = 0; i < row * col * dep; i++) {
        newTensor.pimpl->data[i] *= rhs;
    }

    return newTensor;
}

/**
 * Operator overloading / between a Tensor and a constant
 * 
 * It performs the point-wise division between a Tensor and a constant
 * 
 * result(i,j,k)=this(i,j,k)/rhs
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator/(const float& rhs) const {
    if (!pimpl)
        throw(tensor_not_initialized());

    Tensor newTensor{*this};
    for (int i = 0; i < row * col * dep; i++) {
        newTensor.pimpl->data[i] /= rhs;
    }

    return newTensor;
}

/** 
 * Rows 
 * 
 * @return the number of rows in the tensor
 */
int Tensor::rows() const {
    return row;
}

/** 
 * Cols 
 * 
 * @return the number of columns in the tensor
 */
int Tensor::cols() const {
    return col;
}

/** 
 * Depth 
 * 
 * @return the depth of the tensor
 */
int Tensor::depth() const {
    return dep;
}
=======
>>>>>>> 5994379c2e9a2fef0f3b6192386a473179f024f1

/**
 * Random Initialization
 * 
 * Perform a random initialization of the tensor
 * 
 * @param mean The mean
 * @param std  Standard deviation
 */
<<<<<<< HEAD
void Tensor::init_random(float mean, float std) {
    if (pimpl) {
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean, std);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                for (int k = 0; k < dep; k++) {
                    this->operator()(i, j, k) = distribution(generator);
                }
            }
        }
    } else {
        throw(tensor_not_initialized());
    }
}

/** 
 * showSize
 * 
 * shows the dimensions of the tensor on the standard output.
 * 
 * The format is the following:
 * rows" x "colums" x "depth
 * 
 */
void Tensor::showSize() const {
    cout << this->rows() << " " << this->cols() << " " << this->depth() << endl;
}

/**
 * Reading from file
 * 
 * Load the content of a tensor from a textual file.
 * 
 * The file should have this structure: the first three lines provide the dimensions while 
 * the following lines contains the actual data by channel.
 * 
 * For example, a tensor of size 4x3x2 will have the following structure:
 * 4
 * 3
 * 2
 * data(0,0,0)
 * data(0,1,0)
 * data(0,2,0)
 * data(1,0,0)
 * data(1,1,0)
 * .
 * .
 * .
 * data(3,1,1)
 * data(3,2,1)
 * 
 * if the file is not reachable throw unable_to_read_file()
 * 
 * @param filename the filename where the tensor is stored
 */
void Tensor::read_file(string filename) {
    ifstream ifs{filename};

    if (!ifs) throw(unable_to_read_file());

    ifs >> row;
    ifs >> col;
    ifs >> dep;

    init(row, col, dep);

    for (int k = 0; k < dep; k++) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                ifs >> pimpl->matrix_p[i][j][k];
            }
        }
    }
}

/**
 * Write the tensor to a file
 * 
 * Write the content of a tensor to a textual file.
 * 
 * The file should have this structure: the first three lines provide the dimensions while 
 * the following lines contains the actual data by channel.
 * 
 * For example, a tensor of size 4x3x2 will have the following structure:
 * 4
 * 3
 * 2
 * data(0,0,0)
 * data(0,1,0)
 * data(0,2,0)
 * data(1,0,0)
 * data(1,1,0)
 * .
 * .
 * .
 * data(3,1,1)
 * data(3,2,1)
 * 
 * 
 * @param filename the filename where the tensor is stored
 */
void Tensor::write_file(string filename) const {
    ofstream ost{filename};

    ost << row << "\n"
        << col << "\n"
        << dep << "\n";

    for (int i = 0; i < dep; i++)
        for (int j = 0; j < row; j++)
            for (int k = 0; k < col; k++)
                ost << pimpl->matrix_p[j][k][i] << "\n";
}

/**
 * funzione usata per inizializzare un tensore filtro data una matrice di float
 * usato in daisgram per i filtri per la convolve
 * 
 * @param f 
 * @param row 
 * @param col 
 */
void Tensor::init_filter(float* f, int row, int col) {
    init(row, col, 1);

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            pimpl->matrix_p[i][j][0] = f[i * col + j];
        }
    }
}
=======
void Tensor::init_random(float mean, float std){
    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }    

    }else{
        throw(tensor_not_initialized());
    }
}
>>>>>>> 5994379c2e9a2fef0f3b6192386a473179f024f1
