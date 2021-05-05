#include "tensor.h"

struct Tensor::Impl {
    float* data;
    float** cols_p;
    float*** matrix_p;
};

void Tensor::allocate_matrix(size_t row, size_t col, size_t dep) {
    pimpl->data = new float[row * col * dep];
    pimpl->cols_p = new float*[row * col];
    pimpl->matrix_p = new float**[row];

    for (size_t i = 0; i < row; i++) {
        pimpl->matrix_p[i] = &(pimpl->cols_p[i * col]);

        for (size_t j = 0; j < col; j++)
            pimpl->matrix_p[i][j] = pimpl->data + i * col * dep + j * dep;
    }
}

Tensor::Tensor(const Tensor& that) {
    pimpl = new Impl;

    row = that.row;
    col = that.col;
    dep = that.dep;

    allocate_matrix(row, col, dep);

    for (size_t i = 0; i < row * col * dep; i++) {
        pimpl->data[i] = that.pimpl->data[i];
    }
}

Tensor::Tensor(int r, int c, int d, float v) {
    init(r, c, d, v);
}

Tensor::~Tensor() {
    delete[] pimpl->data;
    delete[] pimpl->cols_p;
    delete[] pimpl->matrix_p;
    delete pimpl;
}

void Tensor::init(int r, int c, int d, float v) {
    if (r < 0 || c < 0 || d < 0)
        throw(invalid_parameter());

    pimpl = new Impl;
    allocate_matrix(r, c, d);

    row = r;
    col = c;
    dep = d;

    for (size_t i = 0; i < row * col * dep; i++) {
        pimpl->data[i] = v;
    }
}

float Tensor::getMin(int k) {
    float min = pimpl->data[0];

    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++) {
            if (pimpl->matrix_p[i][j][k] < min) min = pimpl->matrix_p[i][j][k];
        }
    }

    return min;
}

float Tensor::getMax(int k) {
    float max = pimpl->data[0];

    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++) {
            if (pimpl->matrix_p[i][j][k] > max) max = pimpl->matrix_p[i][j][k];
        }
    }

    return max;
}

float Tensor::operator()(int i, int j, int k) const {
    if (i < 0 || i >= row || j < 0 || j >= col || k < 0 || k >= dep)
        throw(index_out_of_bound());

    return pimpl->matrix_p[i][j][k];
}

float& Tensor::operator()(int i, int j, int k) {
    if (i < 0 || i >= row || j < 0 || j >= col || k < 0 || k >= dep)
        throw(index_out_of_bound());

    return pimpl->matrix_p[i][j][k];
}

ostream& operator<<(ostream& stream, const Tensor& obj) {
    for (size_t k = 0; k < obj.dep; k++) {
        for (size_t i = 0; i < obj.row; i++) {
            for (size_t j = 0; j < obj.col; j++) {
                stream << obj.pimpl->matrix_p[i][j][k] << " ";
            }

            stream << "\n";
        }

        stream << "\n";
    }

    return stream;
}

void Tensor::clamp(float low, float high) {
    for (size_t i = 0; i < row * col * dep; i++) {
        if (pimpl->data[i] < low)
            pimpl->data[i] = low;
        else if (pimpl->data[i] > high)
            pimpl->data[i] = high;
    }
}

void Tensor::rescale(float new_max) {
    for (size_t k = 0; k < dep; k++) {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < col; j++) {
                pimpl->matrix_p[i][j][k] = ((pimpl->matrix_p[i][j][k] - getMin(k)) / (getMax(k) - getMin(k))) * new_max;
            }
        }
    }
}

Tensor Tensor::padding(int pad_h, int pad_w) {
    Tensor new_t{row + pad_h * 2, col + pad_w * 2, dep};

    for (size_t k = 0; k < new_t.dep; k++) {
        for (size_t i = 0; i < new_t.row; i++) {
            for (size_t j = 0; j < new_t.col; j++) {
                if (i >= pad_h && i < new_t.row - pad_h && j >= pad_w && j < new_t.col - pad_w)
                    new_t(i, j, k) = pimpl->matrix_p[i - pad_h][j - pad_w][k];
                else
                    new_t(i, j, k) = 0;
            }
        }
    }

    return new_t;
}

Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end) {
    if (row_start < 0 || row_start > row + 1 || row_end < 0 || row_end > row + 1 ||
        col_start < 0 || col_start > col + 1 || col_end < 0 || col_end > col + 1 ||
        depth_start < 0 || depth_start > dep + 1 || depth_end < 0 || depth_end > dep + 1)
        throw(index_out_of_bound());

    if (row_end >= row_start || col_end >= col_start || depth_end >= depth_start)
        throw(invalid_parameter());

    Tensor new_t{row_end - row_start, col_end - col_start, depth_end - depth_start};

    for (size_t k = depth_start; k < depth_end; k++) {
        for (size_t i = row_start; i < row_end; i++) {
            for (size_t j = col_start; j < col_end; j++) {
                new_t(i - row_start, j - col_start, k - depth_start) = pimpl->matrix_p[i][j][k];
            }
        }
    }

    return new_t;
}

Tensor Tensor::concat(const Tensor& rhs, int axis) {
    Tensor new_t;

    switch (axis) {
        case 0:
            if (col != rhs.col || dep != rhs.dep)
                throw(concat_wrong_dimension());

            new_t.init(row + rhs.row, col, dep);
            for (size_t k = 0; k < dep; k++) {
                for (size_t i = 0; i < row + rhs.row; i++) {
                    for (size_t j = 0; j < col; j++) {
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
            for (size_t k = 0; k < dep; k++) {
                for (size_t i = 0; i < row; i++) {
                    for (size_t j = 0; j < col + rhs.col; j++) {
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
            for (size_t k = 0; k < dep + rhs.dep; k++) {
                for (size_t i = 0; i < row; i++) {
                    for (size_t j = 0; j < col; j++) {
                        if (k < dep)
                            new_t(i, j, k) = pimpl->matrix_p[i][j][k];
                        else
                            new_t(i, j, k) = rhs.pimpl->matrix_p[i][j][k - dep];
                    }
                }
            }

            break;
        default:
            throw(invalid_parameter());
    }

    return new_t;
}

Tensor operator-(Tensor lhs, const Tensor& rhs) {
    if (lhs.row != rhs.row || lhs.col != rhs.col || lhs.dep != rhs.dep)
        throw(dimension_mismatch());

    for (size_t i = 0; i < lhs.row * lhs.col * lhs.dep; i++) {
        lhs.pimpl->data[i] -= rhs.pimpl->data[i];
    }

    return lhs;
}

Tensor operator+(Tensor lhs, const Tensor& rhs) {
    if (lhs.row != rhs.row || lhs.col != rhs.col || lhs.dep != rhs.dep)
        throw(dimension_mismatch());

    for (size_t i = 0; i < lhs.row * lhs.col * lhs.dep; i++) {
        lhs.pimpl->data[i] += rhs.pimpl->data[i];
    }

    return lhs;
}

Tensor operator*(Tensor lhs, const Tensor& rhs) {
    if (lhs.row != rhs.row || lhs.col != rhs.col || lhs.dep != rhs.dep)
        throw(dimension_mismatch());

    for (size_t i = 0; i < lhs.row * lhs.col * lhs.dep; i++) {
        lhs.pimpl->data[i] *= rhs.pimpl->data[i];
    }

    return lhs;
}

Tensor operator/(Tensor lhs, const Tensor& rhs) {
    if (lhs.row != rhs.row || lhs.col != rhs.col || lhs.dep != rhs.dep)
        throw(dimension_mismatch());

    for (size_t i = 0; i < lhs.row * lhs.col * lhs.dep; i++) {
        lhs.pimpl->data[i] /= rhs.pimpl->data[i];
    }

    return lhs;
}

Tensor operator-(Tensor lhs, const float& rhs) {
    for (size_t i = 0; i < lhs.row * lhs.col * lhs.dep; i++) {
        lhs.pimpl->data[i] -= rhs;
    }

    return lhs;
}

Tensor operator+(Tensor lhs, const float& rhs) {
    for (size_t i = 0; i < lhs.row * lhs.col * lhs.dep; i++) {
        lhs.pimpl->data[i] += rhs;
    }

    return lhs;
}

Tensor operator*(Tensor lhs, const float& rhs) {
    for (size_t i = 0; i < lhs.row * lhs.col * lhs.dep; i++) {
        lhs.pimpl->data[i] *= rhs;
    }

    return lhs;
}

Tensor operator/(Tensor lhs, const float& rhs) {
    for (size_t i = 0; i < lhs.row * lhs.col * lhs.dep; i++) {
        lhs.pimpl->data[i] /= rhs;
    }

    return lhs;
}

Tensor& Tensor::operator=(const Tensor& other) {
    if (this != &other) {
        if (pimpl) {
            delete[] pimpl->data;
            delete[] pimpl->cols_p;
            delete[] pimpl->matrix_p;
            delete pimpl;
        }

        init(other.row, other.col, other.dep);
        for (size_t i = 0; i < row * col * dep; i++) {
            pimpl->data[i] = other.pimpl->data[i];
        }
    }

    return (*this);
}

Tensor Tensor::convolve(const Tensor& f) {
    Tensor new_t{row, col, dep};

    int pad_h = (f.row - 1) / 2;
    int pad_w = (f.col - 1) / 2;

    Tensor this_padded = padding(pad_h, pad_w);

    for (size_t k = 0; k < dep; k++)
    {
           
    }

    return new_t;
}

int Tensor::rows() {
    return row;
}

int Tensor::cols() {
    return col;
}

int Tensor::depth() {
    return dep;
}