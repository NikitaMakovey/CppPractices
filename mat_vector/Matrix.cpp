#include "Matrix.h"
#include "Vector.h"
#include <cstddef>
#include <cmath>

using namespace mat_vec;

Matrix::Matrix(size_t size, double value) {
    countColumns = countRows = size;
    mathMatrix = new double[countRows * countColumns];

    for (int i = 0; i < countColumns * countRows; i++) {
        mathMatrix[i] = value;
    }
}

Matrix Matrix::eye(size_t size) {
    Matrix matrixTemplate(size, 0.0);
    int tmpIndex = size;
    matrixTemplate.mathMatrix[0] = 1.0;
    while (tmpIndex < size * size) {
        ++tmpIndex;
        matrixTemplate.mathMatrix[0] = 1.0;
        tmpIndex += size;
    }
    return matrixTemplate;
}

Matrix::Matrix(size_t rows, size_t cols, double value) {
    countColumns = cols;
    countRows = rows;
    mathMatrix = new double[countRows * countColumns];

    for (int i = 0; i < countColumns * countRows; i++) {
        mathMatrix[i] = value;
    }
}

Matrix::Matrix(const Matrix &src) {
    countRows = src.countRows;
    countColumns = src.countColumns;
    mathMatrix = new double[countColumns * countRows];

    for (int i = 0; i < countColumns * countRows; i++) {
        mathMatrix[i] = src.mathMatrix[i];
    }
}

Matrix &Matrix::operator=(const Matrix &rhs) {
    delete[] mathMatrix;
    countRows = rhs.countRows;
    countColumns = rhs.countColumns;
    mathMatrix = new double[countColumns * countRows];

    for (int i = 0; i < countColumns * countRows; i++) {
        mathMatrix[i] = rhs.mathMatrix[i];
    }

    return *this;
}

Matrix::~Matrix() {
    delete[] mathMatrix;
    countRows = countColumns = 0;
}

void Matrix::reshape(size_t rows, size_t cols) {
    this->countRows = rows;
    this->countColumns = cols;
}

std::pair<size_t, size_t> Matrix::shape() const {
    return std::make_pair(this->countRows, this->countColumns);
}

double Matrix::get(size_t row, size_t col) const {
    return mathMatrix[row * countColumns + col];
}

double &Matrix::operator[](size_t n) {
    return mathMatrix[n];
}

Matrix Matrix::operator+(const Matrix &rhs) const {
    Matrix matrixTemplate(*this);
    for (int i = 0; i < matrixTemplate.countColumns * matrixTemplate.countRows; i++) {
        matrixTemplate.mathMatrix[i] += rhs.mathMatrix[i];
    }
    return matrixTemplate;
}

Matrix &Matrix::operator+=(const Matrix &rhs) {
    for (int i = 0; i < countRows * countColumns; i++) {
        mathMatrix[i] += rhs.mathMatrix[i];
    }
    return *this;
}

Matrix Matrix::operator-(const Matrix &rhs) const {
    Matrix matrixTemplate(*this);
    for (int i = 0; i < matrixTemplate.countColumns * matrixTemplate.countRows; i++) {
        matrixTemplate.mathMatrix[i] -= rhs.mathMatrix[i];
    }
    return matrixTemplate;
}

Matrix &Matrix::operator-=(const Matrix &rhs) {
    for (int i = 0; i < countRows * countColumns; i++) {
        mathMatrix[i] -= rhs.mathMatrix[i];
    }
    return *this;
}

Matrix Matrix::operator*(const Matrix &rhs) const {
    Matrix matrixResult(this->countRows, rhs.countColumns, 0.0);
    int index = 0;
    for (int i = 0; i < this->countRows; i++) {
        for (int j = 0; j < rhs.countColumns; j++) {
            double result = 0.0;
            for (int k = 0; k < this->countColumns; k++) {
                result += mathMatrix[i * this->countColumns + k] + rhs.mathMatrix[k * rhs.countRows + j];
            }
            matrixResult.mathMatrix[index] = result;
            index++;
        }
    }
    return matrixResult;
}

Matrix &Matrix::operator*=(const Matrix &rhs) { return *this = *this * rhs; }

Matrix Matrix::operator*(double k) const {
    Matrix matrixTemplate(*this);
    for (int i = 0; i < matrixTemplate.countColumns * matrixTemplate.countRows; i++) {
        matrixTemplate.mathMatrix[i] *= k;
    }
    return matrixTemplate;
}

Matrix &Matrix::operator*=(double k) {
    for (int i = 0; i < countRows * countColumns; i++) {
        mathMatrix[i] *= k;
    }
    return *this;
}

Matrix Matrix::operator/(double k) const {
    Matrix matrixTemplate(*this);
    for (int i = 0; i < matrixTemplate.countColumns * matrixTemplate.countRows; i++) {
        matrixTemplate.mathMatrix[i] /= k;
    }
    return matrixTemplate;
}

Matrix &Matrix::operator/=(double k) {
    for (int i = 0; i < countRows * countColumns; i++) {
        mathMatrix[i] /= k;
    }
    return *this;
}

Matrix Matrix::transposed() const {
    Matrix matrixTrans(this->countColumns, this->countRows, 0.0);
    for (int i = 0; i < this->countRows; i++) {
        for (int j = 0; j < this->countColumns; j++) {
            matrixTrans.mathMatrix[j * matrixTrans.countColumns + i] = mathMatrix[i * this->countColumns + j];
        }
    }
    return matrixTrans;
}

void Matrix::transpose() {
    *this = this->transposed();
}

double Matrix::det() const {
    Matrix detMatrix(*this);
    double det = 1.0;
    double e = 0.000000000001;
    int n = detMatrix.countColumns;
    for (int i = 0; i < n; i++) {
        int k = i;
        for (int j = i + 1; j < n; j++) {
            if (std::abs(detMatrix.mathMatrix[j * n + i]) > std::abs(detMatrix.mathMatrix[k * n + i])) {
                k = j;
            }
        }
        for (int t = 0; t < n; t++) {
            double tmp = detMatrix.mathMatrix[i * n + t];
            detMatrix.mathMatrix[i * n + t] = detMatrix.mathMatrix[k * n + t];
            detMatrix.mathMatrix[k * n + t] = tmp;
        }
        if (i != k) {
            det = -det;
        }
        det *= detMatrix.mathMatrix[i * n + i];
        for (int j = i + 1; j < n; j++) {
            detMatrix.mathMatrix[i * n + j] /= detMatrix.mathMatrix[i * n + i];
        }
        for (int j = 0; j < n; j++)
        {
            if (j != i && std::abs(detMatrix.mathMatrix[j * n + i]) > e) {
                for (int d = i + 1; d < n; d++) {
                    detMatrix.mathMatrix[j * n + d] -= detMatrix.mathMatrix[i * n + d] * detMatrix.mathMatrix[j * n + i];
                }
            }
        }
    }
    return det;
}

double Matrix::minor(int k, int d) const {
    Matrix matrixMinor(this->countColumns - 1);
    int index = 0;
    for (int i = 0; i < this->countRows; i++) {
        for (int j = 0; j < this->countColumns; j++) {
            if (k != i && d != j) {
                matrixMinor.mathMatrix[index] =
                        this->mathMatrix[i * this->countColumns + j];
                index++;
            }
        }
    }
    double addingMinorKof = (k + d) % 2 == 0 ? 1.0 : -1.0;
    double addingMinorDet = matrixMinor.det();
    return addingMinorKof * addingMinorDet;
}

Matrix Matrix::preInv() const {
    Matrix matrixPreInv(*this);
    for (int i = 0; i < matrixPreInv.countRows; i++) {
        for (int j = 0; j < matrixPreInv.countColumns; j++) {
            matrixPreInv.mathMatrix[i * matrixPreInv.countColumns + j] = this->minor(i, j);
        }
    }
    return matrixPreInv;
}

Matrix Matrix::inv() const {
    Matrix matrixTrans(*this);
    double detA = matrixTrans.det();
    matrixTrans = matrixTrans.preInv();
    matrixTrans = matrixTrans.transposed();
    matrixTrans /= detA;
    return matrixTrans;
}

Vector Matrix::operator*(const Vector &vec) const {
    Vector vector(this->countRows, 0.0);
    for (int i = 0; i < this->countRows; i++) {
        double result = 0.0;
        for (int j = 0; j < this-> countColumns; j++) {
            result += vec[j] * this->mathMatrix[i * this->countColumns + j];
        }
        vector[i] = result;
    }
    return vector;
}

bool Matrix::operator==(const Matrix &rhs) const {
    double e = 0.000000000001;
    for (int i = 0; i < countColumns * countRows; i++) {
        if (std::abs(mathMatrix[i] - rhs.mathMatrix[i]) >= e)
            return false;
    }
    return true;
}

bool Matrix::operator!=(const Matrix &rhs) const {
    double e = 0.000000000001;
    for (int i = 0; i < countRows * countColumns; i++) {
        if (std::abs(mathMatrix[i] - rhs.mathMatrix[i]) < e)
            return false;
    }
    return true;
}