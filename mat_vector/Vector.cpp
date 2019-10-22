#include <cstddef>
#include "Vector.h"
#include "Matrix.h"
#include <cmath>

using namespace mat_vec;

Vector mat_vec::operator*(double k, const Vector &v) {
    Vector vectorTemplate(v);
    for (size_t i = 0; i < vectorTemplate.size(); i++) {
        vectorTemplate[i] *= k;
    }
    return vectorTemplate;
}

Vector::Vector(size_t size, double value) {
    sizeMathVector = size;
    mathVector = new double[sizeMathVector];

    for (int i = 0; i < sizeMathVector; i++) {
        mathVector[i] = value;
    }
}

Vector::Vector(const Vector &src) {
    sizeMathVector = src.sizeMathVector;
    mathVector = new double[sizeMathVector];

    for (int i = 0; i < sizeMathVector; i++) {
        mathVector[i] = src.mathVector[i];
    }
}

Vector &Vector::operator=(const Vector &rhs) {
    delete[] mathVector;
    sizeMathVector = rhs.sizeMathVector;
    mathVector = new double[sizeMathVector];

    for (int i = 0; i < rhs.sizeMathVector; i++) {
        mathVector[i] = rhs.mathVector[i];
    }
    return *this;
}

Vector::~Vector() {
    delete[] mathVector;
    sizeMathVector = 0;
}

size_t Vector::size() const {
    return sizeMathVector;
}

double Vector::operator[](size_t n) const {
    return mathVector[n];
}

double &Vector::operator[](size_t n) {
    return mathVector[n];
}

double Vector::norm() const {
    double sum = 0;
    for (int i = 0; i < sizeMathVector; i++) {
        sum += mathVector[i] * mathVector[i];
    }
    return sqrt(sum);
}

Vector Vector::normalized() const {
    Vector vectorTemplate(*this);
    double normValue = vectorTemplate.norm();
    for (int i = 0; i < sizeMathVector; i++) {
        vectorTemplate.mathVector[i] = mathVector[i] / normValue;
    }
    return vectorTemplate;
}

void Vector::normalize() {
    double normValue = norm();
    for (int i = 0; i < sizeMathVector; i++) {
        mathVector[i] /= normValue;
    }
}

Vector Vector::operator+(const Vector &rhs) const {
    Vector vectorTemplate(*this);
    for (int i = 0; i < vectorTemplate.sizeMathVector; i++) {
        vectorTemplate.mathVector[i] += rhs.mathVector[i];
    }
    return vectorTemplate;
}

Vector &Vector::operator+=(const Vector &rhs) {
    for (int i = 0; i < sizeMathVector; i++) {
        mathVector[i] += rhs.mathVector[i];
    }
    return *this;
}

Vector Vector::operator-(const Vector &rhs) const {
    Vector vectorTemplate(*this);
    for (int i = 0; i < sizeMathVector; i++) {
        vectorTemplate.mathVector[i] -= rhs.mathVector[i];
    }
    return vectorTemplate;
}

Vector &Vector::operator-=(const Vector &rhs) {
    for (int i = 0; i < sizeMathVector; i++) {
        mathVector[i] -= rhs.mathVector[i];
    }
    return *this;
}

Vector Vector::operator^(const Vector &rhs) const {
    Vector vectorTemplate(*this);
    for (int i = 0; i < sizeMathVector; i++) {
        vectorTemplate.mathVector[i] *= rhs.mathVector[i];
    }
    return vectorTemplate;
}

Vector &Vector::operator^=(const Vector &rhs) {
    for (int i = 0; i < sizeMathVector; i++) {
        mathVector[i] *= rhs.mathVector[i];
    }
    return *this;
}

double Vector::operator*(const Vector &rhs) const {
    double product = 0;
    for (int i = 0; i < sizeMathVector; i++)
        product += mathVector[i] * rhs.mathVector[i];
    return product;
}

Vector Vector::operator*(double k) const {
    Vector vectorTemplate(*this);
    for (int i = 0; i < sizeMathVector; i++) {
        vectorTemplate.mathVector[i] *= k;
    }
    return vectorTemplate;
}

Vector &Vector::operator*=(double k) {
    for (int i = 0; i < sizeMathVector; i++) {
        mathVector[i] *= k;
    }
    return *this;
}

Vector Vector::operator/(double k) const {
    Vector vectorTemplate(*this);
    for (int i = 0; i < sizeMathVector; i++) {
        vectorTemplate.mathVector[i] /= k;
    }
    return vectorTemplate;
}

Vector &Vector::operator/=(double k) {
    for (int i = 0; i < sizeMathVector; i++) {
        mathVector[i] /= k;
    }
    return *this;
}

Vector Vector::operator*(const Matrix &mat) const {
    Vector vector(mat.shape().second);
    for (int i = 0; i < vector.sizeMathVector; i++) {
        double result = 0.0;
        for (int j = 0; j < this->sizeMathVector; j++) {
            double value = mat.get(j, i);
            result += this->mathVector[j] * value;
        }
        vector[i] = result;
    }
    return vector;
}

Vector &Vector::operator*=(const Matrix &mat) { return *this = *this * mat; }

bool Vector::operator==(const Vector &rhs) const {
    double e = 0.000000000001;
    for (int i = 0; i < sizeMathVector; i++) {
        if (std::abs(mathVector[i] - rhs.mathVector[i]) >= e)
            return false;
    }
    return true;
}

bool Vector::operator!=(const Vector &rhs) const {
    double e = 0.000000000001;
    for (int i = 0; i < sizeMathVector; i++) {
        if (std::abs(mathVector[i] - rhs.mathVector[i]) < e)
            return false;
    }
    return true;
}