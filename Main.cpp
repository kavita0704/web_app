#include <bits/stdc++.h>
using namespace std;

class BigInteger {
private:
    vector<int> digits; // store digits in reverse order

public:
    // Constructor from string
    BigInteger(string num) {
        for (int i = num.size() - 1; i >= 0; i--) {
            if (!isdigit(num[i])) {
                throw invalid_argument("Invalid number");
            }
            digits.push_back(num[i] - '0');
        }
    }

    // Default constructor
    BigInteger() {
        digits.push_back(0);
    }

    // Remove leading zeros
    void normalize() {
        while (digits.size() > 1 && digits.back() == 0) {
            digits.pop_back();
        }
    }

    // Print number
    void print() const {
        for (int i = digits.size() - 1; i >= 0; i--) {
            cout << digits[i];
        }
        cout << endl;
    }

    // Addition
    BigInteger add(const BigInteger &other) {
        BigInteger result;
        result.digits.clear();

        int carry = 0;
        int n = digits.size();
        int m = other.digits.size();

        for (int i = 0; i < max(n, m); i++) {
            int sum = carry;

            if (i < n) sum += digits[i];
            if (i < m) sum += other.digits[i];

            result.digits.push_back(sum % 10);
            carry = sum / 10;
        }

        if (carry) result.digits.push_back(carry);

        return result;
    }

    // Multiplication
    BigInteger multiply(const BigInteger &other) {
        BigInteger result;
        int n = digits.size();
        int m = other.digits.size();

        result.digits.assign(n + m, 0);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result.digits[i + j] += digits[i] * other.digits[j];
            }
        }

        // Handle carry
        int carry = 0;
        for (int i = 0; i < result.digits.size(); i++) {
            int temp = result.digits[i] + carry;
            result.digits[i] = temp % 10;
            carry = temp / 10;
        }

        result.normalize();
        return result;
    }
};

int main() {
    BigInteger num1("12345");
    BigInteger num2("6789");

    cout << "First number: ";
    num1.print();

    cout << "Second number: ";
    num2.print();

    BigInteger sum = num1.add(num2);
    cout << "Sum: ";
    sum.print();

    BigInteger product = num1.multiply(num2);
    cout << "Product: ";
    product.print();

    return 0;
}