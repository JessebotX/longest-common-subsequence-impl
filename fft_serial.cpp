#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <iostream>

// Bit-reverse an integer's lowest n bits
unsigned int bitReverse(unsigned int x, int n) {
    unsigned int result = 0;
    for (int i = 0; i < n; i++) {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

// Helper function to check if number is power of 2
bool isPowerOfTwo(size_t n) {
    return n && !(n & (n - 1));
}

// Calculate log2 of a number that is power of 2
int log2(size_t n) {
    int result = 0;
    while (n >>= 1) ++result;
    return result;
}

// Iterative FFT implementation
void fft_iterative(std::vector<std::complex<float>>& x) {
    const size_t N = x.size();
    const int n = log2(N);

    // Bit-reverse reordering
    for (size_t i = 0; i < N; i++) {
        size_t j = bitReverse(i, n);
        if (i < j) {
            std::swap(x[i], x[j]);
        }
    }

    // Cooley-Tukey FFT
    for (int s = 1; s <= n; s++) {
        size_t m = 1 << s;  // 2^s
        float angle = -2.0f * static_cast<float>(M_PI) / static_cast<float>(m);
        std::complex<float> wm = std::polar(1.0f, angle);
        
        for (size_t k = 0; k < N; k += m) {
            std::complex<float> w = 1.0f;
            for (size_t j = 0; j < m/2; j++) {
                std::complex<float> t = w * x[k + j + m/2];
                std::complex<float> u = x[k + j];
                x[k + j] = u + t;
                x[k + j + m/2] = u - t;
                w *= wm;
            }
        }
    }
}

// Main FFT function that handles the input array
std::vector<std::complex<float>> fft(const std::vector<float>& input) {
    size_t n = input.size();
    
    // Check if input size is power of 2
    if (!isPowerOfTwo(n)) {
        throw std::runtime_error("Input size must be a power of 2");
    }

    // Convert input to complex numbers
    std::vector<std::complex<float>> data(n);
    for (size_t i = 0; i < n; i++) {
        data[i] = std::complex<float>(input[i], 0.0f);
    }

    // Perform FFT
    fft_iterative(data);
    return data;
}

// Helper function to print complex numbers
void printComplex(const std::vector<std::complex<float>>& data) {
    for (const auto& c : data) {
        std::cout << c.real() << " + " << c.imag() << "i\n";
    }
}

int main() {
    // Test with a simple signal
    std::vector<float> input = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
    
    try {
        std::cout << "Original signal:\n";
        for (float f : input) {
            std::cout << f << " ";
        }
        std::cout << "\n\nFFT result:\n";
        
        auto result = fft(input);
        printComplex(result);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    return 0;
}