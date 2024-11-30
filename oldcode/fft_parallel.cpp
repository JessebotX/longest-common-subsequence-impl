#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <thread>
#include <string>
#include <exception>

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

// Find next power of 2
size_t nextPowerOf2(size_t n) {
    size_t power = 1;
    while (power < n) {
        power *= 2;
    }
    return power;
}

// Worker thread task
void worker() {
    while (true) {

    }
}

// Iterative FFT implementation
void fft_parallel(uint nThreads, std::vector<std::complex<float>>& x) {
    

    const size_t N = x.size();
    const int n = log2(N);

    // Bit-reverse reordering
    for (size_t i = 0; i < N; i++) {
        size_t j = bitReverse(i, n);
        if (i < j) {
            std::swap(x[i], x[j]);
        }
    }

    // Create worker threads
    std::vector<std::thread> threads;
    for (int i = 0; i < nThreads; i++) {
        // each worker thread calls worker() and waits for work
    }

    // Start timer

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

// Calculate magnitude of a complex number
float magnitude(const std::complex<float>& c) {
    return std::sqrt(c.real() * c.real() + c.imag() * c.imag());
}

// Calculate phase in degrees of a complex number
float phase_degrees(const std::complex<float>& c) {
    return std::atan2(c.imag(), c.real()) * 180.0f / static_cast<float>(M_PI);
}



// Main FFT function that handles the input array with configurable FFT size
std::vector<std::complex<float>> fft(uint nThreads, const std::vector<float>& input, size_t fft_size = 0) {
    size_t n = input.size();
    
    // If fft_size is 0, use next power of 2 of input size
    // If fft_size is specified, ensure it's a power of 2 and >= input size
    if (fft_size == 0) {
        fft_size = nextPowerOf2(n);
    } else {
        if (!isPowerOfTwo(fft_size)) {
            throw std::runtime_error("FFT size must be a power of 2");
        }
        if (fft_size < n) {
            throw std::runtime_error("FFT size must be >= input size");
        }
    }

    // Convert input to complex numbers and pad with zeros if necessary
    std::vector<std::complex<float>> data(fft_size);
    for (size_t i = 0; i < fft_size; i++) {
        data[i] = std::complex<float>(i < n ? input[i] : 0.0f, 0.0f);
    }

    // Perform FFT
    fft_parallel(nThreads, data);
    return data;
}

// Helper function to print complex numbers with frequency information
void printFFTResults(const std::vector<std::complex<float>>& data, float sampling_freq) {
    const size_t N = data.size();
    const float freq_resolution = sampling_freq / N;
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nFFT Results (N=" << N << ", Resolution=" << freq_resolution << " Hz):\n";
    std::cout << std::setw(15) << "Frequency (Hz)" 
              << std::setw(15) << "Magnitude" 
              << std::setw(15) << "Phase (deg)"
              << std::setw(25) << "Complex Value" << "\n";
    std::cout << std::string(70, '-') << "\n";

    // Only print up to Nyquist frequency (N/2 points)
    for (size_t i = 0; i <= N/2; i++) {
        float frequency = i * freq_resolution;
        float mag = magnitude(data[i]);
        float phase = phase_degrees(data[i]);
        
        std::cout << std::setw(15) << frequency 
                  << std::setw(15) << mag
                  << std::setw(15) << phase
                  << std::setw(25) << data[i] << "\n";
    }
}

int main(int argc, char* argv[]) {

    // The FFT converts a signal into frequencies contained in the signal
    // Imagine an audio file, which is just a bunch of amplitude values over time (a signal)
    // The FFT takes that as an input, and spits out what frequencies of sine waves are present in that window of time

    // Test with a simple signal
    std::vector<float> input = {1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f}; // imagine this is a slice of an audio file
    float sampling_freq = 1000.0f;  // 1 kHz sampling rate (if this was an audio recording, then 1k audio samples recorded per second)
                                    // in real audio this would be 44100 Hz normally
    size_t fft_size = 16;           // Use 16-point FFT (increased frequency resolution) (also known as 'bin size')
    
    // Nyquist frequency = Sampling Frequency / 2
    // the output of the FFT is all the amplitudes and phases of all the frequencies in the signal from 0hz to the Nyquist Frequency in steps of (Sampling Freq / Bin Size)

    /*
    Online FFT calculators may show outputs of FFT above the nyquist frequency, but this code only shows up to the nyquist frequency

    "While my code was only showing up to the Nyquist frequency. The full FFT output is actually symmetric (conjugate symmetric for real inputs).
    For a real input signal, the FFT output has this symmetry property:

    First half (0 to 500 Hz): The frequencies we normally care about
    Second half (500 to 1000 Hz): Mirror of the first half (complex conjugates)"
    - Claude's explanation for why that is */ 

    int nThreads = 1;  // Default number of threads is 1
    if (argc > 1) {
        try {
            nThreads = std::stoi(argv[1]);
            if (nThreads <= 0) {
                throw std::invalid_argument("Number of threads must be positive.");
            }
        } catch (const std::exception& e) {
            std::cout << "Error: " << e.what() << std::endl;
            return 1;
        }
    }

    try {
        std::cout << "Original signal:\n";
        for (float f : input) {
            std::cout << f << " ";
        }
        std::cout << "\n";
        
        auto result = fft(nThreads, input, fft_size);
        printFFTResults(result, sampling_freq);
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    return 0;
}