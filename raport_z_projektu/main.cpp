#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <matplot/matplot.h>
#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <tuple>

constexpr double PI = 3.1415926535897932384626433;

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace matplot;


const int RAND = 10;
matplot::figure_handle fig;
// === Proste operacje matematyczne ===
int add(int i, int j) { return i + j; }
int subtract(int i, int j) { return i - j; }

std::vector<double> prep_x(double start, double end, size_t n_samples)
{
    if (n_samples < 2) return {};
    double x_width = (end - start) / (n_samples - 1);
    std::vector<double> x(n_samples);
    for(size_t i = 0 ; i < n_samples; i++)
        x[i] = start + i * x_width;
    return x;
}

// === Rysowanie: dowolny sygnał ===
void plot(const std::vector<double>& x, const std::vector<double>& y) {
    fig = matplot::figure();
    matplot::plot(x, y);
    matplot::title("Signal");
    matplot::xlabel("X");
    matplot::ylabel("Y");
    matplot::show();
}

void download(const std::string name)
{
    if (name.empty()) return;
    std::string file_name = name + ".png";
    fig->save(file_name);
    std::cout << "Wykres zapisany jako: " << file_name << std::endl;
}

double Amplitude(std::vector<double>& signal)
{
    if (signal.empty()) return 0.0;
    double high = *std::max_element(signal.begin(), signal.end());
    double low = *std::min_element(signal.begin(), signal.end());
    return (high - low) / 2.0;
}

std::vector<double> noise(std::vector<double>& signal, double noise_level)
{
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    std::vector<double> y;
    y.reserve(signal.size());
    double amplitude = Amplitude(signal);
    for(size_t i = 0; i < signal.size(); i++) {
        double level = (5 - std::rand()) % RAND;
        while(level == 0.0)
        {
             level = (5 - std::rand()) % RAND;
        }
        double noise_val = (static_cast<double>(amplitude) * noise_level / level);
        y.push_back(signal[i] + noise_val);
    }
    return y;
}


std::vector<double> gen_sin(const std::vector<double>& x, double freq, double phase, double shift_y, double amplitude)
{
    
    size_t n_samples = x.size();
    std::vector<double> y;
    y.reserve(n_samples);
    for(size_t i = 0; i < n_samples; i++) {
        double t = x[i];
        y.push_back(amplitude * std::sin(2 * PI * freq * t - phase * PI) + shift_y);
    }
    return y;
}


std::vector<double> gen_cos(const std::vector<double>& x, double freq, double phase, double shift_y, double amplitude)
{
    size_t n_samples = x.size();
    std::vector<double> y;
    y.reserve(n_samples);
    for(size_t i = 0; i < n_samples; i++) {
        double t = x[i];
        y.push_back(amplitude * std::cos(2 * PI * freq * t - phase * PI) + shift_y);
    }
    return y;
}   




std::vector<double> gen_parabola(const std::vector<double>& x, double a, double p, double q)
{
    size_t n_samples = x.size();
    std::vector<double> y;
    y.reserve(n_samples);
    for(size_t i = 0; i < n_samples; i++)
        y.push_back(a * (x[i] - p) * (x[i] - p) + q);
    return y;
}

std::vector<double> gen_sqr(const std::vector<double>& x, double fill, double freq, double low, double pick)
{
    if(freq == 0.0) return {}; // Prevent division by zero
    if (x.empty()) return {}; // zabezpieczenie przed pustym wektorem
    size_t n_samples = x.size();
    std::vector<double> y;
    y.reserve(n_samples);
    fill = std::clamp(fill, 0.0, 1.0); // zapewnienie, że fill jest w zakresie [0, 1]
    for(size_t i = 0; i < n_samples; i++)
    {
        double t = x[i];
        double period = 1.0 / freq; // okres sygnału
        // obliczenie fazy w zakresie [0, 1)
        double phase = std::fmod(t, period) / period;
        if(phase < 0) phase += 1.0; // zapewnienie, że phase jest w zakresie [0, 1)
        if (phase >= 0 && phase < fill) {
            y.push_back(pick);
        } else {
            y.push_back(low);
        }
    }
    if (y.empty()) return {}; // zabezpieczenie przed pustym wektorem
    return y;
}



std::vector<double> gen_saw(const std::vector<double>& x, double freq, double grade, double low, double pick) {
    if (freq == 0.0) return {}; // Prevent division by zero
    size_t n_samples = x.size();
    std::vector<double> y;
    y.reserve(n_samples);
    grade = std::clamp(grade, 0.0, 1.0); // zapewnienie, że grade jest w zakresie [0, 1]
    for (size_t i = 0; i < n_samples; i++) {
        double t = x[i];
        double period = 1.0 / freq; // okres sygnału
        double up_period = period * grade; // okres narastania
        double down_period = period * (1.0 - grade); // okres opadania
        double phase = std::fmod(t, period) / period; // faza w zakresie [0, 1)
        if (phase < 0) phase += 1.0; // zapewnienie, że phase jest w zakresie [0, 1)
        if (phase < grade && up_period > 0) {
            // Narastanie
            double up_grade = (pick - low) / up_period;
            y.push_back(low + up_grade * (phase * period));
        } else if (down_period > 0) {
            // Opadanie
            double down_grade = (pick - low) / down_period;
            y.push_back(pick - down_grade * ((phase - grade) * period));
        } else {
            y.push_back(low);
        }
    }
    return y;
}




std::vector<std::complex<double>> DFT(const std::vector<double>& signal) {
    int N = signal.size();
    if (N == 0) return {};
    std::vector<std::complex<double>> result(N);

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0, 0);
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * PI * k * n / N;
            sum += signal[n] * std::polar(1.0, angle);
        }
        result[k] = sum;
    }
    return result;
}

void plot_dft(const std::vector<double>& signal, double sampling_rate, const std::string name) {
   std::string full_name = "Widmo amplitudowe (DFT)" + name;
    auto X = DFT(signal);
    size_t N = X.size();
    if (N == 0) return;
    std::vector<double> freqs(N), mag(N);
    for (size_t k = 0; k < N; ++k) {
        freqs[k] = k * sampling_rate / static_cast<double>(N);
        mag[k] = std::abs(X[k]);
    }
    plot(freqs, mag);
    title(full_name);
    xlabel("Częstotliwość [Hz]");
    ylabel("|X(f)|");
    show();
}


std::vector<double> idft(const std::vector<std::complex<double>>& spectrum) {
    int N = spectrum.size();
    std::vector<double> result(N);

    for (int n = 0; n < N; ++n) {
        std::complex<double> sum(0, 0);
        for (int k = 0; k < N; ++k) {
            double angle = 2.0 * PI * k * n / N;
            sum += spectrum[k] * std::polar(1.0, angle);  // e^(j*angle)
        }
        result[n] = sum.real() / N; // Normalizacja wyniku
    }

    return result;
}

// Gaussowskie ważenie w filtrze B-spline
std::vector<double> gaussianWeights(int degree) {
    if (degree <= 0) return {1.0};
    std::vector<double> weights(2 * degree + 1);
    double sigma = std::max(degree / 2.0, 1e-6); // avoid division by zero
    double sum = 0.0;

    for (int i = 0; i < static_cast<int>(weights.size()); ++i) {
        double x = i - degree;
        weights[i] = std::exp(-x * x / (2 * sigma * sigma));
        sum += weights[i];
    }

    if (sum == 0.0) sum = 1.0; // avoid division by zero
    for (double& w : weights) {
        w /= sum; // Normalize weights
    }

    return weights;
}


std::vector<double> filter_signal(const std::vector<double>& signal, int degree) {
    int n = signal.size();
    if (degree < 0 || 2 * degree + 1 > n) {
        // Invalid degree, return original signal
        return signal;
    }
    std::vector<double> smoothedSignal(n, 0.0);
    std::vector<double> weights = gaussianWeights(degree);

    // Przesuwne okno wygładzające z wagami Gaussa
    for (int i = degree; i < n - degree; ++i) {
        double sum = 0.0;
        for (int j = -degree; j <= degree; ++j) {
            sum += signal[i + j] * weights[j + degree];
        }
        smoothedSignal[i] = sum;
    }

    // Utrzymanie wartości brzegowych
    for (int i = 0; i < degree; ++i) {
        smoothedSignal[i] = signal[i];
        smoothedSignal[n - 1 - i] = signal[n - 1 - i];
    }

    return smoothedSignal;
}

// Struktura reprezentująca piksel
struct Pixel {
    unsigned char r, g, b;
};

// Funkcja wczytująca obraz PPM===============================================================================================


std::tuple<std::vector<std::vector<Pixel>>, int, int> loadPPM(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Nie udało się otworzyć pliku!" << std::endl;
        exit(0);
    }

    std::string format;
    int width, height;
    file >> format >> width >> height;
    int maxVal;
    file >> maxVal;
    file.ignore(); // Pominięcie znaku nowej linii

    std::vector<std::vector<Pixel>> image(height, std::vector<Pixel>(width));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            file.read(reinterpret_cast<char*>(&image[i][j]), sizeof(Pixel));
        }
    }

    return {image, width, height};
}



std::vector<std::vector<Pixel>> filterImage(const std::vector<std::vector<Pixel>>& image, int width, int height, int filtrationLevel) {
    std::vector<std::vector<Pixel>> filteredImage(height, std::vector<Pixel>(width));

    for (int i = 0; i < height - filtrationLevel; ++i) {
        for (int j = 0; j < width - filtrationLevel; ++j) {
        
                
            int sumR = 0, sumG = 0, sumB = 0;
            int count = 0;

            for (int fi = -filtrationLevel; fi <= filtrationLevel; ++fi) {
                for (int fj = -filtrationLevel; fj <= filtrationLevel; ++fj) {
                    int ni = i + fi;
                    int nj = j + fj;
                    if (ni >= 0 && ni < height && nj >= 0 && nj < width) {
                        sumR += image[ni][nj].r;
                        sumG += image[ni][nj].g;
                        sumB += image[ni][nj].b;
                        count++;
                    }
                }
            }

            filteredImage[i][j].r = sumR / count;
            filteredImage[i][j].g = sumG / count;
            filteredImage[i][j].b = sumB / count;
        }
    }

    // Copy border pixels from the original image to the filtered image
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            if (i < filtrationLevel || i >= height - filtrationLevel ||
                j < filtrationLevel || j >= width - filtrationLevel) {
                filteredImage[i][j] = image[i][j];
            }
        }
    }

    return filteredImage;
}

// Funkcja zapisująca obraz PPM
void savePPM(const std::string& filename, const std::vector<std::vector<Pixel>>& image, int width, int height) {
    std::ofstream file(filename, std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";

    for (const auto& row : image) {
        for (const auto& pixel : row) {
            file.write(reinterpret_cast<const char*>(&pixel), sizeof(Pixel));
        }
    }
}

// Funkcja wyświetlająca obraz w Matplot++
void displayImageMatplot(const std::vector<std::vector<Pixel>>& image, int width, int height, const std::string& window_title = "Obraz") {
    std::vector<std::vector<double>> grayscaleImage(height, std::vector<double>(width));

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            grayscaleImage[i][j] = (image[i][j].r + image[i][j].g + image[i][j].b) / 3.0;
        }
    }

    auto fig = figure();
    imagesc(grayscaleImage);
    title(window_title);
    show();
}

// Funkcja filtrująca obraz 2D
void filter2D(const std::string& input_filename, const std::string& output_filename, int filtrationLevel) {
    auto [image, width, height] = loadPPM(input_filename);
    auto filteredImage = filterImage(image, width, height, filtrationLevel);
    savePPM(output_filename, filteredImage, width, height);
    displayImageMatplot(filteredImage, width, height, "Filtered Image");
}

//=====================================================================================================================================================

// === Rejestracja w module _core ===
PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Moduł C++ do generowania, analizowania (DFT/IDFT) 
        i wizualizacji sygnałów z pybind11 + matplot++
    )pbdoc";

    // Expose Pixel struct to Python
    py::class_<Pixel>(m, "Pixel")
        .def(py::init<>())
        .def_readwrite("r", &Pixel::r)
        .def_readwrite("g", &Pixel::g)
        .def_readwrite("b", &Pixel::b);

    // podstawowe
    m.def("add", &add, "Dodaje dwie liczby całkowite");
    m.def("subtract", &subtract, "Odejmuje dwie liczby całkowite");
    m.def("prep_x", &prep_x, py::arg("start"), py::arg("end"), py::arg("n_samples"),
          "Przygotowuje wektor x z wartościami od start do end z n_samples próbkami");
    m.def("noise", &noise, py::arg("signal"), py::arg("noise_level"),
          "Generuje szum dla sygnału z określonym poziomem szumu");
    m.def("Amplitude", &Amplitude, py::arg("signal"),
          "Oblicza amplitudę sygnału");
    m.def("download", &download, py::arg("name") = "",
          "Zapisuje wykres do pliku o nazwie name (bez rozszerzenia) w formacie PNG");
    // generatory
    m.def("gen_sin", &gen_sin, py::arg("x"), py::arg("freq"), py::arg("phase"), py::arg("shift_y"), py::arg("amplitude"),
          "Generuje wartości sinusa z wektora x i parametrami");    
    m.def("gen_cos", &gen_cos, py::arg("x"), py::arg("freq"), py::arg("phase"), py::arg("shift_y"), py::arg("amplitude"),
          "Generuje wartości cosinusa z wektora x i parametrami");
    m.def("gen_parabola", &gen_parabola, py::arg("x"), py::arg("a"), py::arg("p"), py::arg("q"),
          "Generuje wartości paraboli z wektora x i parametrami: y = a*(x - p)^2 + q");
    m.def("gen_sqr", &gen_sqr, py::arg("x"), py::arg("fill"), py::arg("freq"), py::arg("low"), py::arg("pick"),
          "Generuje wartości sygnału prostokątnego z wektora x i parametrami: wypełnienie, częstotliwość, niskie i wysokie wartości");
    m.def("gen_saw", &gen_saw, py::arg("x"), py::arg("freq"), py::arg("grade"), py::arg("low"), py::arg("pick"),
          "Generuje wartości sygnału piłokształtnego z wektora x i parametrami: częstotliwość, stopień narastania, niskie i wysokie wartości");
    
    // wykresy sygnałów
    m.def("plot", &plot, py::arg("x"), py::arg("y"),
          "Rysuje wykres sygnału z wektorów x i y");    



    // DFT plotting functions
    m.def("DFT", &DFT, py::arg("signal"), "Oblicza DFT");
     m.def("plot_dft", &plot_dft, py::arg("signal"), py::arg("sampling_rate"), py::arg("name") = "",
          "Rysuje widmo amplitudowe sygnału (DFT)");

         m.def("idft", &idft, py::arg("spectrum"), "Oblicza IDFT");

    m.def("filter_signal", &filter_signal, py::arg("signal"), py::arg("degree"),
          "Filtruje sygnał przez uśrednianie z n_samples próbkami");
    m.def("gaussianWeights", &gaussianWeights, py::arg("degree"),
          "Generuje wagi Gaussa dla filtru B-spline o stopniu degree");
    // filtracja 2D
    m.def("filter2D", &filter2D, py::arg("input_filename"), py::arg("output_filename"), py::arg("filtrationLevel"),
          "Filtruje obraz 2D z pliku PPM i zapisuje wynik do innego pliku PPM");
           
  m.def("loadPPM",
      static_cast<std::tuple<std::vector<std::vector<Pixel>>, int, int>(*)(const std::string&)>(&loadPPM),
      py::arg("filename"),
      "Wczytuje obraz PPM z pliku i zwraca (obraz, szerokość, wysokość)");

   
    m.def("savePPM", &savePPM, py::arg("filename"), py::arg("image"), py::arg("width"), py::arg("height"),
          "Zapisuje wektor pikseli do pliku PPM");
   m.def("displayImageMatplot", &displayImageMatplot, 
      py::arg("image"), py::arg("width"), py::arg("height"), py::arg("window_title") = "Obraz",
      "Wyświetla obraz w Matplot++ z wektora pikseli");
    
#ifdef VERSION_INFO
    m.attr("version") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("version") = "dev";
#endif
}
