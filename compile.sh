#!/bash/bin
g++ -O3 -Wall -shared -std=c++17 -fPIC -I/path/to/pybind11/include `python3 -m pybind11 --includes` line.cpp binding.cpp -o binding`python3-config --extension-suffix`

