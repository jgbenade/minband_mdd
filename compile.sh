#!bin/bash
g++ -pg -std=gnu++11 main_mb.cpp minband_*.cpp -I/c/boost_1_58_0 -O3 -o main_mb.exe
