#!bin/bash
g++ -pg -std=gnu++11 main_mb.cpp minband_*.cpp -I/c/boost_1_58_0 -ID:/Users/jbenade/Documents/armadillo-6.200.4/include -ID:/Users/jbenade/Documents/eigen-eigen/ -O3 -o main_mb.exe
