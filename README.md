Cloth sim using SDL2

Builds on Windows with mingw-w64 using:
g++ main.cpp -o build/main.exe -std=c++17 -mwindows -lmingw32 -lSDL2main -lSDL2 -O3

SSE is not efficiently implemented, uncomment '// #define use_SSE' in main.cpp
Also need to add '-mavx' build argument to g++.