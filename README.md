Cloth sim using SDL2 on Windows.

![cloth.png](https://github.com/doug-h/cloth_sim/blob/master/cloth.png)

Builds on with mingw-w64 using the following arguments:
```
g++ main.cpp  -std=c++17 -mwindows -lmingw32 -lSDL2main -lSDL2
```


SSE is not efficiently implemented, to use it uncomment '// #define use_SSE' in main.cpp.

Also need to add '-mavx' build argument to g++.
