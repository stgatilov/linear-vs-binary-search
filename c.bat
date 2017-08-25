cl search.cpp /O2 /W2 /EHsc /D _CRT_SECURE_NO_DEPRECATE /D NDEBUG /arch:AVX2 /FAs
rem g++ search.cpp --std=c++11 -fno-strict-overflow -O3 -D NDEBUG -mavx2 -o search_gcc.exe
