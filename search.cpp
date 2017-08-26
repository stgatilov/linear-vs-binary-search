#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdint.h>
#include <algorithm>
#include <random>
#include <immintrin.h>
#include <intrin.h>

//little-endian _MM_SHUFFLE
#define SHUF(i0, i1, i2, i3) (i0 + i1*4 + i2*16 + i3*64)

#ifdef _MSC_VER
    #define FORCEINLINE __forceinline
    #define NOINLINE __declspec(noinline)
    #define ALIGN(n) __declspec(align(n))
    FORCEINLINE uint32_t bsr(uint32_t x) {
        unsigned long res;
        _BitScanReverse(&res, x);
        return res;
    }
    FORCEINLINE uint32_t bsf(uint32_t x) {
        unsigned long res;
        _BitScanForward(&res, x);
        return res;
    }
#else
    #define FORCEINLINE __attribute__((always_inline)) inline
    #define NOINLINE __attribute__((noinline))
    #define ALIGN(n) __attribute__((aligned(n)))
    FORCEINLINE uint32_t bsr(uint32_t x) {
        return 31 - __builtin_clz(x);
    }
    FORCEINLINE uint32_t bsf(uint32_t x) {
        return __builtin_ctz(x);
    }
#endif

//if true, then average latency of one search is measured
//if false, then average throughput performance of one search is measured
//implementation-wise, setting it makes the index of the next array/key
//to be searched dependent on the answer of the current search
#define MEASURE_LATENCY false

//controls inlining of all search functions being tested
//NOINLINE means that inlining is forbidden:
//in this case benchmarking code is less likely to influence search performance
#define TESTINLINE NOINLINE

//======================= search implementations =======================

static TESTINLINE int binary_search_std (const int *arr, int n, int key) {
    return std::lower_bound(arr, arr + n, key) - arr;
}

static TESTINLINE int binary_search_simple (const int *arr, int n, int key) {
    intptr_t left = -1;
    intptr_t right = n;
    while (right - left > 1) {
        intptr_t middle = (left + right) >> 1;
        if (arr[middle] < key)
            left = middle;
        else
            right = middle;
    }
    return right;
}

intptr_t MINUS_ONE = -1;
static TESTINLINE int binary_search_branchless (const int *arr, int n, int key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 0) {
        pos = (arr[pos + step] < key ? pos + step : pos);
        step >>= 1;
    }
    return pos + 1;
}

template<intptr_t MAXN> static TESTINLINE int binary_search_branchless_UR (const int *arr, int n, int key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    assert(n+1 == MAXN);
    assert(n < (1<<20));

    //intptr_t pos = -1;
    intptr_t pos = MINUS_ONE;
    #define STEP(logstep) \
        if ((1<<logstep) < MAXN) pos = (arr[pos + (1<<logstep)] < key ? pos + (1<<logstep) : pos);
    STEP(19)
    STEP(18)
    STEP(17)
    STEP(16)
    STEP(15)
    STEP(14)
    STEP(13)
    STEP(12)
    STEP(11)
    STEP(10)
    STEP(9)
    STEP(8)
    STEP(7)
    STEP(6)
    STEP(5)
    STEP(4)
    STEP(3)
    STEP(2)
    STEP(1)
    STEP(0)
    #undef STEP
    return pos + 1;
}

static TESTINLINE int linear_search_scalar (const int *arr, int n, int key) {
    int cnt = 0;
    for (int i = 0; i < n; i++)
        cnt += (arr[i] < key);
    return cnt;
}

static TESTINLINE int linear_search_sse (const int *arr, int n, int key) {
    assert(size_t(arr) % 16 == 0);
    __m128i vkey = _mm_set1_epi32(key);
    __m128i cnt = _mm_setzero_si128();
    for (int i = 0; i < n; i += 16) {
        __m128i mask0 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[i+0]), vkey);
        __m128i mask1 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[i+4]), vkey);
        __m128i mask2 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[i+8]), vkey);
        __m128i mask3 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[i+12]), vkey);
        __m128i sum = _mm_add_epi32(_mm_add_epi32(mask0, mask1), _mm_add_epi32(mask2, mask3));
        cnt = _mm_sub_epi32(cnt, sum);
    }
    cnt = _mm_add_epi32(cnt, _mm_shuffle_epi32(cnt, SHUF(2, 3, 0, 1)));
    cnt = _mm_add_epi32(cnt, _mm_shuffle_epi32(cnt, SHUF(1, 0, 3, 2)));
    return _mm_cvtsi128_si32(cnt);
}

template<intptr_t MAXN> static TESTINLINE int linear_search_sse_UR (const int *arr, int n, int key) {
    assert(size_t(arr) % 16 == 0);
    assert(n <= 1024);
    __m128i vkey = _mm_set1_epi32(key);
    __m128i cnt = _mm_setzero_si128();
    intptr_t i = 0;
    #define STEP \
    if (i < MAXN) {\
        __m128i mask0 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[i+0]), vkey); \
        __m128i mask1 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[i+4]), vkey); \
        __m128i mask2 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[i+8]), vkey); \
        __m128i mask3 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[i+12]), vkey); \
        __m128i sum = _mm_add_epi32(_mm_add_epi32(mask0, mask1), _mm_add_epi32(mask2, mask3)); \
        cnt = _mm_sub_epi32(cnt, sum); \
    } i += 16;
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    #undef STEP
    cnt = _mm_add_epi32(cnt, _mm_shuffle_epi32(cnt, SHUF(2, 3, 0, 1)));
    cnt = _mm_add_epi32(cnt, _mm_shuffle_epi32(cnt, SHUF(1, 0, 3, 2)));
    return _mm_cvtsi128_si32(cnt);
}

static TESTINLINE int linear_search_avx (const int *arr, int n, int key) {
    assert(size_t(arr) % 32 == 0);
    __m256i vkey = _mm256_set1_epi32(key);
    __m256i cnt = _mm256_setzero_si256();
    for (int i = 0; i < n; i += 16) {
        __m256i mask0 = _mm256_cmpgt_epi32(vkey, _mm256_load_si256((__m256i *)&arr[i+0]));
        __m256i mask1 = _mm256_cmpgt_epi32(vkey, _mm256_load_si256((__m256i *)&arr[i+8]));
        __m256i sum = _mm256_add_epi32(mask0, mask1);
        cnt = _mm256_sub_epi32(cnt, sum);
    }
    __m128i xcnt = _mm_add_epi32(_mm256_extracti128_si256(cnt, 1), _mm256_castsi256_si128(cnt));
    xcnt = _mm_add_epi32(xcnt, _mm_shuffle_epi32(xcnt, SHUF(2, 3, 0, 1)));
    xcnt = _mm_add_epi32(xcnt, _mm_shuffle_epi32(xcnt, SHUF(1, 0, 3, 2)));
    return _mm_cvtsi128_si32(xcnt);
}

template<intptr_t MAXN> static TESTINLINE int linear_search_avx_UR (const int *arr, int n, int key) {
    assert(size_t(arr) % 32 == 0);
    assert(n <= 1024);
    __m256i vkey = _mm256_set1_epi32(key);
    __m256i cnt = _mm256_setzero_si256();
    intptr_t i = 0;
    #define STEP \
    if (i < MAXN) {\
        __m256i mask0 = _mm256_cmpgt_epi32(vkey, _mm256_load_si256((__m256i *)&arr[i+0])); \
        __m256i mask1 = _mm256_cmpgt_epi32(vkey, _mm256_load_si256((__m256i *)&arr[i+8])); \
        __m256i sum = _mm256_add_epi32(mask0, mask1); \
        cnt = _mm256_sub_epi32(cnt, sum); \
    } i += 16;
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    STEP STEP STEP STEP STEP STEP STEP STEP
    #undef STEP
    __m128i xcnt = _mm_add_epi32(_mm256_extracti128_si256(cnt, 1), _mm256_castsi256_si128(cnt));
    xcnt = _mm_add_epi32(xcnt, _mm_shuffle_epi32(xcnt, SHUF(2, 3, 0, 1)));
    xcnt = _mm_add_epi32(xcnt, _mm_shuffle_epi32(xcnt, SHUF(1, 0, 3, 2)));
    return _mm_cvtsi128_si32(xcnt);
}

// from https://stackoverflow.com/questions/2741859/how-fast-can-you-make-linear-search
// type of "i" changed from int to intptr_t
static TESTINLINE int linearX_search_scalar (const int *arr, int n, int key) {
    intptr_t i = 0;
    while (i < n) {
        if (arr[i] >= key)
            break;
        ++i;
    }
    return i;
}

// from https://schani.wordpress.com/2010/04/30/linear-vs-binary-search/
// type of "i" changed from int to intptr_t
static TESTINLINE int linearX_search_sse (const int *arr, int n, int key) {
    __m128i *in_data = (__m128i*)arr;
    __m128i key4 = _mm_set1_epi32(key);
    intptr_t i = 0;
    int res;
    for (;;) {
        __m128i cmp0 = _mm_cmpgt_epi32 (key4, in_data [i + 0]);
        __m128i cmp1 = _mm_cmpgt_epi32 (key4, in_data [i + 1]);
        __m128i cmp2 = _mm_cmpgt_epi32 (key4, in_data [i + 2]);
        __m128i cmp3 = _mm_cmpgt_epi32 (key4, in_data [i + 3]);
        __m128i pack01 = _mm_packs_epi32 (cmp0, cmp1);
        __m128i pack23 = _mm_packs_epi32 (cmp2, cmp3);
        __m128i pack0123 = _mm_packs_epi16 (pack01, pack23);
        res = _mm_movemask_epi8 (pack0123);
        if (res != 0xFFFF)
            break;
        i += 4;
    }
    return i * 4 + bsf(~res);
}

//additional versions (experimental)

static TESTINLINE int hybrid_search (const int *arr, int n, int key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    intptr_t pos = MINUS_ONE;
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 8) {
        pos = (arr[pos + step] < key ? pos + step : pos);
        step >>= 1;
    }
    pos++;
    step <<= 1;

    assert(size_t(arr) % 16 == 0);
    assert(pos % 16 == 0);
    __m128i vkey = _mm_set1_epi32(key);
    __m128i mask0 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[pos+0]), vkey);
    __m128i mask1 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[pos+4]), vkey);
    __m128i mask2 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[pos+8]), vkey);
    __m128i mask3 = _mm_cmplt_epi32(_mm_load_si128((__m128i *)&arr[pos+12]), vkey);
    __m128i cnt = _mm_add_epi32(_mm_add_epi32(mask0, mask1), _mm_add_epi32(mask2, mask3));
    cnt = _mm_add_epi32(cnt, _mm_shuffle_epi32(cnt, SHUF(2, 3, 0, 1)));
    cnt = _mm_add_epi32(cnt, _mm_shuffle_epi32(cnt, SHUF(1, 0, 3, 2)));

    return pos - _mm_cvtsi128_si32(cnt);
}

static TESTINLINE int hybridX_search (const int *arr, int n, int key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    intptr_t pos = MINUS_ONE;
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 8) {
        pos = (arr[pos + step] < key ? pos + step : pos);
        step >>= 1;
    }
    pos++;
    step <<= 1;

    assert(size_t(arr) % 16 == 0);
    assert(pos % 16 == 0);
    __m128i vkey = _mm_set1_epi32(key);
    __m128i cnt = _mm_setzero_si128();
    __m128i cmp0 = _mm_cmpgt_epi32 (vkey, _mm_load_si128((__m128i *)&arr[pos+0]));
    __m128i cmp1 = _mm_cmpgt_epi32 (vkey, _mm_load_si128((__m128i *)&arr[pos+4]));
    __m128i cmp2 = _mm_cmpgt_epi32 (vkey, _mm_load_si128((__m128i *)&arr[pos+8]));
    __m128i cmp3 = _mm_cmpgt_epi32 (vkey, _mm_load_si128((__m128i *)&arr[pos+12]));
    __m128i pack01 = _mm_packs_epi32 (cmp0, cmp1);
    __m128i pack23 = _mm_packs_epi32 (cmp2, cmp3);
    __m128i pack0123 = _mm_packs_epi16 (pack01, pack23);
    uint32_t res = _mm_movemask_epi8 (pack0123);

    return pos + bsf(~res);
}

static TESTINLINE int binary_search_branchlessM (const int *arr, int n, int key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 0) {
        pos += (arr[pos + step] < key) * step;
        step >>= 1;
    }
    return pos + 1;
}

static TESTINLINE int binary_search_branchlessA (const int *arr, int n, int key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 0) {
        pos += (-(arr[pos + step] < key)) & step;
        step >>= 1;
    }
    return pos + 1;
}

static TESTINLINE int binary_search_branchlessS (const int *arr, int n, int key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 0) {
        pos += ((arr[pos + step] - key) >> 31) & step;
        step >>= 1;
    }
    return pos + 1;
}

//======================= testing code =======================

//length of each input array (including one sentinel element)
//must be power of two
static const int SIZE = 64;
//how many searches are done in benchmark
static const int TRIES = (1<<30) / SIZE;
//number of pre-generated input arrays (rotated cyclically)
static const int ARR_SAMPLES = (4<<10) / SIZE;
//number of pre-generated keys for search (rotated cyclically)
static const int KEY_SAMPLES = (4<<10);

//number of elements in every array to be searched (excluding sentinel element)
int n = SIZE - 1;

//input arrays for search (aligned for AVX)
ALIGN(32) int input[ARR_SAMPLES][SIZE];
//keys to be searched
int keys[KEY_SAMPLES];


int main() {
    //used RNG
    std::mt19937 rnd;
    std::uniform_int_distribution<int> distr(0, SIZE);

    //generate all input arrays
    for (int s = 0; s < ARR_SAMPLES; s++) {
        for (int i = 0; i < n; i++)
            input[s][i] = distr(rnd);
        std::sort(input[s], input[s] + n);
        //set sentinel element to INT_MAX
        for (int i = n; i < (n+15)/16*16; i++)
            input[s][i] = INT_MAX;
    }

    //test correctness of searches AND generate all keys to be searched
    const int ITERS = std::max(TRIES / 10, std::max(KEY_SAMPLES, ARR_SAMPLES)) + 10;
    for (int t = 0; t < ITERS; t++) {
        const int *arr = input[t % ARR_SAMPLES];
        int key = keys[t % KEY_SAMPLES] = distr(rnd);

        int res[16], sk = 0;
        res[sk++] = binary_search_std(arr, n, key);
        res[sk++] = binary_search_simple(arr, n, key);
        res[sk++] = binary_search_branchless(arr, n, key);
        res[sk++] = binary_search_branchless_UR<SIZE>(arr, n, key);
        res[sk++] = linearX_search_scalar(arr, n, key);
        res[sk++] = linear_search_scalar(arr, n, key);
        res[sk++] = linearX_search_sse(arr, n, key);
        res[sk++] = linear_search_sse(arr, n, key);
        res[sk++] = linear_search_sse_UR<SIZE>(arr, n, key);
        res[sk++] = linear_search_avx(arr, n, key);
        res[sk++] = linear_search_avx_UR<SIZE>(arr, n, key);
        //some experimental implementations:
        //res[sk++] = hybrid_search(arr, n, key);
        //res[sk++] = hybridX_search(arr, n, key);
        //res[sk++] = binary_search_branchlessM(arr, n, key);
        //res[sk++] = binary_search_branchlessA(arr, n, key);
        //res[sk++] = binary_search_branchlessS(arr, n, key);

        //program terminates if any search gives different answer
        for (int i = 1; i < sk; i++)
            if (res[i-1] != res[i]) {
                printf("ERROR: ");
                for (int j = 0; j < sk; j++)
                    printf(" %d", res[j]);
                printf("\n");
                exit(0);
            }
    }

    //print some info about current benchmark parameters
    static const int DARR = 1779033703 & (ARR_SAMPLES - 1);
    static const int DKEY = 2654435769 & (KEY_SAMPLES - 1);
    printf("Arrays: %d x %d   (+%d)\n", ARR_SAMPLES, SIZE, DARR);
    printf("Keys: %d        (+%d)\n", KEY_SAMPLES, DKEY);
    printf("Memory: %d B\n", int(sizeof(input) + sizeof(keys)));

    //benchmark environment
    //note: it could had been a function instead of macro
    //but originally I wanted to benchmark inlined code of every search
    #define TEST_SEARCH(func) \
    { \
        int start = clock(); \
        int check = 0; \
        for (int t = 0; t < TRIES; t++) { \
            int i = (t * DARR + (MEASURE_LATENCY ? check&1 : 0)) & (ARR_SAMPLES - 1); \
            int j = (t * DKEY + (MEASURE_LATENCY ? check&1 : 0)) & (KEY_SAMPLES - 1); \
            const int *arr = input[i]; \
            int key = keys[j]; \
            int res = func(arr, n, key); \
            check += res; \
        } \
        double elapsed = double(clock() - start) / CLOCKS_PER_SEC; \
        printf("%8.1lf ns : %40s   (%d)\n", 1e+9 * elapsed / TRIES, #func, check); \
    }

    //run performance benchmark and print formatted results
    TEST_SEARCH(binary_search_std);
    TEST_SEARCH(binary_search_simple);
    TEST_SEARCH(binary_search_branchless);
    TEST_SEARCH(binary_search_branchless_UR<SIZE>);

    TEST_SEARCH(linearX_search_scalar);
    TEST_SEARCH(linear_search_scalar);

    TEST_SEARCH(linearX_search_sse);
    TEST_SEARCH(linear_search_sse);
    TEST_SEARCH(linear_search_sse_UR<SIZE>);
    TEST_SEARCH(linear_search_avx);
    TEST_SEARCH(linear_search_avx_UR<SIZE>);

    //some experimental implementations:
    //TEST_SEARCH(hybrid_search);
    //TEST_SEARCH(hybridX_search);

    //TEST_SEARCH(binary_search_branchlessM);
    //TEST_SEARCH(binary_search_branchlessA);
    //TEST_SEARCH(binary_search_branchlessS);

    return 0;
}
