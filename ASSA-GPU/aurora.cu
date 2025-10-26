// ===================================================================
// ==   ASSA-NI-ATOM: "Aurora" (Final arXiv Version)                ==
// ===================================================================
// vFinal (SRE Audit Passed):

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdint>
#include <chrono>
#include <stdexcept>
#include <algorithm>
#include <omp.h>
#include <cuda_runtime.h>
#include <cmath> // для sqrt

using u64 = uint64_t;
using u32 = uint32_t;

#define CUDA_CHECK(call)                                             \
    do {                                                             \
        cudaError_t err = call;                                      \
        if (err != cudaSuccess) {                                    \
            fprintf(stderr, "FATAL CUDA Error in %s at line %d: %s\n",\
                    __FILE__, __LINE__, cudaGetErrorString(err));    \
            exit(EXIT_FAILURE);                                      \
        }                                                            \
    } while (0)

// ИСПРАВЛЕНО: Кросс-платформенный макрос для POPCOUNT
#if defined(_MSC_VER)
#include <intrin.h>
#define POPCOUNT __popcnt64
#else
#define POPCOUNT __builtin_popcountll
#endif

// Прототипы
std::vector<u32> buildFilters(u64 limit);
void gpuSieve(u64 n, u64 n2, u64 maxX, const std::vector<u32>& h_filters, u64& total_primes);

// Ядро GPU с корректной логикой для СЕГМЕНТОВ
__global__ void sieve_kernel_segmented(
    u64 n, u64 n2,
    const u32* __restrict__ d_filters,
    u32 num_filters,
    u64 seg_base_x, u64 seg_size_x,
    u64* __restrict__ d_sieve_mask)
{
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_filters) return;

    u32 d = d_filters[tid];
    if (d == 0) return;
    u64 r = n2 % d;
    
    u64 first_x_in_seg;
    if (r >= seg_base_x) {
        first_x_in_seg = r;
    } else {
        // ИСПРАВЛЕНО: Защита от переполнения
        u64 diff = seg_base_x - r;
        u64 steps = (diff + d - 1) / d;
        if (d > 1 && steps > (UINT64_MAX - r) / d) return; 
        first_x_in_seg = r + steps * d;
    }
    if ((first_x_in_seg & 1) != 0) {
        first_x_in_seg += d;
    }

    for (u64 x = first_x_in_seg; x < seg_base_x + seg_size_x; x += 2 * d) {
        u64 relative_x = x - seg_base_x;
        u64 bit_index = relative_x / 2;
        
        if (bit_index < (seg_size_x / 2)) {
            atomicOr(reinterpret_cast<unsigned long long*>(&d_sieve_mask[bit_index / 64]), (1ULL << (bit_index % 64)));
        }
    }
}

int main(int argc, char *argv[]) {
    #ifdef _MSC_VER
        system("chcp 65001 > nul");
    #endif
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <n>\n";
        return 1;
    }
    
    auto t0 = std::chrono::high_resolution_clock::now();
    try {
        u64 n = std::stoull(argv[1]);
        if (n < 3 || n > 2147483647ULL) {
            std::cerr << "Error: n must be in [3, 2147483647]\n";
            return 1;
        }
        if ((n & 1) == 0) ++n;

        const u64 n2 = n * n;
        const u64 maxX = 4 * n - 4;
        
        std::cout << "==========================================\n";
        std::cout << "Aurora Sieve (arXiv Ready)   n = " << n << "\n";
        std::cout << "==========================================\n";

        auto h_filters = buildFilters(n);
        u64 total_primes = 0;
        gpuSieve(n, n2, maxX, h_filters, total_primes);
        
        auto t1 = std::chrono::high_resolution_clock::now();
        double sec = std::chrono::duration<double>(t1 - t0).count();
        
        std::cout << "------------------------------------------\n";
        std::cout << "Primes found : " << total_primes << "\n";
        std::cout << "Time         : " << std::fixed << std::setprecision(4) << sec << " s\n";
        std::cout << "==========================================\n";

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

std::vector<u32> buildFilters(u64 limit) {
    std::cout << "Generating filters up to n=" << limit << "..." << std::flush;
    auto t_start = std::chrono::high_resolution_clock::now();
    std::vector<char> raw(limit + 1, 1);
    raw[0] = raw[1] = 0;
    
    // ИСПРАВЛЕНО: Канонический OpenMP цикл
    u64 sqrt_limit = static_cast<u64>(std::sqrt((long double)limit)) + 1;
    #pragma omp parallel for schedule(static)
    for (u64 p = 2; p <= sqrt_limit; ++p) {
        if (raw[p]) {
            for (u64 i = p * p; i <= limit; i += p) raw[i] = 0;
        }
    }

    std::vector<u32> out;
    out.reserve(limit / 10);
    for (u64 p = 3; p <= limit; p += 2) if (raw[p]) out.push_back(p);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << " Done (" << std::fixed << std::setprecision(2) 
              << std::chrono::duration<double>(t_end - t_start).count() << "s, " 
              << out.size() << " filters)\n";
    return out;
}

void gpuSieve(u64 n, u64 n2, u64 maxX, const std::vector<u32>& h_filters, u64& total_primes) {
    u32* d_filters_ptr;
    CUDA_CHECK(cudaMalloc(&d_filters_ptr, h_filters.size() * sizeof(u32)));
    CUDA_CHECK(cudaMemcpy(d_filters_ptr, h_filters.data(), h_filters.size() * sizeof(u32), cudaMemcpyHostToDevice));

    const u64 SEG_SIZE_X = 1ULL << 24;
    const u64 SEG_WORDS = (SEG_SIZE_X / 2 + 63) / 64;
    u64* d_seg_mask;
    CUDA_CHECK(cudaMalloc(&d_seg_mask, SEG_WORDS * sizeof(u64)));
    
    total_primes = 0;

    for (u64 seg_base_x = 2; seg_base_x <= maxX; seg_base_x += SEG_SIZE_X) {
        u64 current_seg_size_x = std::min(SEG_SIZE_X, maxX - seg_base_x + 2);
        u64 num_candidates_in_seg = current_seg_size_x / 2;
        u64 current_seg_words = (num_candidates_in_seg + 63) / 64;
        
        CUDA_CHECK(cudaMemset(d_seg_mask, 0, current_seg_words * sizeof(u64)));

        int threads = 256;
        int blocks = (h_filters.size() + threads - 1) / threads;

        if (blocks > 0) {
            sieve_kernel_segmented<<<blocks, threads>>>(
                n, n2, d_filters_ptr, h_filters.size(),
                seg_base_x, current_seg_size_x,
                d_seg_mask
            );
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        std::vector<u64> h_seg(current_seg_words);
        CUDA_CHECK(cudaMemcpy(h_seg.data(), d_seg_mask, current_seg_words * sizeof(u64), cudaMemcpyDeviceToHost));
        
        u64 seg_primes = 0;
        #pragma omp parallel for reduction(+:seg_primes)
        for(int64_t i = 0; i < static_cast<int64_t>(current_seg_words); ++i) {
             if (i == static_cast<int64_t>(current_seg_words) - 1) {
                // ИСПРАВЛЕНО: Безопасная маска для последнего слова
                u64 last_mask = (num_candidates_in_seg % 64 == 0) ? ~0ULL : ((1ULL << (num_candidates_in_seg % 64)) - 1);
                seg_primes += POPCOUNT(~h_seg[i] & last_mask);
            } else {
                seg_primes += POPCOUNT(~h_seg[i]);
            }
        }
        total_primes += seg_primes;
    }

    CUDA_CHECK(cudaFree(d_filters_ptr));
    CUDA_CHECK(cudaFree(d_seg_mask));
}
