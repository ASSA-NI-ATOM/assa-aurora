cmake_minimum_required(VERSION 3.19)
project(aurora_arxiv_final LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 89)

find_package(OpenMP REQUIRED)

# ИСПРАВЛЕНО: Имя цели теперь едино
add_executable(aurora_arxiv_final aurora_arxiv_final.cu)

target_compile_options(aurora_arxiv_final PRIVATE -O3)
target_link_libraries(aurora_arxiv_final PRIVATE OpenMP::OpenMP_CXX)