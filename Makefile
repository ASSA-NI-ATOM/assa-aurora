# =============================================================
# ==  Canonical Makefile for ASSA-NI Aurora Sieve            ==
# =============================================================
TARGET=aurora
SOURCE=aurora.cu

# Используем g++ как хост-компилятор
CXX=g++
NVCC=nvcc

# Флаги C++ (включая OpenMP для CPU-части)
CXXFLAGS=-O3 -fopenmp -std=c++17
# Флаги CUDA (только оптимизация и архитектура)
NVCCFLAGS=-O3 -arch=sm_89 # <-- ЗАМЕНИТЕ 89 НА ВАШУ АРХИТЕКТУРУ

all: ${TARGET}

${TARGET}: ${SOURCE}
	${NVCC} ${NVCCFLAGS} -o $@ $< -Xcompiler "${CXXFLAGS}"

# Команда для быстрой проверки работоспособности
check: ${TARGET}
	./${TARGET} 134217728

# Команда для очистки
clean:
	rm -f ${TARGET} *.o

.PHONY: all check clean