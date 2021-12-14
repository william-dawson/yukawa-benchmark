CC=g++
CFLAGS=-O3
INCLUDE=-I$(CONDA_PREFIX)/include/ -I$(CONDA_PREFIX)/include/eigen3
LIBS=-lint2

prog: main.cc
	$(CC) $^ $(INCLUDE) -o $@ $(LIBS)

test:
	./prog input/H2O.xyz cc-pVDZ 1.0 1e-6 output/H2O.mtx
	python viz.py output/H2O.mtx output/H2O.png

clean:
	rm prog
