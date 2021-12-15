CC=g++
CFLAGS=-O2
INCLUDE=-I$(CONDA_PREFIX)/include/ -I$(CONDA_PREFIX)/include/eigen3
LIBS=-lint2

prog: main.cc
	$(CC) $^ $(CFLAGS) $(INCLUDE) -o $@ $(LIBS)

matrices=H2O Lysozyme 6lu7
$(matrices):
	./prog input/$@.xyz def2-SVP 1.0 1e-8 output/$@.block output/$@.mtx
	python viz.py output/$@.mtx output/$@.png

clean:
	rm prog
