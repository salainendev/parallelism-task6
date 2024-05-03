ARGS = -Minfo=all -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/ -lboost_program_options

all:
	make serial
	make multi

serial:
	pgc++ -acc=host $(ARGS) -o $@ cpu.cpp

multi:
	pgc++ -acc=multicore $(ARGS) -o $@ cpu.cpp

clean:all
	rm serial multi