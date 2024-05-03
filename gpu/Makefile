ARGS = -Minfo=all -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/ -lboost_program_options

all:
	make gpu
gpu:
	pgc++ -acc=gpu $(ARGS) -o $@ gpu.cpp

clean:all
	rm gpu