ARGS = -Minfo=all -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/ -lboost_program_options

all:
	make cpu
	make gpu
	make multi

cpu:
	pgc++ -acc=host $(ARGS) -o $@ task6.cpp
gpu:
	pgc++ -acc=gpu $(ARGS) -o $@ task6.cpp
multi:
	pgc++ -acc=multicore $(ARGS) -o $@ task6.cpp

clean:all
	rm cpu gpu multi