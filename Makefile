cpu:
	pgc++ -acc=host -Minfo=accel -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/ -o task6Host task6.cpp
gpu:
	pgc++ -acc=gpu -Minfo=accel -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/ -o task6Gpu task6.cpp
multi:
	pgc++ -acc=multicore -Minfo=accel -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/ -o task6Multi task6.cpp