cpu:
	pgc++ -acc=host -Minfo=accel -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/ -o task6_host task6.cpp
gpu:
	pgc++ -acc=gpu -Minfo=accel -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/ -o task6_gpu task6.cpp