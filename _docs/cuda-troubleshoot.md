---
title: CUDA Troubleshooting
---

## Compiling Issues

#### `make: nvcc: Command not found` error when I compile with `make cuda=1`

Make sure that the NVIDIA CUDA toolkit is installed. See instruction at the [official installation guide](https://docs.nvidia.com/cuda/).
If you still get this error after the toolkit installation, then `nvcc` is probably not in your PATH. In that case, either add the nvcc location to your PATH or manually specify the nvcc location through a Makefile variable.

Example :

If you installed the CUDA toolkit through `apt` in Ubuntu,
```sh
make cuda=1 NVCC=/usr/local/cuda/bin/nvcc
```

If you did the [Runfile installation](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#runfile) on Ubuntu,
```sh
make cuda=1 NVCC=/usr/local/cuda-<toolkit-version>/bin/nvcc
```
Note that the location of `nvcc` might be different depending on your distribution and the installation method.

#### Cannot find `-lcudart_static` error.

The default CUDA library path in the Makefile is set to be `/usr/local/cuda/lib64`.

While this is the default path for an Ubuntu 64-bit system with the CUDA toolkit installed using the package manager `apt`, it might be different on your system. You can manually specify the path to the cuda library when compiling.

Example:

If you did the [Runfile installation](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#runfile) on Ubuntu,
```sh
make cuda=1 CUDA_LIB=/usr/local/cuda-<toolkit-version>/lib64/
```

If you are using Ubuntu 32-bit,
```sh
make cuda=1 CUDA_LIB=/usr/local/cuda/lib/
```

Note that the location of the CUDA library path might be different depending on your distribution and the installation method.

#### memcpy was not declared in this scope error

If you get an error like this:

```console
$ make cuda=1
nvcc -x cu -g -O2 -std=c++11 -lineinfo -Xcompiler -Wall -I./htslib -DHAVE_CUDA=1 -rdc=true -c src/f5c.cu -o build/f5c_cuda.o
/usr/include/string.h: In function ‘void* __mempcpy_inline(void*, const void*, size_t)’:
/usr/include/string.h:652:42: error: ‘memcpy’ was not declared in this scope
return (char *) memcpy (__dest, __src, __n) + __n;
^
Makefile:76: recipe for target 'build/f5c_cuda.o' failed
make: *** [build/f5c_cuda.o] Error 1
```

Compile with `-D_FORCE_INLINES` appended to `CUDA_CFLAGS` when calling `make`:

```sh
CUDA_CFLAGS+="-D_FORCE_INLINES" make cuda=1
```

The issue is reported in [#36](https://github.com/hasindu2008/f5c/issues/36). And this fix is suggested in [opencv/opencv#6500](https://github.com/opencv/opencv/issues/6500).

## Runtime Errors

#### CUDA driver version is insufficient for CUDA runtime version

Check the following in order:
1. Do you have an NVIDIA GPU / is your NVIDIA GPU recognised by the system?
On most distributions you can use the following command to verify:
```sh
lspci | grep -i "vga\|3d\|display"
```
It should list the NVIDIA GPU:
```
01:00.0 3D controller: NVIDIA Corporation GP107M [GeForce GTX 1050 Ti Mobile] (rev a1)
```
2. Have you installed the NVIDIA driver (not the open source nouveau driver)? On most distributions you can check your graphics card driver using
```sh
lspci -nnk | grep -iA2 "vga\|3d\|display"
```
If the kernel driver output contains `nvidia`, then you are using the correct driver.
```
01:00.0 3D controller [0302]: NVIDIA Corporation GP107M [GeForce GTX 1050 Ti Mobile] [10de:1c8c] (rev a1)
       Kernel driver in use: nvidia
       Kernel modules: nvidiafb, nouveau, nvidia_396, nvidia_396_drm
```
3. If you are using a Tegra GPU (e.g. Jetson TX2), does the current user belong to the "video" user group?
Check the current group names with `group [user]`.
4. Is the CUDA driver version too old for the toolkit that is used to compile with? See <https://docs.nvidia.com/deploy/cuda-compatibility/index.html#binary-compatibility>.
The release CUDA binary that we provide is compiled using CUDA toolkit 6.5. The CUDA runtime library is statically linked and therefore release CUDA binaries work on driver version >= 340.21.
Use `nvidia-smi` to check your driver version.
```console
$ nvidia-smi
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 396.44                 Driver Version: 396.44                    |
```
If you compiled the binary yourself, see <https://docs.nvidia.com/deploy/cuda-compatibility/index.html#binary-compatibility> to check if your toolkit version and driver version match.

## Cuda error: named symbol not found 

If you get an error like `[gpu_assert::ERROR] Cuda error: named symbol not found`, this is likely to be that your GPU is old and it is deprecated in nvcc. First find the compute capability for your GPU from [here](https://developer.nvidia.com/cuda-gpus). For example, if my GPU is Quadro K620, compute capability is 5.0. Now compile f5c for your GPU architecture by passing the `CUDA_ARCH=-arch=sm_xy` to the make file. In our example, as the compute capability is 5.0. sm_xy is sm_50:

```
make cuda=1 CUDA_ARCH=-arch=sm_50
```
Remeber to change sm_50 to yours.

