# Mandelbrot Set Renderer (MPI + OpenMP)
This project implements a high-performance Mandelbrot Set generator using hybrid parallelism with MPI and OpenMP. It visualizes the Mandelbrot fractal by distributing image row computations across multiple compute nodes (MPI ranks) and parallelizing pixel calculations within each node using OpenMP threads.

## 🚀 Key Features
- Master-worker architecture using MPI to coordinate distributed work
- Dynamic task assignment for adaptive load balancing
- Fine-grained multithreading using OpenMP
- Dynamic scheduling to reduce thread idle time
- Final image generation using custom image rendering functions

## 🖥️ Technologies
- C
- MPI (Message Passing Interface)
- OpenMP
- Custom image output via `matToImageColor`

## ⚙️ How It Works
- The master process (rank 0) dynamically assigns chunks of image rows to worker processes.
- Each worker uses OpenMP to parallelize the computation of each pixel in its assigned rows.
- Workers return computed image data back to the master, which assembles the final image.
- The resulting Mandelbrot fractal is saved as a `.jpg` file.

## 🧮 Performance
- Supports tuning of:
  - MPI work chunk sizes
  - OpenMP thread scheduling strategies
- Performance improves with increased thread count and properly balanced chunk sizes.
- Designed to scale across multiple nodes with 6 threads per node.

## 📸 Output
The result is a colored Mandelbrot fractal image (`mandelbrot.jpg`) based on a pixel-by-pixel iteration count across a complex number grid.

## 📁 File Structure
- `main.c` — core computation and parallel logic
- `matToImage.c` / `.h` — external image output functions (provided)
- `Makefile` — compilation script

## 🧪 Example Usage
```bash
mpicc -fopenmp -o mandelbrot main.c matToImage.c -lm
mpirun -np 6 ./mandelbrot

