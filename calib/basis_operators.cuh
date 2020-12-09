// # expected kernels.

__global__ void kernel_basis_cannonical_reordering (int * idx_, int n_, int in_orientation_, int * orientation_)
{
    __shared__ int orientation [1];
    __shared__ bool sorted [1];

    orientation [0] = in_orientation_;
    sorted [0] = false;

    int pos = threadIdx. x + (blockDim. x * blockIdx. x);

    while (!sorted [0])
    {
        sorted [0] = true;

        if (pos % 2 == 0 && pos + 1 < n_)
        {
            if (idx_ [pos] == idx_ [pos+1])
            {
                // # eliminação.
                idx_ [pos] = -1;
                idx_ [pos+1] = -1;
            } else {
                if (idx_ [pos] > idx_ [pos+1])
                {
                    int aux = idx_ [pos];
                    idx_ [pos] = idx_ [pos+1];
                    idx_ [pos+1] = aux;
                    sorted [0] = false;
                    orientation [0] *= -1;
                }
            }
        }

        __syncthreads ();
        
        if (pos % 2 != 0 && pos + 1 < n_)
        {
            if (idx_ [pos] == idx_ [pos+1])
            {
                // # eliminação.
                idx_ [pos] = -1;
                idx_ [pos+1] = -1;
            } else {
                if (idx_ [pos] > idx_ [pos+1])
                {
                    int aux = idx_ [pos];
                    idx_ [pos] = idx_ [pos+1];
                    idx_ [pos+1] = aux;
                    sorted [0] = false;
                    orientation [0] *= -1;
                }
            }
        }

        __syncthreads ();
    }

    orientation_ [0] = orientation [0];
}

// # bridge operations between host code and device code.

void bridge_cannonical_reordering (calib::basis& base)
{
    int n_ = base. base_index. size ();

    int * idx_;
    int * orientation;
    int * orientation_;

    orientation = (int *) std::malloc (sizeof (int));

    orientation [0] = base. orientation;

    cudaMalloc ((void **) &idx_, n_ * sizeof (int));
    cudaMalloc ((void **) &orientation_, sizeof (int));

    cudaMemcpy (idx_, base. base_index. data (), n_ * sizeof (int), cudaMemcpyHostToDevice);

    kernel_basis_cannonical_reordering <<<ceil (n_ / (float) CALIB_N_THREADS_PER_BLOCK), CALIB_N_THREADS_PER_BLOCK>>> (idx_, n_, orientation [0], orientation_);

    cudaDeviceSynchronize ();

    cudaMemcpy (base. base_index. data (), idx_, n_ * sizeof (int), cudaMemcpyDeviceToHost);
    cudaMemcpy (orientation, orientation_, sizeof (int), cudaMemcpyDeviceToHost);

    base. orientation = orientation [0];

    cudaFree (idx_);
    cudaFree (orientation_);

    std::free (orientation);
}