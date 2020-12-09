// # expected kernels.

__constant__ float fact_ [CALIB_FACTORIAL_TABLE_SIZE];

__device__ float fact (int x)
{
    return (x < CALIB_FACTORIAL_TABLE_SIZE) ? fact_ [x] : x * fact (x-1);
}

__device__ void find_in_combinations (int k, int ini, int fim, int p, int q, unsigned long long int idx, int start_index, int & relative_pos, int & global_pos, int * all_indices, int * aux_combinations, int aux_combinations_size, int & result)
{
    if (k == 0)
    {
        unsigned long long int index = aux_combinations [0];

        for (int i = 1; i < aux_combinations_size; i++)
        {
            index = pair_to_uid (index, aux_combinations [i]);
        }

        if (index == idx)
        {
            global_pos = (int) ((float) relative_pos);
            if (global_pos <= p) result = 1;
            if (global_pos > q) result = 0;
            if (global_pos > p && global_pos <= q) result = -1;
        }

        relative_pos += 1;

        return;
    }

    for (int i = start_index; i <= DEFAULT_SPACE_DIM - k; i++)
    {
        aux_combinations [aux_combinations_size++] = all_indices [i];
        find_in_combinations (k - 1, ini, fim, p, q, idx, i + 1, relative_pos, global_pos, all_indices, aux_combinations, aux_combinations_size, result);
        aux_combinations_size--;
    }
}

__device__ int find_in_combinations2 (int degree, int ini, int fim, int p, int q, unsigned long long int idx)
{
    int relative_pos = ini;
    int global_pos = ini + 1;
    int result = 0;

    int all_indices [DEFAULT_SPACE_DIM];
    int aux_combinations [DEFAULT_SPACE_DIM];

    for (int i = 1; i <= DEFAULT_SPACE_DIM; i++)
    {
        all_indices [i] = i;
    }

    find_in_combinations (degree, ini, fim, p, q, idx, 0, relative_pos, global_pos, all_indices, aux_combinations, 0, result);

    return result;
}

__device__ int device_metric_inner_value (int degree, unsigned long long int idx)
{
    int p = 5;
    int q = 8;

    int ini = 0;
    int fim = 0;

    for (int k = 0; k <= degree - 1; k++)
    {
        ini += fact (DEFAULT_SPACE_DIM) / (fact (DEFAULT_SPACE_DIM - k) * fact (k));
    }

    if (ini >= q) return 0;   // # o intervalo inicia-se após q.

    fim = ini + fact (DEFAULT_SPACE_DIM) / (fact (DEFAULT_SPACE_DIM - degree) * fact (degree));

    if (fim < p) return 1;   // # o intervalo finaliza-se antes de p.

    // # o intervalo cruza (p, q).

    return find_in_combinations2 (degree, ini, fim, p, q, idx);
}

__device__ int device_binary_search (int n, int deg, unsigned long long int uid, int * deg_, unsigned long long int * uid_)
{
    int ini = 0, end = n;
    int mid = n >> 1;

    while (ini != end)
    {
        if (deg_ [mid] == deg && uid_ [mid] == uid)
            return mid;

        if ((deg < deg_ [mid]) || (deg == deg_ [mid]) && (uid < uid_ [mid]))
            end = mid;

        if ((deg > deg_ [mid]) || (deg == deg_ [mid]) && (uid > uid_ [mid]))
            ini = mid + 1;

        mid = (ini + end) >> 1;
    }
    
    return -1;
}

__global__ void kernel_metric_sum_reduct (float * array)
{
    __shared__ float s_array [CALIB_N_THREADS_PER_BLOCK];

    int tid = threadIdx. x;
    int pos = tid + blockDim. x * blockIdx. x;

    s_array [tid] = array [pos];
    
    __syncthreads ();

    for (unsigned int shift = blockDim. x / 2; shift > 0; shift >>= 1)
    {
        if (tid < shift)
            s_array [tid] += s_array [tid + shift];
    }

    __syncthreads ();

    if (tid == 0)
    {
        array [0] = 0;
        __syncthreads ();
        atomicAdd (&array [0], s_array [0]);
    }
}

__global__ void kernel_metric_inner_prd (int n_, int m_, int * deg_mv1_, float * dir_mv1_, unsigned long long int * idx_mv1_, int * deg_mv2_, float * dir_mv2_, unsigned long long int * idx_mv2_, float * inner_result_)
{
    int pos = threadIdx. x + (blockDim. x * blockIdx. x);

    if (pos < n_)
    {
        int index_in_mv2 = device_binary_search (m_, deg_mv1_ [pos], idx_mv1_ [pos], deg_mv2_, idx_mv2_);

        if (index_in_mv2 >= 0)
        {
            int inner_value = device_metric_inner_value (deg_mv1_ [pos], idx_mv1_ [pos]);

            inner_result_ [pos] = dir_mv1_ [pos] * inner_value * dir_mv2_ [index_in_mv2];
        }
    }
}

// # bridge operations between host code and device code.

double bridge_metric_inner_prd (calib::multivector& mv1, calib::multivector& mv2)
{
    // # get dimension of each multivector.
    int n_ = mv1. elems. size ();
    int m_ = mv2. elems. size ();

    // # partition size sent to device for each time.
    int partition_ = 1024;

    int times_mv1_ = std::ceil (n_ / (float) partition_);
    int times_mv2_ = std::ceil (m_ / (float) partition_);

    // # host-side memory allocation [mv #1].
    int * deg_mv1;
    float * dir_mv1;
    unsigned long long int * idx_mv1;

    deg_mv1 = (int *) std::malloc (n_ * sizeof (int));
    dir_mv1 = (float *) std::malloc (n_ * sizeof (float));
    idx_mv1 = (unsigned long long int *) std::malloc (n_ * sizeof (unsigned long long int));

    for (int i = 0; i < n_; i++)
    {
        calib::basis base = mv1. elems [i];

        deg_mv1 [i] = base. degree ();
        dir_mv1 [i] = base. direction ();
        idx_mv1 [i] = base. unique_index ();
    }

    // # host-side memory allocation [mv #2].
    int * deg_mv2;
    float * dir_mv2;
    unsigned long long int * idx_mv2;

    deg_mv2 = (int *) std::malloc (m_ * sizeof (int));
    dir_mv2 = (float *) std::malloc (m_ * sizeof (float));
    idx_mv2 = (unsigned long long int *) std::malloc (m_ * sizeof (unsigned long long int));

    for (int i = 0; i < m_; i++)
    {
        calib::basis base = mv2. elems [i];

        deg_mv2 [i] = base. degree ();
        dir_mv2 [i] = base. direction ();
        idx_mv2 [i] = base. unique_index ();
    }

    // # host-side memory allocation [mv #3].
    float * inner_result;

    inner_result = (float *) std::malloc ((n_ + m_) * sizeof (float));

    // # device-side memory allocation [mv #1].
    int * deg_mv1_;
    float * dir_mv1_;
    unsigned long long int * idx_mv1_;

    cudaMalloc ((void **) &deg_mv1_, partition_ * sizeof (int));
    cudaMalloc ((void **) &dir_mv1_, partition_ * sizeof (float));
    cudaMalloc ((void **) &idx_mv1_, partition_ * sizeof (unsigned long long int));

    // # device-side memory allocation [mv #2].
    int * deg_mv2_;
    float * dir_mv2_;
    unsigned long long int * idx_mv2_;

    cudaMalloc ((void **) &deg_mv2_, partition_ * sizeof (int));
    cudaMalloc ((void **) &dir_mv2_, partition_ * sizeof (float));
    cudaMalloc ((void **) &idx_mv2_, partition_ * sizeof (unsigned long long int));

    // # device-side memory allocation [mv #3].
    float acc_result = 0.0f;
    float * inner_result_;

    cudaMalloc ((void **) &inner_result_, partition_ * sizeof (float));

    // # control variables.
    int shift = 0;

    for (int t_i = 0; t_i < times_mv1_; t_i++)
    {
        int step_mv1 = (t_i) * partition_;
        int step_n_elems = partition_;

        if (step_mv1 + step_n_elems > n_)
        {
            step_n_elems = n_ - step_mv1;
        }

        // # sends packet from first multivector to device.

        cudaMemcpy (deg_mv1_, deg_mv1 + step_mv1, step_n_elems * sizeof (int), cudaMemcpyHostToDevice);
        cudaMemcpy (dir_mv1_, dir_mv1 + step_mv1, step_n_elems * sizeof (float), cudaMemcpyHostToDevice);
        cudaMemcpy (idx_mv1_, idx_mv1 + step_mv1, step_n_elems * sizeof (unsigned long long int), cudaMemcpyHostToDevice);

        for (int t_j = 0; t_j < times_mv2_; t_j++)
        {
            int step_mv2 = (t_j) * partition_;
            int step_m_elems = partition_;

            if (step_mv2 + step_m_elems > m_)
            {
                step_m_elems = m_ - step_mv2;
            }

            // # sends packet from second multivector to device.
            
            cudaMemcpy (deg_mv2_, deg_mv2 + step_mv2, step_m_elems * sizeof (int), cudaMemcpyHostToDevice);
            cudaMemcpy (dir_mv2_, dir_mv2 + step_mv2, step_m_elems * sizeof (float), cudaMemcpyHostToDevice);
            cudaMemcpy (idx_mv2_, idx_mv2 + step_mv2, step_m_elems * sizeof (unsigned long long int), cudaMemcpyHostToDevice);

            // # total of elments to be computed.

            int step_n_m_elems = step_n_elems * step_m_elems;

            // # operate them as desired.

            kernel_metric_inner_prd <<<ceil (step_n_elems / (float) CALIB_N_THREADS_PER_BLOCK), CALIB_N_THREADS_PER_BLOCK>>> (step_n_elems, step_m_elems, deg_mv1_, dir_mv1_, idx_mv1_, deg_mv2_, dir_mv2_, idx_mv2_, inner_result_);

            kernel_metric_sum_reduct <<<ceil (step_n_elems / (float) CALIB_N_THREADS_PER_BLOCK), CALIB_N_THREADS_PER_BLOCK>>> (inner_result_);

            cudaDeviceSynchronize ();

            // # recover result packet from device.

            cudaMemcpy (inner_result, inner_result_, step_n_elems * sizeof (int), cudaMemcpyDeviceToHost);

            acc_result += inner_result [0];

            shift += step_n_m_elems;
        }
    }

    // # free memory.
    cudaFree (&deg_mv1_); cudaFree (&deg_mv2_);
    cudaFree (&dir_mv1_); cudaFree (&dir_mv2_);
    cudaFree (&idx_mv1_); cudaFree (&idx_mv2_);
    
    cudaFree (&inner_result_);
    
    std::free (deg_mv1); std::free (deg_mv2);
    std::free (dir_mv1); std::free (dir_mv2);
    std::free (idx_mv1); std::free (idx_mv2);
    
    std::free (inner_result);

    return acc_result;
}