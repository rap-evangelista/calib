// # expected kernels.

__global__ void kernel_cannonical_reordering_multivectors (int n_, int * deg_, unsigned long long int * idx_, float * mgn_, int * ori_)
{
    __shared__ bool sorted [1];

    sorted [0] = false;

    int pos = threadIdx. x + (blockDim. x * blockIdx. x);

    while (!sorted [0])
    {
        sorted [0] = true;

        if (pos % 2 == 0 && pos + 1 < n_)
        {
            if (deg_ [pos] == deg_ [pos+1] && idx_ [pos] == idx_ [pos+1])
            {
                // # aglomera.
                auto scalar = (ori_ [pos] * mgn_ [pos]) + (ori_ [pos+1] * mgn_ [pos+1]); //mv. elems [i]. direction () + mv. elems [i+1]. direction ();

                mgn_ [pos] = fabsf (scalar);
                ori_ [pos] *= signbit (scalar) ? -1 : 1;
            } else {
                if ((deg_ [pos] > deg_ [pos+1]) || (deg_ [pos] == deg_ [pos+1]) && (idx_ [pos] > idx_ [pos+1]))
                {
                    int aux_deg = deg_ [pos];
                    int aux_idx = idx_ [pos];
                    int aux_mgn = mgn_ [pos];
                    int aux_ori = ori_ [pos];

                    deg_ [pos] = deg_ [pos+1];
                    idx_ [pos] = idx_ [pos+1];
                    mgn_ [pos] = mgn_ [pos+1];
                    ori_ [pos] = ori_ [pos+1];

                    deg_ [pos+1] = aux_deg;
                    idx_ [pos+1] = aux_idx;
                    mgn_ [pos+1] = aux_mgn;
                    ori_ [pos+1] = aux_ori;
                }
            }
        }

        __syncthreads ();
        
        if (pos % 2 != 0 && pos + 1 < n_)
        {
            if (deg_ [pos] == deg_ [pos+1] && idx_ [pos] == idx_ [pos+1])
            {
                // # aglomera.
                auto scalar = (ori_ [pos] * mgn_ [pos]) + (ori_ [pos+1] * mgn_ [pos+1]); //mv. elems [i]. direction () + mv. elems [i+1]. direction ();

                mgn_ [pos] = fabsf (scalar);
                ori_ [pos] *= signbit (scalar) ? -1 : 1;
            } else {
                if ((deg_ [pos] > deg_ [pos+1]) || (deg_ [pos] == deg_ [pos+1]) && (idx_ [pos] > idx_ [pos+1]))
                {
                    int aux_deg = deg_ [pos];
                    int aux_idx = idx_ [pos];
                    int aux_mgn = mgn_ [pos];
                    int aux_ori = ori_ [pos];

                    deg_ [pos] = deg_ [pos+1];
                    idx_ [pos] = idx_ [pos+1];
                    mgn_ [pos] = mgn_ [pos+1];
                    ori_ [pos] = ori_ [pos+1];

                    deg_ [pos+1] = aux_deg;
                    idx_ [pos+1] = aux_idx;
                    mgn_ [pos+1] = aux_mgn;
                    ori_ [pos+1] = aux_ori;
                }
            }
        }

        __syncthreads ();
    }
}

__device__ void device_cannonical_reordering (int deg, int & deg_out, int * base_index, int & orientation)
{
    deg_out = deg;

    bool sorted = false;

    while (!sorted)
    {
        sorted = true;
        for (int i = 0; i < deg - 1; i++)
        {
            // # se bases iguais, elas se eliminam
            if (base_index [i] == base_index [i+1])
            {
                base_index [i] = -1;
                base_index [i+1] = -1;
                deg_out -= 2;

                if (i < 0) break;
                continue;
            }

            // # se base maior, troca as bases e inverte a orientação.
            if (base_index [i] > base_index [i+1])
            {
                sorted = false;
                int aux = base_index [i];
                base_index [i] = base_index [i+1];
                base_index [i+1] = aux;
                orientation *= -1;
            }
        }
    }
}

__device__ unsigned long long int pair_to_uid (unsigned long long int x, unsigned long long int y)
{
    return (x * x + 3 * x + 2 * x * y + y + y * y) / 2;
}

__device__ void uid_to_pair (unsigned long long int uid, unsigned long long int & x, unsigned long long int & y)
{
    unsigned long long int aux = floorf ((-1 + sqrtf (1 + 8 * uid)) / 2);

    x = uid - (aux * (1 + aux) / 2);
    y = (aux * (3 + aux) / 2) - uid;
}

__device__ void device_outer_product (const int deg1, unsigned long long int idx1, const int deg2, unsigned long long int idx2, int & deg3, unsigned long long int & idx3, int & orientation)
{
    // # unroll both indices.

    int smem_elems [2 * DEFAULT_SPACE_DIM];

    unsigned long long int elem_aux = idx1;

    if (deg1 == 1) smem_elems [deg1 - 1] = elem_aux;

    if (deg1 > 1)
    {
        for (int i = 0; i < deg1; i++)
        {
            unsigned long long int x = -1LLU;
            unsigned long long int y = -1LLU;

            uid_to_pair (elem_aux, x, y);

            elem_aux = x;
            smem_elems [deg1 - 1 - i] = y;
        }
    }

    elem_aux = idx2;

    if (deg2 == 1) smem_elems [deg2 - 1] = elem_aux;

    if (deg2 > 1)
    {
        for (int i = 0; i < deg2; i++)
        {
            unsigned long long int x = -1LLU;
            unsigned long long int y = -1LLU;

            uid_to_pair (elem_aux, x, y);

            elem_aux = x;
            smem_elems [deg1 + deg2 - 1 - i] = y;
        }
    }

    // # cannonical reordering on the shared memory.

    __syncthreads ();

    orientation = 1;

    device_cannonical_reordering (deg1 + deg2, deg3, smem_elems, orientation);

    __syncthreads ();

    // # re-calculate the 'idx'.

    idx3 = (deg3 <= 0) ? 0 : smem_elems [deg1 + deg2 - deg3];

    if (deg3 <= 1) return;

    for (int i = 1; i < deg3; i++)
    {
        idx3 = pair_to_uid (idx3, smem_elems [i + (deg1 + deg2 - deg3)]);
    }
}

__device__ void device_regr_product (const int deg1, unsigned long long int idx1, const int deg2, unsigned long long int idx2, int & deg3, unsigned long long int & idx3, int & orientation)
{
    // # unroll both indices.

    unsigned long long int idx1_aux = idx1;
    unsigned long long int idx2_aux = idx2;
    
    if (deg1 <= 0 || deg2 <= 0)
    {
        orientation = 1;
        deg3 = 0;
        idx3 = 0;
        return;
    }

    if (deg1 == 1 && deg2 == 1)
    {
        if (idx1_aux == idx2_aux)
        {
            orientation = 1;
            deg3 = 1;
            idx3 = idx1_aux;
        } else {
            orientation = 1;
            deg3 = 0;
            idx3 = 0;
        }

        return;
    }

    int smem_elems1 [DEFAULT_SPACE_DIM];
    int smem_elems2 [DEFAULT_SPACE_DIM];
    int smem_elems3 [DEFAULT_SPACE_DIM];

    for (int i = 0; i < deg1; i++)
    {
        if (i == deg1 - 1)
        {
            smem_elems1 [deg1 - 1 - i] = idx1_aux;
            break;
        }

        unsigned long long int x = -1LLU;
        unsigned long long int y = -1LLU;

        uid_to_pair (idx1_aux, x, y);

        idx1_aux = x;
        smem_elems1 [deg1 - 1 - i] = y;
    }

    for (int i = 0; i < deg2; i++)
    {
        if (i == deg2 - 1)
        {
            smem_elems2 [deg2 - 1 - i] = idx2_aux;
            break;
        }

        unsigned long long int x = -1LLU;
        unsigned long long int y = -1LLU;

        uid_to_pair (idx2_aux, x, y);

        idx2_aux = x;
        smem_elems2 [deg2 - 1 - i] = y;
    }

    int i = 0, j = 0, k = 0;
    while (i < deg1 && j < deg2)
    {
        if (smem_elems1 [i] == smem_elems2 [j])
        {
            smem_elems3 [k++] = smem_elems1 [i];
            i++; j++;
        } else {
            if (smem_elems1 [i] < smem_elems2 [j])
            {
                i++;
            } else {
                j++;
            }
        }
    }

    deg3 = k;

    __syncthreads ();

    // # re-calculate the 'idx'.

    idx3 = (deg3 <= 0) ? 0 : smem_elems3 [0];

    if (deg3 <= 1) return;

    for (int i = 1; i < deg3; i++)
    {
        idx3 = pair_to_uid (idx3, smem_elems3 [i]);
    }
}

__global__ void kernel_outer_prd_multivectors (int n_, int m_, int * deg_mv1_, float * dir_mv1_, unsigned long long int * idx_mv1_, int * deg_mv2_, float * dir_mv2_, unsigned long long int * idx_mv2_, int * deg_mv3_, float * dir_mv3_, unsigned long long int * idx_mv3_)
{
    int pos = threadIdx. x + (blockDim. x * blockIdx. x);

    // # partition size sent to device for each time.
    
    const int partition_ = 512;

    int times_ = std::ceil (n_ / (float) partition_);

    // # shared memory allocation.
    int smem_deg [partition_];
    float smem_dir [partition_];
    unsigned long long int smem_idx [partition_];

    // # iterate over the mv#1, filling 'smem' each time.

    if (pos < m_)
    {
        for (int i = 0; i < times_; i++)
        {
            int smem_step = (i * partition_);
            int smem_n_elems = partition_;
        
            if (smem_step + smem_n_elems > n_)
            {
                smem_n_elems = n_ - smem_step;
            }

            // # fill shared memory with current partition of mv#1.

            for (int x = 0; x < smem_n_elems; x++)
            {
                if (smem_step + x < n_)
                {
                    smem_deg [x] = deg_mv1_ [smem_step + x];
                    smem_dir [x] = dir_mv1_ [smem_step + x];
                    smem_idx [x] = idx_mv1_ [smem_step + x];
                }
            }

            // # cross operate mv#2 with the shared memory.
            for (int j = 0; j < n_; j++)
            {
                int stride = m_ * j;

                int deg_aux = -1;
                unsigned long long int idx_aux = -1LLU;

                int orientation;

                device_outer_product (deg_mv2_ [pos], idx_mv2_ [pos], smem_deg [j % smem_n_elems], smem_idx [j % smem_n_elems], deg_aux, idx_aux, orientation);

                deg_mv3_ [pos + stride] = deg_aux;
                idx_mv3_ [pos + stride] = idx_aux;
                dir_mv3_ [pos + stride] = orientation * smem_dir [j % smem_n_elems] * dir_mv2_ [pos];
            }
        }

        __syncthreads ();
    }
}

__global__ void kernel_regr_prd_multivectors (int n_, int m_, int * deg_mv1_, float * dir_mv1_, unsigned long long int * idx_mv1_, int * deg_mv2_, float * dir_mv2_, unsigned long long int * idx_mv2_, int * deg_mv3_, float * dir_mv3_, unsigned long long int * idx_mv3_)
{
    int pos = threadIdx. x + (blockDim. x * blockIdx. x);

    // # partition size sent to device for each time.
    
    const int partition_ = 512;

    int times_ = std::ceil (n_ / (float) partition_);

    // # shared memory allocation.
    int smem_deg [partition_];
    float smem_dir [partition_];
    unsigned long long int smem_idx [partition_];

    // # iterate over the mv#1, filling 'smem' each time.

    if (pos < m_)
    {
        for (int i = 0; i < times_; i++)
        {
            int smem_step = (i * partition_);
            int smem_n_elems = partition_;
        
            if (smem_step + smem_n_elems > n_)
            {
                smem_n_elems = n_ - smem_step;
            }

            // # fill shared memory with current partition of mv#1.

            for (int x = 0; x < smem_n_elems; x++)
            {
                if (smem_step + x < n_)
                {
                    smem_deg [x] = deg_mv1_ [smem_step + x];
                    smem_dir [x] = dir_mv1_ [smem_step + x];
                    smem_idx [x] = idx_mv1_ [smem_step + x];
                }
            }

            // # cross operate mv#2 with the shared memory.
            for (int j = 0; j < n_; j++)
            {
                int stride = m_ * j;

                int deg_aux = -1;
                unsigned long long int idx_aux = -1LLU;

                int orientation;

                device_regr_product (deg_mv2_ [pos], idx_mv2_ [pos], smem_deg [j % smem_n_elems], smem_idx [j % smem_n_elems], deg_aux, idx_aux, orientation);

                deg_mv3_ [pos + stride] = deg_aux;
                idx_mv3_ [pos + stride] = idx_aux;
                dir_mv3_ [pos + stride] = orientation * smem_dir [j % smem_n_elems] * dir_mv2_ [pos];
            }
        }

        __syncthreads ();
    }
}

__global__ void kernel_inner_prd_multivectors ()
{
    printf ("'kernel_inner_prd_multivectors' not implemented yet.\n");
}

__global__ void kernel_left_contr_multivectors ()
{
    printf ("'kernel_left_contr_multivectors' not implemented yet.\n");
}

__global__ void kernel_geom_prd_multivectors ()
{
    printf ("'kernel_geom_prd_multivectors' not implemented yet.\n");
}


namespace calib
{
    // # bridge operations between host code and device code.

    void bridge_cannonical_reordering_multivectors (multivector& mv)
    {
        int n_ = mv. elems. size ();

        // # host-side variables.
        int * deg;
        unsigned long long int * idx;
        float * mgn;
        int * ori;

        deg = (int *) std::malloc (n_ * sizeof (int));
        idx = (unsigned long long int *) std::malloc (n_ * sizeof (unsigned long long int));
        mgn = (float *) std::malloc (n_ * sizeof (float));
        ori = (int *) std::malloc (n_ * sizeof (int));

        // # device-side variables.
        int * deg_;
        unsigned long long int * idx_;
        float * mgn_;
        int * ori_;

        cudaMalloc ((void**) &deg_, n_ * sizeof (int));
        cudaMalloc ((void**) &idx_, n_ * sizeof (unsigned long long int));
        cudaMalloc ((void**) &mgn_, n_ * sizeof (float));
        cudaMalloc ((void**) &ori_, n_ * sizeof (int));

        // # prepare host data and copy them to device.
        
        for (int i = 0; i < n_; i++)
        {
            basis base = mv. elems [i];

            deg [i] = base. degree ();
            idx [i] = base. unique_index ();
            mgn [i] = base. magnitude;
            ori [i] = base. orientation;
        }

        cudaMemcpy (deg_, deg, n_ * sizeof (int), cudaMemcpyHostToDevice);
        cudaMemcpy (idx_, idx, n_ * sizeof (unsigned long long int), cudaMemcpyHostToDevice);
        cudaMemcpy (mgn_, mgn, n_ * sizeof (float), cudaMemcpyHostToDevice);
        cudaMemcpy (ori_, ori, n_ * sizeof (int), cudaMemcpyHostToDevice);

        kernel_cannonical_reordering_multivectors <<<ceil (n_ / (float) CALIB_N_THREADS_PER_BLOCK), CALIB_N_THREADS_PER_BLOCK>>> (n_, deg_, idx_, mgn_, ori_);

        cudaDeviceSynchronize ();

        cudaMemcpy (deg, deg_, n_ * sizeof (int), cudaMemcpyDeviceToHost);
        cudaMemcpy (idx, idx_, n_ * sizeof (unsigned long long int), cudaMemcpyDeviceToHost);
        cudaMemcpy (mgn, mgn_, n_ * sizeof (float), cudaMemcpyDeviceToHost);
        cudaMemcpy (ori, ori_, n_ * sizeof (int), cudaMemcpyDeviceToHost);

        // # copy result data back to 'mv'.

        mv. clear ();

        for (int i = 0; i < n_; i++)
        {
            basis base = basis ();

            base. recover_from_degree_and_unique_index (deg [i], idx [i]);

            base. magnitude = mgn [i];
            base. orientation = ori [i];

            mv. add_elem (base);
        }

        cudaFree (&deg_);
        cudaFree (&idx_);
        cudaFree (&mgn_);
        cudaFree (&ori_);

        std::free (deg);
        std::free (idx);
        std::free (mgn);
        std::free (ori);
    }

    multivector bridge_outer_prd_multivectors (multivector& mv1, multivector& mv2)
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
            basis base = mv1. elems [i];

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
            basis base = mv2. elems [i];

            deg_mv2 [i] = base. degree ();
            dir_mv2 [i] = base. direction ();
            idx_mv2 [i] = base. unique_index ();
        }

        // # host-side memory allocation [mv #3].
        int * deg_mv3;
        float * dir_mv3;
        unsigned long long int * idx_mv3;

        deg_mv3 = (int *) std::malloc (n_ * m_ * sizeof (int));
        dir_mv3 = (float *) std::malloc (n_ * m_ * sizeof (float));
        idx_mv3 = (unsigned long long int *) std::malloc (n_ * m_ * sizeof (unsigned long long int));

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
        int * deg_mv3_;
        float * dir_mv3_;
        unsigned long long int * idx_mv3_;

        cudaMalloc ((void **) &deg_mv3_, partition_ * sizeof (int));
        cudaMalloc ((void **) &dir_mv3_, partition_ * sizeof (float));
        cudaMalloc ((void **) &idx_mv3_, partition_ * sizeof (unsigned long long int));

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

                kernel_outer_prd_multivectors <<<ceil (step_n_m_elems / (float) CALIB_N_THREADS_PER_BLOCK), CALIB_N_THREADS_PER_BLOCK>>> (step_n_elems, step_m_elems, deg_mv1_, dir_mv1_, idx_mv1_, deg_mv2_, dir_mv2_, idx_mv2_, deg_mv3_, dir_mv3_, idx_mv3_);

                cudaDeviceSynchronize ();

                // # recover result packet from device.

                cudaMemcpy (deg_mv3 + shift, deg_mv3_, step_n_m_elems * sizeof (int), cudaMemcpyDeviceToHost);
                cudaMemcpy (dir_mv3 + shift, dir_mv3_, step_n_m_elems * sizeof (float), cudaMemcpyDeviceToHost);
                cudaMemcpy (idx_mv3 + shift, idx_mv3_, step_n_m_elems * sizeof (unsigned long long int), cudaMemcpyDeviceToHost);

                shift += step_n_m_elems;
            }
        }

        // # copy result array.

        multivector mv = multivector ();

        for (int i = 0; i < n_ * m_; i++)
        {
            basis base = basis ();
            
            base. recover_from_degree_and_unique_index (deg_mv3 [i], idx_mv3 [i]);

            base. magnitude = std::abs (dir_mv3 [i]);
            base. orientation = std::copysign (1, dir_mv3 [i]);
            
            mv. add_elem (base);
        }

        // # free memory.

        cudaFree (&deg_mv1_); cudaFree (&deg_mv2_); cudaFree (&deg_mv3_);
        cudaFree (&dir_mv1_); cudaFree (&dir_mv2_); cudaFree (&dir_mv3_);
        cudaFree (&idx_mv1_); cudaFree (&idx_mv2_); cudaFree (&idx_mv3_);
        
        std::free (deg_mv1); std::free (deg_mv2); std::free (deg_mv3);
        std::free (dir_mv1); std::free (dir_mv2); std::free (dir_mv3);
        std::free (idx_mv1); std::free (idx_mv2); std::free (idx_mv3);

        return mv;
    }

    multivector bridge_regr_prd_multivectors (multivector& mv1, multivector& mv2)
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
            basis base = mv1. elems [i];

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
            basis base = mv2. elems [i];

            deg_mv2 [i] = base. degree ();
            dir_mv2 [i] = base. direction ();
            idx_mv2 [i] = base. unique_index ();
        }

        // # host-side memory allocation [mv #3].
        int * deg_mv3;
        float * dir_mv3;
        unsigned long long int * idx_mv3;

        deg_mv3 = (int *) std::malloc (n_ * m_ * sizeof (int));
        dir_mv3 = (float *) std::malloc (n_ * m_ * sizeof (float));
        idx_mv3 = (unsigned long long int *) std::malloc (n_ * m_ * sizeof (unsigned long long int));

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
        int * deg_mv3_;
        float * dir_mv3_;
        unsigned long long int * idx_mv3_;

        cudaMalloc ((void **) &deg_mv3_, partition_ * sizeof (int));
        cudaMalloc ((void **) &dir_mv3_, partition_ * sizeof (float));
        cudaMalloc ((void **) &idx_mv3_, partition_ * sizeof (unsigned long long int));

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

                kernel_regr_prd_multivectors <<<ceil (step_n_m_elems / (float) CALIB_N_THREADS_PER_BLOCK), CALIB_N_THREADS_PER_BLOCK>>> (step_n_elems, step_m_elems, deg_mv1_, dir_mv1_, idx_mv1_, deg_mv2_, dir_mv2_, idx_mv2_, deg_mv3_, dir_mv3_, idx_mv3_);

                cudaDeviceSynchronize ();

                // # recover result packet from device.

                cudaMemcpy (deg_mv3 + shift, deg_mv3_, step_n_m_elems * sizeof (int), cudaMemcpyDeviceToHost);
                cudaMemcpy (dir_mv3 + shift, dir_mv3_, step_n_m_elems * sizeof (float), cudaMemcpyDeviceToHost);
                cudaMemcpy (idx_mv3 + shift, idx_mv3_, step_n_m_elems * sizeof (unsigned long long int), cudaMemcpyDeviceToHost);

                shift += step_n_m_elems;
            }
        }

        // # copy result array.

        multivector mv = multivector ();

        for (int i = 0; i < n_ * m_; i++)
        {
            basis base = basis ();
            
            base. recover_from_degree_and_unique_index (deg_mv3 [i], idx_mv3 [i]);

            base. magnitude = std::abs (dir_mv3 [i]);
            base. orientation = std::copysign (1, dir_mv3 [i]);
            
            mv. add_elem (base);
        }

        // # free memory.

        cudaFree (&deg_mv1_); cudaFree (&deg_mv2_); cudaFree (&deg_mv3_);
        cudaFree (&dir_mv1_); cudaFree (&dir_mv2_); cudaFree (&dir_mv3_);
        cudaFree (&idx_mv1_); cudaFree (&idx_mv2_); cudaFree (&idx_mv3_);
        
        std::free (deg_mv1); std::free (deg_mv2); std::free (deg_mv3);
        std::free (dir_mv1); std::free (dir_mv2); std::free (dir_mv3);
        std::free (idx_mv1); std::free (idx_mv2); std::free (idx_mv3);

        return mv;
    }

    multivector bridge_left_contr_multivectors (multivector& mv1, multivector& mv2)
    {
        kernel_left_contr_multivectors <<<1, 1>>> ();
        return multivector ();
    }
}