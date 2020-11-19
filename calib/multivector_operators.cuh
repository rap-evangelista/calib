// # expected kernels.

__global__ void kernel_sum_multivectors ()
{
    printf ("'kernel_sum_multivectors' not implemented yet.\n");
}

__global__ void kernel_outer_prd_multivectors ()
{
    printf ("'kernel_outer_prd_multivectors' not implemented yet.\n");
}

__global__ void kernel_regr_prd_multivectors ()
{
    printf ("'kernel_regr_prd_multivectors' not implemented yet.\n");
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

// # bridge operations between host code and device code.

const calib::multivector bridge_sum_multivectors (calib::multivector& mv1, calib::multivector& mv2)
{
    kernel_sum_multivectors <<<1, 1>>> ();
    return calib::multivector ();
}

const calib::multivector bridge_outer_prd_multivectors (calib::multivector& v1, calib::multivector& v2)
{
    kernel_outer_prd_multivectors <<<1, 1>>> ();
    return calib::multivector ();
}

const calib::multivector bridge_regr_prd_multivectors (calib::multivector& v1, calib::multivector& v2)
{
    kernel_regr_prd_multivectors <<<1, 1>>> ();
    return calib::multivector ();
}

const calib::multivector bridge_inner_prd_multivectors (calib::multivector& v1, calib::multivector& v2)
{
    kernel_inner_prd_multivectors <<<1, 1>>> ();
    return calib::multivector ();
}

const calib::multivector bridge_left_contr_multivectors (calib::multivector& v1, calib::multivector& v2)
{
    kernel_left_contr_multivectors <<<1, 1>>> ();
    return calib::multivector ();
}

const calib::multivector bridge_geom_prd_multivectors (calib::multivector& v1, calib::multivector& v2)
{
    kernel_geom_prd_multivectors <<<1, 1>>> ();
    return calib::multivector ();
}