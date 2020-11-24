// # expected kernels.

__global__ void kernel_metric_inner_prd ()
{
    printf ("'kernel_cannonical_reordering' not implemented yet.\n");
}

// # bridge operations between host code and device code.

double bridge_metric_inner_prd (calib::multivector& m1, calib::multivector& m2, std::vector <double> & metric)
{
    kernel_cannonical_reordering <<<1, 1>>> ();
    return 0;
}