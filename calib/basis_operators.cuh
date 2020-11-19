// # expected kernels.

__global__ void kernel_cannonical_reordering ()
{
    printf ("'kernel_cannonical_reordering' not implemented yet.\n");
}

// # bridge operations between host code and device code.

void bridge_cannonical_reordering (calib::basis& base)
{
    kernel_cannonical_reordering <<<1, 1>>> ();
}