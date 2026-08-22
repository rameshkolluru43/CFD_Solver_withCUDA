void solvePressureCorrection(
    double** p, double** u_star, double** v_star, int nx, int ny, double dx, double dy, double dt, double rho
) {
    double** p_prime = new double*[nx];
    for (int i = 0; i < nx; ++i)
        p_prime[i] = new double[ny]();

    for (int iter = 0; iter < 100; ++iter) { // Example fixed iteration count
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                // Divergence of intermediate velocity
                double div_velocity = ((u_star[i+1][j] - u_star[i][j]) / dx) + ((v_star[i][j+1] - v_star[i][j]) / dy);

                // Pressure correction (example stencil)
                p_prime[i][j] = (div_velocity * rho * dx * dy) / (2 * (dx + dy));
            }
        }

        // Update pressure field
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                p[i][j] += p_prime[i][j];
            }
        }
    }

    // Free memory
    for (int i = 0; i < nx; ++i)
        delete[] p_prime[i];
    delete[] p_prime;
}