void solveMomentumSemiImplicit(
    double** u, double** v, double** u_old, double** v_old,
    double** p, int nx, int ny, double dx, double dy, double dt, double rho, double mu
) {
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            // Explicit convective terms for u
            double convection_u = -(u_old[i][j] * (u_old[i+1][j] - u_old[i-1][j]) / (2 * dx))
                                - (v_old[i][j] * (u_old[i][j+1] - u_old[i][j-1]) / (2 * dy));

            // Implicit diffusion and pressure terms for u
            double diffusion_u = mu * (u_old[i+1][j] - 2 * u_old[i][j] + u_old[i-1][j]) / (dx * dx)
                               + mu * (u_old[i][j+1] - 2 * u_old[i][j] + u_old[i][j-1]) / (dy * dy);

            double pressure_u = (p[i][j] - p[i-1][j]) / dx;

            // Semi-implicit update for u-velocity
            u[i][j] = u_old[i][j] + dt * (convection_u + diffusion_u - pressure_u);
        }
    }

    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            // Explicit convective terms for v
            double convection_v = -(u_old[i][j] * (v_old[i+1][j] - v_old[i-1][j]) / (2 * dx))
                                - (v_old[i][j] * (v_old[i][j+1] - v_old[i][j-1]) / (2 * dy));

            // Implicit diffusion and pressure terms for v
            double diffusion_v = mu * (v_old[i+1][j] - 2 * v_old[i][j] + v_old[i-1][j]) / (dx * dx)
                               + mu * (v_old[i][j+1] - 2 * v_old[i][j] + v_old[i][j-1]) / (dy * dy);

            double pressure_v = (p[i][j] - p[i][j-1]) / dy;

            // Semi-implicit update for v-velocity
            v[i][j] = v_old[i][j] + dt * (convection_v + diffusion_v - pressure_v);
        }
    }
}