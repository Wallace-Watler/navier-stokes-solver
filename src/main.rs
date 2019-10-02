fn main() {
    const PHYSICAL_X: f64 = 1.0;
    const PHYSICAL_Y: f64 = 1.0;
    const GRID_SIZE: usize = 7;
    const DX: f64 = PHYSICAL_X / (GRID_SIZE - 1) as f64;
    const DY: f64 = PHYSICAL_Y / (GRID_SIZE - 1) as f64;
    const DT: f64 = 0.001;
    const MU: f64 = 0.01;
    const DELTA: f64 = 4.5;
    const ACCEPTABLE_ERROR: f64 = 1e-8;
    const MAX_ITERATIONS: u32 = 1000000;

    let mut u_vel = [[0.0; GRID_SIZE + 1]; GRID_SIZE];
    let mut u_vel_new = [[0.0; GRID_SIZE + 1]; GRID_SIZE];
    let mut v_vel = [[0.0; GRID_SIZE]; GRID_SIZE + 1];
    let mut v_vel_new = [[0.0; GRID_SIZE]; GRID_SIZE + 1];
    let mut pressure = [[1.0; GRID_SIZE + 1]; GRID_SIZE + 1];

    // Initialize velocities at boundaries
    for x in 0..GRID_SIZE {
        u_vel[x][GRID_SIZE] = 2.0;
    }

    let mut error = 1.0;
    let mut iterations: u32 = 1;

    while iterations < MAX_ITERATIONS && error > ACCEPTABLE_ERROR {
        // Solve u-velocity equation
        for y in 1..GRID_SIZE {
            for x in 1..(GRID_SIZE - 1) {
                let up = u_vel[x][y]; let ue = u_vel[x + 1][y]; let uw = u_vel[x - 1][y]; let un = u_vel[x][y + 1]; let us = u_vel[x][y - 1];
                let v1 = v_vel[x][y - 1]; let v2 = v_vel[x + 1][y - 1]; let v3 = v_vel[x + 1][y]; let v4 = v_vel[x][y];
                let pe = pressure[x + 1][y]; let pw = pressure[x][y];
                let viscous_term = MU * ((ue - 2.0 * up + uw) / (DX * DX) + (un - 2.0 * up + us) / (DY * DY));
                let pressure_diff = (pe - pw) / DX;
                let u_divergence = (ue * ue - uw * uw) / (2.0 * DX);
                let v_divergence = ((up + un) * (v3 + v4) - (up + us) * (v1 + v2)) / (4.0 * DY);
                u_vel_new[x][y] = up + (viscous_term - pressure_diff - u_divergence - v_divergence) * DT;
            }
        }

        // Set u boundary conditions
        for y in 1..GRID_SIZE {
            // Velocities lie on the boundary, so simply set the condition
            u_vel_new[0][y] = 0.0;
            u_vel_new[GRID_SIZE - 1][y] = 0.0;
        }

        for x in 0..GRID_SIZE {
            // Boundary lies between velocities; ensure that the outside velocity is such that
            // its average with the inside velocity equals the desired boundary condition
            u_vel_new[x][0] = 2.0 * 0.0 - u_vel_new[x][1];
            u_vel_new[x][GRID_SIZE] = 2.0 * 1.0 - u_vel_new[x][GRID_SIZE - 1];
        }

        // Solve v-velocity equation
        for y in 1..(GRID_SIZE - 1) {
            for x in 1..GRID_SIZE {
                let vp = v_vel[x][y]; let ve = v_vel[x + 1][y]; let vw = v_vel[x - 1][y]; let vn = v_vel[x][y + 1]; let vs = v_vel[x][y - 1];
                let u1 = u_vel[x - 1][y]; let u2 = u_vel[x][y]; let u3 = u_vel[x][y + 1]; let u4 = u_vel[x - 1][y + 1];
                let pn = pressure[x][y + 1]; let ps = pressure[x][y];
                let viscous_term = MU * ((ve - 2.0 * vp + vw) / (DX * DX) + (vn - 2.0 * vp + vs) / (DY * DY));
                let pressure_diff = (pn - ps) / DY;
                let v_divergence = (vn * vn - vs * vs) / (2.0 * DY);
                let u_divergence = ((vp + vw) * (u1 + u4) - (vp + ve) * (u2 + u3)) / (4.0 * DX);
                v_vel_new[x][y] = vp + (viscous_term - pressure_diff - v_divergence - u_divergence) * DT;
            }
        }

        // Set v boundary conditions
        for x in 1..GRID_SIZE {
            // Velocities lie on the boundary, so simply set the condition
            v_vel_new[x][0] = 0.0;
            v_vel_new[x][GRID_SIZE - 1] = 0.0;
        }

        for y in 0..GRID_SIZE {
            // Boundary lies between velocities; ensure that the outside velocity is such that
            // its average with the inside velocity equals the desired boundary condition
            v_vel_new[0][y] = 2.0 * 0.0 - v_vel_new[1][y];
            v_vel_new[GRID_SIZE][y] = 2.0 * 0.0 - v_vel_new[GRID_SIZE - 1][y];
        }

        error = 0.0;

        // Solve continuity equation using new velocities
        for y in 1..GRID_SIZE {
            for x in 1..GRID_SIZE {
                let ue = u_vel_new[x][y]; let uw = u_vel_new[x - 1][y];
                let vn = v_vel_new[x][y]; let vs = v_vel_new[x][y - 1];
                let divergence = (ue - uw) / DX + (vn - vs) / DY;
                pressure[x][y] -= DELTA * divergence * DT;

                // Error is the sum of velocity divergences across grid
                error += divergence.abs();
            }
        }

        // Set pressure boundary conditions
        // No need, all boundaries are set to velocities

        // Update velocities
        u_vel = u_vel_new;
        v_vel = v_vel_new;

        println!("Error: {}", error);

        iterations += 1;
    }
}
