// =========================================================================
// CubeSat Attitude Control Simulation with Magnetorquers
// =========================================================================
// This simulation models a 3U CubeSat with magnetorquer-based attitude
// control using PID controllers. The simulation includes:
//  - Full quaternion-based attitude dynamics
//  - Realistic orbital magnetic field model
//  - PID control implementation with anti-windup
//  - Proper magnetorquer constraints (perpendicular torque only)
//  - Professional visualization of results
//  - ODE solver enhancements for stiff systems
// =========================================================================

// Clear workspace and close all figures
clear;
clc;
close(winsid());

// =========================================================================
// SATELLITE PARAMETERS
// =========================================================================
// 3U CubeSat parameters (realistic values)
J = diag([0.0035, 0.0342, 0.0342]); // Inertia matrix (kg·m²)
max_dipole = 0.2;                    // Maximum magnetic dipole moment (A·m²)

// =========================================================================
// CONTROLLER PARAMETERS
// =========================================================================
// PID controller gains (reduced for numerical stability)
Kp = [0.1; 0.1; 0.1];      // Proportional gain (reduced from 0.5)
Ki = [0.005; 0.005; 0.005]; // Integral gain (reduced from 0.01)
Kd = [0.1; 0.1; 0.1];      // Derivative gain (reduced from 0.3)

// =========================================================================
// SIMULATION PARAMETERS
// =========================================================================
// Orbit parameters
altitude = 500e3;          // Orbit altitude (m)
inclination = 51.6*%pi/180; // Orbit inclination (ISS-like, rad)
orbital_period = 2*%pi*sqrt((6371e3 + altitude)^3/(3.986e14)); // Seconds

// Initial conditions (reduced tumbling for stability)
q0 = [0.1; 0.05; 0.05; 0.99]; // Initial quaternion [x,y,z,w] (less deviation)
q0 = q0/norm(q0);           // Normalize quaternion
omega0 = [0.01; -0.01; 0.005]; // Initial angular velocity (rad/s) (reduced)

// Reference attitude (desired pointing)
q_ref = [0; 0; 0; 1];       // Identity quaternion (aligned with orbit frame)

// Simulation time
sim_duration = 1.5*orbital_period; // Simulate for 1.5 orbits
dt = 0.5;                   // Time step (s) (reduced from 1.0)
t = 0:dt:sim_duration;      // Time vector

// =========================================================================
// HELPER FUNCTIONS
// =========================================================================

// Magnetic field model (IGRF-based simplified model)
function B = magnetic_field(t)
    // Orbital position (simplified circular orbit)
    theta = 2*%pi*t/orbital_period;

    // Earth's magnetic field components in orbit frame (typical LEO values)
    // Using smoother field model to reduce stiffness
    B_orbit = [2e-5*cos(theta);
               2e-5*sin(theta)*sin(inclination);
               -3e-5*sin(inclination)];

    // Reduced random variations for stability
    B_orbit = B_orbit .* (1 + 0.01*rand(3,1,"normal"));

    B = B_orbit;
endfunction

// Quaternion multiplication
function q = quat_mult(q1, q2)
    // Extract scalar and vector parts
    v1 = q1(1:3); s1 = q1(4);
    v2 = q2(1:3); s2 = q2(4);

    // Quaternion product
    s = s1*s2 - sum(v1.*v2);
    v = s1*v2 + s2*v1 + cross(v1, v2);

    q = [v; s];
endfunction

// Quaternion conjugate
function q_conj = quat_conj(q)
    q_conj = [-q(1:3); q(4)];
endfunction

// Rotate vector from orbit to body frame using quaternion
function v_body = quat_rotate(q, v_orbit)
    // Create pure quaternion from vector
    qv = [v_orbit; 0];

    // Perform rotation: q * v * q^(-1)
    q_inv = quat_conj(q);
    temp = quat_mult(q, qv);
    result = quat_mult(temp, q_inv);

    v_body = result(1:3);
endfunction

// Skew-symmetric matrix for cross product
function S = skew(v)
    S = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0];
endfunction

// =========================================================================
// ATTITUDE DYNAMICS
// =========================================================================
function dX = attitude_dynamics(t, X)
    // State vector X = [q(4); omega(3); integral_error(3)]
    q = X(1:4);
    q = q/norm(q);  // Normalize quaternion
    omega = X(5:7);
    integral_error = X(8:10);

    // Get magnetic field in orbit frame
    B_orbit = magnetic_field(t);

    // Rotate magnetic field to body frame
    B_body = quat_rotate(q, B_orbit);

    // Calculate attitude error (q_err = q_ref * q^(-1))
    q_inv = quat_conj(q);
    q_err = quat_mult(q_ref, q_inv);

    // Extract error vector (small angle approximation)
    error = q_err(1:3);

    // Ensure shortest rotation path (handle double-covering of SO(3))
    if q_err(4) < 0 then
        error = -error;
    end

    // Angular velocity error
    omega_err = omega; // Reference is zero

    // Update integral error with anti-windup
    integral_error = integral_error + error*dt;
    max_integral = 0.5; // Reduced from 1.0 for stability
    for i = 1:3
        if abs(integral_error(i)) > max_integral then
            integral_error(i) = sign(integral_error(i))*max_integral;
        end
    end

    // PID control law
    m_desired = Kp.*error + Ki.*integral_error + Kd.*omega_err;

    // Limit magnetic dipole moment (hardware constraint)
    for i = 1:3
        if abs(m_desired(i)) > max_dipole then
            m_desired(i) = sign(m_desired(i))*max_dipole;
        end
    end

    // Project control onto plane perpendicular to B
    // (magnetorquers can only produce torque perpendicular to B)
    B_unit = B_body/norm(B_body);
    m_parallel = (sum(m_desired.*B_unit))*B_unit;
    m_perpendicular = m_desired - m_parallel;

    // Calculate control torque (τ = m × B)
    tau_control = cross(m_perpendicular, B_body);

    // External disturbance torques (reduced for stability)
    tau_disturbance = [5e-8*sin(0.01*t);
                       5e-8*cos(0.02*t);
                       2.5e-8*sin(0.015*t)];

    // Total torque
    tau = tau_control + tau_disturbance;

    // Quaternion kinematics (q_dot = 0.5*q⊗[ω,0])
    omega_quat = [omega; 0];
    q_dot = 0.5*quat_mult(q, omega_quat);

    // Angular velocity dynamics (Euler's equation)
    omega_dot = inv(J)*(tau - cross(omega, J*omega));

    // Return state derivatives
    dX = [q_dot; omega_dot; error];
endfunction

// =========================================================================
// JACOBIAN FUNCTION FOR STIFF SOLVER
// =========================================================================
function J = attitude_jacobian(t, X)
    // Analytical Jacobian matrix for the system
    // State vector X = [q(4); omega(3); integral_error(3)]
    // Size of Jacobian is 10x10 (matching state vector size)

    // Initialize Jacobian matrix
    J = zeros(10, 10);

    // This is a simplified Jacobian that helps the solver
    // A full analytical Jacobian would be more complex

    // Quaternion kinematics part (∂q_dot/∂q and ∂q_dot/∂ω)
    // Diagonal elements for quaternion derivatives
    J(1:4, 1:4) = eye(4, 4) * 0.01;

    // Coupling between quaternion and angular velocity
    J(1:3, 5:7) = eye(3, 3) * 0.5;

    // Angular velocity dynamics (∂ω_dot/∂ω)
    // Diagonal elements for angular velocity
    J(5:7, 5:7) = eye(3, 3) * (-0.1);

    // Integral error dynamics
    J(8:10, 1:3) = eye(3, 3) * 0.1;

    // **FIX:** Do NOT use "return J;" here. Let the function implicitly return J.
endfunction

// =========================================================================
// MAIN SIMULATION
// =========================================================================
// Initial state vector
X0 = [q0; omega0; zeros(3,1)];

// Run simulation with stiff solver and Jacobian
disp("Starting CubeSat attitude control simulation...");
tic();

// Set ODE solver tolerances
rtol = 1e-4; // Relative tolerance (less strict)
atol = 1e-6; // Absolute tolerance (less strict)

// Use stiff solver with Jacobian and tolerances
// **Ensure this line correctly calls the Jacobian function**
X = ode("stiff", X0, 0, t, rtol, atol, attitude_dynamics, attitude_jacobian);

simulation_time = toc();
disp("Simulation completed in " + string(simulation_time) + " seconds");

// =========================================================================
// POST-PROCESSING
// =========================================================================
// Extract results
q_history = X(1:4,:);
omega_history = X(5:7,:);
integral_error_history = X(8:10,:); // Store integral error for analysis if needed

// Calculate pointing error over time
pointing_error = zeros(1, length(t));
for i = 1:length(t)
    q_current = q_history(:,i);
    q_current = q_current/norm(q_current);
    q_err = quat_mult(q_ref, quat_conj(q_current));

    // Convert to angle using q = [sin(θ/2)*axis, cos(θ/2)]
    pointing_error(i) = 2*acos(abs(q_err(4)))*180/%pi;  // in degrees
end

// Calculate control effort (magnetic dipole moment)
m_history = zeros(3, length(t));
tau_history = zeros(3, length(t));

for i = 1:length(t)
    // Get current state
    q = q_history(:,i);
    omega = omega_history(:,i);
    integral_error = integral_error_history(:,i); // Use stored integral error

    // Get magnetic field
    B_orbit = magnetic_field(t(i));
    B_body = quat_rotate(q, B_orbit);

    // Calculate error vector
    q_err = quat_mult(q_ref, quat_conj(q));
    error = q_err(1:3);
    if q_err(4) < 0 then
        error = -error;
    end
    omega_err = omega;

    // Recompute PID control output for logging (using stored integral error)
    m_desired = Kp.*error + Ki.*integral_error + Kd.*omega_err;

    // Limit magnetic dipole
    for j = 1:3
        if abs(m_desired(j)) > max_dipole then
            m_desired(j) = sign(m_desired(j))*max_dipole;
        end
    end

    // Project control
    B_unit = B_body/norm(B_body);
    m_parallel = (sum(m_desired.*B_unit))*B_unit;
    m_perpendicular = m_desired - m_parallel;

    // Store results
    m_history(:,i) = m_perpendicular;
    tau_history(:,i) = cross(m_perpendicular, B_body);
end

// =========================================================================
// VISUALIZATION
// =========================================================================
// Figure 1: Quaternion and Angular Velocity
scf(1); clf();
subplot(3,1,1);
plot(t/60, q_history');
xgrid();
title("Quaternion Components", "fontsize", 3);
xlabel("Time (minutes)", "fontsize", 2);
ylabel("Value", "fontsize", 2);
legend(["q_x", "q_y", "q_z", "q_w"], 2);

subplot(3,1,2);
plot(t/60, omega_history'*180/%pi);
xgrid();
title("Angular Velocity", "fontsize", 3);
xlabel("Time (minutes)", "fontsize", 2);
ylabel("Angular Rate (deg/s)", "fontsize", 2);
legend(["ω_x", "ω_y", "ω_z"], 2);

subplot(3,1,3);
plot(t/60, pointing_error, 'r');
xgrid();
title("Pointing Error", "fontsize", 3);
xlabel("Time (minutes)", "fontsize", 2);
ylabel("Error (degrees)", "fontsize", 2);

// Figure 2: Control Effort
scf(2); clf();
subplot(2,1,1);
plot(t/60, m_history');
xgrid();
title("Magnetic Dipole Moment", "fontsize", 3);
xlabel("Time (minutes)", "fontsize", 2);
ylabel("Dipole (A·m²)", "fontsize", 2);
legend(["m_x", "m_y", "m_z"], 2);

subplot(2,1,2);
plot(t/60, tau_history'*1e6);
xgrid();
title("Control Torque", "fontsize", 3);
xlabel("Time (minutes)", "fontsize", 2);
ylabel("Torque (μN·m)", "fontsize", 2);
legend(["τ_x", "τ_y", "τ_z"], 2);

// =========================================================================
// PERFORMANCE METRICS
// =========================================================================
// Calculate settling time (time to reach within 5% of final value)
threshold = 5.0; // 5 degrees
settling_indices = find(pointing_error < threshold);
if ~isempty(settling_indices) then
    settling_time = t(settling_indices(1))/60; // minutes
    disp("Settling time (to " + string(threshold) + "°): " + string(settling_time) + " minutes");
else
    disp("System did not settle below " + string(threshold) + "° threshold");
end

// Calculate steady-state error
steady_state_window = max(1, length(t) - round(100/dt)):length(t); // Last ~100 seconds
steady_state_error = mean(pointing_error(steady_state_window));
disp("Steady-state pointing error: " + string(steady_state_error) + " degrees");

// Calculate control effort
mean_dipole = mean(sqrt(sum(m_history.^2,1)));
peak_dipole = max(sqrt(sum(m_history.^2,1)));
disp("Mean magnetic dipole: " + string(mean_dipole) + " A·m²");
disp("Peak magnetic dipole: " + string(peak_dipole) + " A·m²");

disp("Simulation complete!");
