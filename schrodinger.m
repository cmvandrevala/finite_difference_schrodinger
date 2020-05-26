%% Prepare Workspace
% These are just a few lines to clear all variables, clear all figures, and 
% clear the workspace.

clc
clf
clear all

%% Assign Constants
% The variable hbar and the mass of the electron are set to one by default 
% (Hartree Units).  You can change them to their SI values, but you will 
% have to adjust step_size accordingly.  The variable step_size is the 
% distance between each position.  The variable number_of_data_points
% gives the total number of data points in the system.  The system spans
% from x = 0 to x = step_size*number_of_data_points.  Try to keep
% number_of_data_points less than or equal to 10000 in order to keep the 
% simulation running time reasonable.

hbar = 1;
mass = 1;
step_size = 0.05;
number_of_data_points = 1000;

%% Assign Simulation Parameters
% If you want periodic boundary conditions, mark this variable as true.
% Otherwise, leave it at its default value of false.

periodic_boundary_conditions = false;

%% Select Potential Energy
% Six different potential energy configurations are included in this 
% simulation.  You can select one of the potential energy distributions
% below by marking the variable as true.  You can mark multiple potentials
% as true and get a superimposed potential energy makeup.

finite_square_well = false;
finite_square_barrier = false;
coulomb_potential = false;
step_potential = false;
harmonic_oscillator = true;
random_potential = false;

%% Calculate Kinetic Energy Matrix
% This section uses the finite difference method to calculate the kinetic 
% energy matrix within the Hamiltonian.  It takes into account the option
% for periodic boundary conditions from above.

kinetic = zeros(number_of_data_points);

for i = 1:number_of_data_points
    kinetic(i,i) = -2;
    if i > 1
        kinetic(i, i-1) = 1;
        kinetic(i-1, i) = 1;
    end
end

if periodic_boundary_conditions
    kinetic(number_of_data_points,1) = 1;
    kinetic(1,number_of_data_points) = 1;
end

kinetic_multiplier = hbar^2/(2*mass*step_size^2);
kinetic = kinetic*kinetic_multiplier;

%% Calculate Potential Energy Matrix
% This section calculates the potential energy matrix within the
% Hamiltonian.  It creates a different matrix for each variable from the
% section "Select Potential Energy".  If you want to create your own custom
% potential energy configuration, at the end of this section, simply fill
% in the diagonal terms of the matrix 'potential' with your own custom
% values.  Each value represents the potential energy of the system at a
% point in space.

potential = zeros(number_of_data_points);

if finite_square_well
    for i = 1:0.4*number_of_data_points
        potential(i,i) = 1 + potential(i,i);
    end
    for i = 0.4*number_of_data_points+1:0.6*number_of_data_points
        potential(i,i) = 0.1 + potential(i,i);
    end
    for i = 0.6*number_of_data_points+1:number_of_data_points
        potential(i,i) = 1 + potential(i,i);
    end
end

if finite_square_barrier
    for i = 1:0.4*number_of_data_points
        potential(i,i) = 0.05 + potential(i,i);
    end
    for i = 0.4*number_of_data_points+1:0.6*number_of_data_points
        potential(i,i) = 1 + potential(i,i);
    end
    for i = 0.6*number_of_data_points+1:number_of_data_points
        potential(i,i) = 0.05 + potential(i,i);
    end
end

if coulomb_potential
    for i = 1:number_of_data_points
        potential(i,i) = 1 - 1/(i*step_size) + potential(i,i);
        if potential(i,i) < 0
            potential(i,i) = 0;
        end
    end
end

if step_potential
    for i = 1:0.5*number_of_data_points
        potential(i,i) = 0.1 + potential(i,i);
    end
    for i = 0.5*number_of_data_points + 1:number_of_data_points
        potential(i,i) = 1 + potential(i,i);
    end
end

if harmonic_oscillator
    for i = 1:number_of_data_points
        potential(i,i) = 0.0002*(i - number_of_data_points/2)^2 + potential(i,i);
        if potential(i,i) >= 5;
            potential(i,i) = 5;
        end
    end
end

if random_potential
    potential(1,1) = 0;
    for i = 2:number_of_data_points
        potential(i,i) = 0.1*(2*rand-1) + potential(i,i);
        if potential(i,i) < 0
            potential(i,i) = 0;
        end
    end
    for i = 2:number_of_data_points-1
        potential(i,i) = (potential(i-1,i-1)+2*potential(i,i)+potential(i+1,i+1))/4;
    end
end

%% Calculate Hamiltonian and Eigenvectors
% This section uses Matlab's built in eig function to calculate the
% eigenvectors and eigenvalues of the Hamiltonian.  As mentioned in section
% "Assign Constants", try to keep number_of_data_points less than or equal
% to 10000 in order to keep the simulation run time reasonable.

hamiltonian = kinetic + potential;
[eigenvectors, eigenvalues] = eig(hamiltonian);

%% Plot Eigenvectors
% This final section plots the results.  It finds the lowest three
% eigenvalues (corresponding to the ground state, first excited state, and
% second excited state) and plots each associated wavefunction.  The plot
% of each wavefunction is placed so that it is centered on its respective
% energy level (eigenvalue).  It also plots the potential energy of the 
% system as a function of position.

x_positions = linspace(0,step_size*number_of_data_points,number_of_data_points);

ground_state_wavefunction = eigenvectors(:,1) + eigenvalues(1,1)*ones(number_of_data_points,1);
first_excited_state = eigenvectors(:,2) + eigenvalues(2,2)*ones(number_of_data_points,1);
second_excited_state = eigenvectors(:,3) + eigenvalues(3,3)*ones(number_of_data_points,1);

area(x_positions, diag(potential))
hold on
plot(x_positions, ground_state_wavefunction, 'g', 'LineWidth', 2)
hold on
plot(x_positions, first_excited_state, 'c', 'LineWidth', 2)
hold on
plot(x_positions, second_excited_state, 'r', 'LineWidth', 2)
legend('Potential Energy', 'Ground State', 'First Excited State', 'Second Excited State')
xlabel('Position')
ylabel('Energy (hbar and mass set to one)')
title('Wavefunction vs. Position for the First Three Energy Levels')
