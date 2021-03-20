set(0, 'DefaultFigureWindowStyle', 'Docked')

clc
clear
close all

% ELEC 4700 - Assignment 3
% Brandon Schelhaas
% 101036851
% Programmed on MATLAB 2020a

global const

% Add constants to the constants struct
const.m_0 = 9.10938356e-31; % kg
const.m_n = 0.26*const.m_0;
const.k = 1.38064852e-23; 
const.q = 1.60217662e-19;

doAll = 1;

% ---------------------------
%           Question 1
% ---------------------------
% Assignment 1 Parameter Definitions
T = 300;    % lattice temperature
v_th = sqrt((2 * const.k * T)/const.m_n);   % thermal velocity
tau_mn = 0.2e-12;   % mean time before collisions
lambda = tau_mn * v_th; % mean free path
regionLength = 200e-9; % nm
regionWidth = 100e-9; % nm

doPart1 = 0;
if doPart1 | doAll
    % Assignment 2 Parameter Definitions

    % a) If a voltage of 0.1V is applied across the x dimension of the semiconductor, what
    % is the electric field on the electrons? You can assume it is constant over the
    % semiconductor

    % For const. E-field, E*d = Delta V => E = V/d [V/m]
    % http://hyperphysics.phy-astr.gsu.edu/hbase/electric/elewor.html
    volt_app.x = 0.1;   % Applied Voltage in x-direction
    volt_app.y = 0;     % Applied Voltage in y-direction
    E.x = volt_app.x/regionLength;  % Calculate electric field from applied voltage - x-direction
    E.y = volt_app.y/regionWidth; % Calculate electric field from applied voltage - y-direction

    % b) What is the force on each electron?
    % F = qE
    F.x = const.q * E.x; % Calculate force from electric field - x-direction
    F.y = const.q * E.y; % Calculate force from electric field - y-direction

    % c) Calculate the acceleration on the electrons and use this in your model to update
    % the velocity of each electron at each time step. Add the capability for the electrons
    % to respond to a static electric field with both an x and a y component. 
    % F = ma --> a = F/m
    acc.x = F.x/const.m_n;
    acc.y = F.y/const.m_n;

    % d) What is the relationship between the electron drift current density and average
    % carrier velocity? Assume an electron concentration of 10^15 cm−2
    % J = e*n*v_d, e = elec charge, n = concentration, v_d = drift velocity
    n_conc = 10e15;
    % J will have to be calculated mid-simulation, because drift is average velocity

    % Simulation setup
    numElectrons = 1000; % Set the amount of electrons to simulate
    numToPlot = 10; % Set the amount of electrons to plot 
    dt = ((0.01)*regionLength)/v_th; % set a time step to move electrons 1/100 of region length
    numSteps = 200;
    % t = zeros(1, numSteps);

    % Setup initial electron positions
    pos.x = zeros(numElectrons, 2); % two columns for old and new positions
    pos.y = zeros(numElectrons, 2); % two columns for old and new positions
    pos.x(:,1) = regionLength .* rand(numElectrons, 1); % fill initial position with random val
    pos.y(:,1) = regionWidth .* rand(numElectrons, 1); % fill initial position with random val
    eColours = hsv(numToPlot); % get colours for each electron

    % Initialize electron velocity with a random thermal velocity
    vel.x = zeros(numElectrons, 1);
    vel.y = zeros(numElectrons, 1);
    vel.x(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
    vel.y(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
    %            v_th/sqrt(2) = sqrt(2kT/m)/sqrt(2) = sqrt(kT/m) = sigma

    % Calculate the probability for scattering
    p_scat = 1 - exp(-dt/tau_mn);

    % Run collisions with MFP script
    part_1
end

% ---------------------------
%           Question 2
% ---------------------------

doPart2 = 1;
if doPart2 | doAll
    volt_app.x = 0.1;   % Applied Voltage in x-direction
    volt_app.y = 0;     % Applied Voltage in y-direction

    % Simulation setup
    numElectrons = 1000; % Set the amount of electrons to simulate
    numToPlot = 10; % Set the amount of electrons to plot 
    dt = ((0.01)*regionLength)/v_th; % set a time step to move electrons 1/100 of region length
    numSteps = 200;
    
    % Set ny and nx, in nanometers, due to regionLength and regionWidth
    ny = 100;
    nx = 200;

    G = sparse(ny*ny, nx*nx);
    Bv = zeros(nx*ny, 1); % voltage boundary

    % Type 1 = specular
    % Type 2 = diffusive
    boxes{1}.x = [80e-9 120e-9];
    boxes{1}.y = [60e-9 100e-9];
    boxes{1}.sigma = 0.01;
    boxes{1}.type = 1;
    boxes{2}.x = [80e-9 120e-9];
    boxes{2}.y = [0 40e-9];
    boxes{2}.sigma = 0.01;
    boxes{2}.type = 1;

    sigma = zeros(ny, nx);

    % Run assignment part 2a code (with bottlenecks)
    part_2 % Calculates voltage and E-field

    % Have the E-field now, so re-run part 1 code with this E-field
    numToPlot = 20;
    numElectrons = 500;
    numSteps = 100;
    F.x = const.q * Ex;
    F.y = const.q * Ey;
    acc.x = F.x/const.m_n;
    acc.y = F.y/const.m_n;

    % Setup initial electron positions
    pos.x = zeros(numElectrons, 2); % two columns for old and new positions
    pos.y = zeros(numElectrons, 2); % two columns for old and new positions
    for i = 1:numElectrons
        randx = regionLength * rand();
        randy = regionWidth * rand();
        flag1 = (randx >= boxes{1}.x(1)) && (randx <= boxes{1}.x(2));
        flag2 = (randy >= boxes{1}.y(1)) || (randy <= boxes{2}.y(2));

        while (flag1 && flag2)
            randx = regionLength * rand();
            randy = regionWidth * rand();
            flag1 = (randx >= boxes{1}.x(1)) && (randx <= boxes{1}.x(2)); % check box1 and box2 x coordinate
            flag2 = (randy >= boxes{1}.y(1)) || (randy <= boxes{2}.y(2)); % check box1 and box2 y coordinate
        end

        pos.x(i,1) = randx;
        pos.y(i,1) = randy;
    end
    eColours = hsv(numToPlot); % get colours for each electron

    % Initialize electron velocity with a random thermal velocity
    vel.x = zeros(numElectrons, 1);
    vel.y = zeros(numElectrons, 1);
    vel.x(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
    vel.y(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
    %            v_th/sqrt(2) = sqrt(2kT/m)/sqrt(2) = sqrt(kT/m) = sigma

    % Calculate the probability for scattering
    p_scat = 1 - exp(-dt/tau_mn);
    
    part_2c
end

% ---------------------------
%           Question 3
% ---------------------------

% Part A - plot density, observe at 0.8V

doPart3a = 0;
if doPart3a | doAll
    volt_app.x = 1e9;
    volt_app.y = 0;

    numElectrons = 2000;

    % Type 1 = specular
    % Type 2 = diffusive
    boxes{1}.x = [80e-9 120e-9];
    boxes{1}.y = [60e-9 100e-9];
    boxes{1}.sigma = 0.01;
    boxes{1}.type = 1;
    boxes{2}.x = [80e-9 120e-9];
    boxes{2}.y = [0 40e-9];
    boxes{2}.sigma = 0.01;
    boxes{2}.type = 1;

    % Setup initial electron positions
    pos.x = zeros(numElectrons, 2); % two columns for old and new positions
    pos.y = zeros(numElectrons, 2); % two columns for old and new positions
    for i = 1:numElectrons
        randx = regionLength * rand();
        randy = regionWidth * rand();
        flag1 = (randx >= boxes{1}.x(1)) && (randx <= boxes{1}.x(2));
        flag2 = (randy >= boxes{1}.y(1)) || (randy <= boxes{2}.y(2));

        while (flag1 && flag2)
            randx = regionLength * rand();
            randy = regionWidth * rand();
            flag1 = (randx >= boxes{1}.x(1)) && (randx <= boxes{1}.x(2)); % check box1 and box2 x coordinate
            flag2 = (randy >= boxes{1}.y(1)) || (randy <= boxes{2}.y(2)); % check box1 and box2 y coordinate
        end

        pos.x(i,1) = randx;
        pos.y(i,1) = randy;
    end
    eColours = hsv(numToPlot); % get colours for each electron

    vel.x = zeros(numElectrons, 1);
    vel.y = zeros(numElectrons, 1);
    vel.x(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
    vel.y(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);

    part_3a
end

% Part B -current density vs. bottleneck width
doPart3b = 0;
if doPart3b | doAll
    numElectrons = 500;

    % Type 1 = specular
    % Type 2 = diffusive
    boxes{1}.x = [80e-9 120e-9];
    boxes{1}.sigma = 0.01;
    boxes{1}.type = 1;
    boxes{2}.x = [80e-9 120e-9];
    boxes{2}.sigma = 0.01;
    boxes{2}.type = 1;

    for count = 1:2
        % Set different bottleneck widths
        if count == 1
            figure
            boxes{1}.y = [70e-9 100e-9];
            boxes{2}.y = [0 30e-9];
        else
            boxes{1}.y = [80e-9 100e-9];
            boxes{2}.y = [0 20e-9];
        end
        percent_open = (boxes{1}.y(1) - boxes{2}.y(2))/(regionWidth) * 100;

        % Setup initial electron positions
        pos.x = zeros(numElectrons, 2); % two columns for old and new positions
        pos.y = zeros(numElectrons, 2); % two columns for old and new positions
        for i = 1:numElectrons
            randx = regionLength * rand();
            randy = regionWidth * rand();
            flag1 = (randx >= boxes{1}.x(1)) && (randx <= boxes{1}.x(2));
            flag2 = (randy >= boxes{1}.y(1)) || (randy <= boxes{2}.y(2));

            while (flag1 && flag2)
                randx = regionLength * rand();
                randy = regionWidth * rand();
                flag1 = (randx >= boxes{1}.x(1)) && (randx <= boxes{1}.x(2)); % check box1 and box2 x coordinate
                flag2 = (randy >= boxes{1}.y(1)) || (randy <= boxes{2}.y(2)); % check box1 and box2 y coordinate
            end

            pos.x(i,1) = randx;
            pos.y(i,1) = randy;
        end
        eColours = hsv(numToPlot); % get colours for each electron

        % Initialize velocities
        vel.x = zeros(numElectrons, 1);
        vel.y = zeros(numElectrons, 1);
        vel.x(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
        vel.y(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);

        part_3b
    end
end