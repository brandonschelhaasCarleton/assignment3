% Boundary Voltage Conditions
leftBC = volt_app.x;
rightBC = 0;

% Setup conductivities for each mesh element
sigmaOutsideBox = 1;
for col = 1:nx
    for row = 1:ny 
        if (col*1e-9 >= boxes{1}.x(1) & col*1e-9 <= boxes{1}.x(2) & row*1e-9 >= boxes{1}.y(1) & row*1e-9 <= boxes{1}.y(2))
            sigma(row,col) = boxes{1}.sigma;
        elseif (col*1e-9 >= boxes{2}.x(1) & col*1e-9 <= boxes{2}.x(2) & row*1e-9 >= boxes{2}.y(1) & row*1e-9 <= boxes{2}.y(2))
            sigma(row,col) = boxes{2}.sigma;
        else
            sigma(row,col) = sigmaOutsideBox;
        end
    end
end

% Setup G-Matrix (resistor mesh)
for col = 1:nx
    for row = 1:ny
        n = row + (col-1)*ny;
        nxm = row + (col-2)*ny;
        nxp = row + (col)*ny;
        nym = (row-1) + (col-1)*ny;
        nyp = (row+1) + (col-1)*ny;
        
        if col == 1
            % V(n) = Boundary Voltage = 1V
            G(n,n) = 1;
            Bv(n) = leftBC;
        elseif col == nx
            % V(n) = Boundary Voltage = 1V
            G(n,n) = 1;
            Bv(n) = rightBC;
        elseif row == 1
            % Calculate resistors between nodes
            rxm = (sigma(row,col) + sigma(row,col-1)) / 2.0;
            rxp = (sigma(row,col) + sigma(row,col+1)) / 2.0;
            ryp = (sigma(row,col) + sigma(row+1,col)) / 2.0;
            
            % Change G-matrix to use different conductivities
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        elseif row == ny
            % Calculate resistors between nodes
            rxm = (sigma(row,col) + sigma(row,col-1)) / 2.0;
            rxp = (sigma(row,col) + sigma(row,col+1)) / 2.0;
            rym = (sigma(row,col) + sigma(row-1,col)) / 2.0;

            % Change G-matrix to use different conductivities
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            % Calculate resistors between nodes
            rxm = (sigma(row,col) + sigma(row,col-1)) / 2.0;
            rxp = (sigma(row,col) + sigma(row,col+1)) / 2.0;
            rym = (sigma(row,col) + sigma(row-1,col)) / 2.0;
            ryp = (sigma(row,col) + sigma(row+1,col)) / 2.0;

            % Change G-matrix to use different conductivities
            G(n,n) = -(rxm+rxp+rym+ryp); %FDM implementation (-4V(n))
            G(n,nxm) = rxm; %(+1V(nxm))...
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end
    end
end

%Calculate Voltage in region
V = G\Bv; %GV = F ... GV = Bv

volt_2d = zeros(nx, ny);
for col = 1:nx
    for row = 1:ny
        n = row + (col-1)*ny;
        volt_2d(col, row) = V(n);
    end
end

% Calculate E field as -gradient of voltage (-ve in plot)
[Ex, Ey] = gradient(volt_2d');
Ex = -Ex;
Ey = -Ey;

% Calculate force and acceleration from E-field
F.x = const.q * Ex;
F.y = const.q * Ey;
acc.x = F.x/const.m_n;
acc.y = F.y/const.m_n;

% Dynamics simulation
J_sum = 0;
for i = 1:numSteps
    % Check for scattering
    randNums = rand(numElectrons,1);
    scattered = find(randNums(:,1) < p_scat);
    if(length(scattered)) %if at least one position is beyond region length
        for k = 1:length(scattered)
            vel.x(scattered(k),1) = v_th/sqrt(2) * randn(); % assign new velocity
            vel.y(scattered(k),1) = v_th/sqrt(2) * randn();
        end
    end
    
    % To update velocities and positions, find electrons within each mesh element and apply the force to those electrons
    % Necessary since electric field may be different for each mesh element, so different forces
    for col = 1:nx
       xhigh = col * 1e-9;
       xlow = xhigh - 1e-9;
       for row = 1:ny
           yhigh = row * 1e-9;
           ylow = yhigh - 1e-9;
           elecsInElement = pos.x(:,1) >= xlow & pos.x(:,1) <= xhigh & pos.y(:,1) >= ylow & pos.y(:,1) <= yhigh;
           vel.x(elecsInElement) = vel.x(elecsInElement) + acc.x(row, col) * dt;
           vel.y(elecsInElement) = vel.y(elecsInElement) + acc.y(row, col) * dt;
           pos.x(elecsInElement,2) = pos.x(elecsInElement,1) + vel.x(elecsInElement,1) * dt + (1/2)*acc.x(row, col)*dt*dt;
           pos.y(elecsInElement,2) = pos.y(elecsInElement,1) + vel.y(elecsInElement,1) * dt + (1/2)*acc.y(row, col)*dt*dt;
       end
    end
    
    % Accommodate x=regionLength BC
    pastRightBC = find(pos.x(:,2) > regionLength); % element = 1 if beyond BC
    if(length(pastRightBC)) %if at least one position is beyond region length
        for k = 1:length(pastRightBC)
            pos.x(pastRightBC(k),1) = pos.x(pastRightBC(k),1) - regionLength;
            pos.x(pastRightBC(k),2) = pos.x(pastRightBC(k),2) - regionLength;
        end
    end
    
    % Accommodate x=0 BC
    pastLeftBC = find(pos.x(:,2) < 0); % element = 1 if beyond BC
    if(length(pastLeftBC)) %if at least one position is beyond region length
        for k = 1:length(pastLeftBC)
            pos.x(pastLeftBC(k),1) = pos.x(pastLeftBC(k),1) + regionLength;
            pos.x(pastLeftBC(k),2) = pos.x(pastLeftBC(k),2) + regionLength;
        end
    end
    
    % Accommodate y=0 or y=regionWidth BC's
    pastTopOrBottomBC = find((pos.y(:,2) > regionWidth) | (pos.y(:,2) < 0)); % find electrons past top/bottom BC
    if(length(pastTopOrBottomBC)) %if at least one position is beyond region length
        for k = 1:length(pastTopOrBottomBC)
            vel.y(pastTopOrBottomBC(k)) = -1 * vel.y(pastTopOrBottomBC(k)); %angle in = angle out
        end
    end
    
    % Accommodate the boxes' BC's (hard-coded boundaries are easier)
    boxesCheck = find( ((pos.x(:,2) > boxes{1}.x(1)) & (pos.x(:,2) < boxes{1}.x(2))) & ((pos.y(:,2) < boxes{2}.y(2)) | (pos.y(:,2) > boxes{1}.y(1))));
    if(length(boxesCheck))
        for k = 1:length(boxesCheck)
            if(pos.y(boxesCheck(k),2) >= boxes{1}.y(1)) % Inside box 1
                if(boxes{1}.type == 1) % Specular Case
                    % Left side
                    if(pos.x(boxesCheck(k),1) < 80e-9)
                        vel.x(boxesCheck(k),1) = -1 * vel.x(boxesCheck(k),1);
                    % Right side
                    elseif(pos.x(boxesCheck(k),1) > 120e-9)
                        vel.x(boxesCheck(k),1) = -1 * vel.x(boxesCheck(k),1);
                    % Bottom
                    else
                        vel.y(boxesCheck(k),1) = -1 * vel.y(boxesCheck(k),1);
                    end
                else % Diffusive case
                    if(pos.x(boxesCheck(k),1) < 80e-9)
                        % Only way to hit left wall is with positive x vel
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.x(boxesCheck(k),1) > 0)
                            % Ensure new x velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate negative random vel
                        end
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    % Right side
                    elseif(pos.x(boxesCheck(k),1) > 120e-9)
                        % Must have negative x vel to hit the right side
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.x(boxesCheck(k),1) < 0)
                            % Ensure new x velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate positive random vel
                        end
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    % Bottom
                    else
                        % Must have positive y vel to hit the bottom side
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.y(boxesCheck(k),1) > 0)
                            % Ensure new y velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate negative random vel
                        end
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    end
                end
            else % Inside box 2
                if(boxes{2}.type == 1) % Specular Case
                    % Left side
                    if(pos.x(boxesCheck(k),1) < 80e-9)
                        vel.x(boxesCheck(k),1) = -1 * vel.x(boxesCheck(k),1);
                    % Right side
                    elseif(pos.x(boxesCheck(k),1) > 120e-9)
                        vel.x(boxesCheck(k),1) = -1 * vel.x(boxesCheck(k),1);
                    % Top
                    else
                        vel.y(boxesCheck(k),1) = -1 * vel.y(boxesCheck(k),1);
                    end
                else % Diffusive case
                    if(pos.x(boxesCheck(k),1) < 80e-9)
                        % Only way to hit left wall is with positive x vel
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.x(boxesCheck(k),1) > 0)
                            % Ensure new x velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate negative random vel
                        end
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    % Right side
                    elseif(pos.x(boxesCheck(k),1) > 120e-9)
                        % Must have negative x vel to hit the right side
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.x(boxesCheck(k),1) < 0)
                            % Ensure new x velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate positive random vel
                        end
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    % Bottom
                    else
                        % Must have positive y vel to hit the bottom side
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.y(boxesCheck(k),1) > 0)
                            % Ensure new y velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate negative random vel
                        end
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    end
                end
            end
        end
    end
    
    % Re-update positions to ensure new velocities are used for new pos's
    for col = 1:nx
       xhigh = col * 1e-9;
       xlow = xhigh - 1e-9;
       for row = 1:ny
           yhigh = row * 1e-9;
           ylow = yhigh - 1e-9;
           elecsInElement = pos.x(:,1) >= xlow & pos.x(:,1) <= xhigh & pos.y(:,1) >= ylow & pos.y(:,1) <= yhigh;
           pos.x(elecsInElement,2) = pos.x(elecsInElement,1) + vel.x(elecsInElement,1) * dt + (1/2)*acc.x(row, col)*dt*dt;
           pos.y(elecsInElement,2) = pos.y(elecsInElement,1) + vel.y(elecsInElement,1) * dt + (1/2)*acc.y(row, col)*dt*dt;
       end
    end
    
    % Calculate current density
    v_dx = mean(vel.x);
    J.x = const.q * n_conc * v_dx; % only need Jx because volt_app.y = 0
    J_sum = J_sum + J.x;
    J_average = J_sum/i;
    
    % Plot current vs. bottle neck width
    subplot(2,1,count)
    plot(i, J_average, 'b.');
    if i == 1
        title(['Average Current | Bottleneck = ', num2str(percent_open), '% Open']);
        xlabel('Step'); ylabel('Current [A/cm^2]');
    end
    hold on;
    
    % Old positions now equal new positions
    pos.x(:,1) = pos.x(:,2);
    pos.y(:,1) = pos.y(:,2);
    
    pause(0.001)
end