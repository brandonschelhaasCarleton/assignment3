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

% Plot V(x,y) in region
figure
p = surf(volt_2d');
set(p,'LineStyle','none')
title('Voltage Distribution within the Region');
xlabel('x [units]'); ylabel('y [units]'); zlabel('Voltage [V]');
colorbar;

% Calculate E field as -gradient of voltage (-ve in plot)
[Ex, Ey] = gradient(volt_2d');
Ex = -Ex;
Ey = -Ey;

% Plot E field as vectors
figure
quiver(Ex, Ey) % Plot Ex and Ey with Quiver function
title('E-field E(x,y)');
xlabel('x'); ylabel('y');
xlim([0 nx]); ylim([0 ny]);
rectangle('Position', [boxes{1}.x(1) boxes{1}.y(1) (boxes{1}.x(2)-boxes{1}.x(1)) (boxes{1}.y(2)-boxes{1}.y(1))]);
rectangle('Position', [boxes{2}.x(1) boxes{2}.y(1) (boxes{2}.x(2)-boxes{2}.x(1)) (boxes{2}.y(2)-boxes{2}.y(1))]);
