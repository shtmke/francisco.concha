function [U]=solveFractionalSystemCRL(alpha, beta,t_steps,z_steps,dt,dz,F)

a2=1/5;       % coefficient from the diffusion equation

% Number of spatial steps + 1 is:
% m = 21; % 11, 21
% 
% % Number of steps in time + 1 is:
% n =148; % 37, 148
% 
% h = L / (m-1);          % spatial step
% tau = h^2 / (6*a2);     % time step


% generating the matrix for approximation
B1 = ban(alpha,t_steps-1,dt)';       % alpha-th order derivative with respect to time
% disp('B1')
% disp(size(B1))
TD = kron(sparse(B1), speye(z_steps));          % time derivative matrix
% disp('TD')
% disp(size(TD))

B2 = ransym(beta,z_steps,dz);          % beta-th order derivative with respect to X
SD2 = kron(speye(t_steps-1), sparse(B2));        % spatial derivative matrix

B3 = ransym(alpha,z_steps,dz);          % beta-th order derivative with respect to X
SD1 = kron(speye(t_steps-1), sparse(B3));        % spatial derivative matrix

SystemMatrix = TD - a2*SD2 + SD1;   % matrix corresponding to discretization
% in space and time

% Initial conditions 
k = 1:z_steps;
U0 = 4*(k-1).*((z_steps-1) - k + 1)*dz^2;


% remove columns with '1' and 'm' from SystemMatrix
S = eliminator (z_steps, [1 z_steps]);
SK = kron(eye(t_steps-1), S);
SystemMatrix_without_columns_1_m = SystemMatrix * SK';

% remove rows with '1' and 'm' from SystemMatrix_without_columns_1_m
S = eliminator (z_steps, [1 z_steps]);
SK = kron(eye(t_steps-1), S);
SystemMatrix_without_rows_columns_1_m = SK * SystemMatrix_without_columns_1_m;
% Solution of the system
% disp('Sistema')
% disp(size(SystemMatrix_without_rows_columns_1_m))
% disp('F')
% disp(size(F))
Y = SystemMatrix_without_rows_columns_1_m\F;
% Reshape solution array -- values for k-th time step
% are in the k-th column of YS:
YS = reshape(Y,z_steps-2,t_steps-1);
YS = fliplr(YS);

for k = 1:(t_steps-1)
    U (:,k) = YS(:,k);
end



