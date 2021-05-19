%singval.m
%computing the largest singular values of the flow 
% Define independent variable
y = chebfun('y');

% Define the parameters
Re = 3250; Re_tau = 182; % Reynolds number
kxval = 1.13;        % streamwise wave-number
omega  = -0.113; % temporal frequency

N = 100;               % number of Chebyshev points for plotting
yd = chebpts(N);

% Define base flow
data = importdata('re180.txt');  
yp = flip(data(:,1)); up = data(:,3)./1.8283888e+01;
U = zeros(100,1);
for i=1:100
    glp = abs(yd(i));
    glv = interp1(yp,up,glp);
    U(i) = glv;
end 

U = chebfun(U);
Uy = diff(U);
Uyy = diff(U,2);

% Looping over different values of kx

kx = kxval; kx2 = kx*kx; kx4 = kx2*kx2;

% Define operator A
A = chebop([-1 1]);
a2 =  -( 2*kx^2/Re  +  1i*kx*U + 1i*omega );
a0 =  kx^4/Re  +  1i*kx^3*U +  1i*kx*Uyy + 1i*omega*kx^2;
A.op = @(y,phi) (diff(phi,4)/Re + a2*diff(phi,2) + a0*phi);

% Specify boundary conditions
A.lbc = @(phi) [phi;diff(phi)];
A.rbc = @(phi) [phi;diff(phi)];

% B operator
B = chebop([-1 1]);
B.op = @(y,d1,d2) (-diff(d1)  +  1i*kx*d2);

% C operator
C = chebop([-1 1]);
C.op = @(y,phi) [diff(phi);-1i*kx*phi];

% Solving for the principal singular value
svals = abs(svdfr(A,B,C,20,'LM'));   % Compute the 20 largest singular values

svals = flip(svals);


plot(svals, 'ob','MarkerSize',8)
xlabel('j')
ylabel('\sigma_j')