%velfluc.m 
%plotting the streamwise fluctuation streaks
% Define independent variable
y = chebfun('y');

% Define the parameters
Re = 3250; Re_tau = 182; % Reynolds number
kxval = 1.13*[1 -1];        % streamwise wave-number
omval  = 0.113*[-1 1]; % temporal frequency

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



% Start computations
for n = 1:2

    omega = omval(n);
    kx = kxval(n);
    kx2 = kx*kx;
    kx4 = kx2*kx2;

    % Define operator A
    A = chebop([-1 1]);
    a2 = -( 2*kx^2/Re + 1i*kx*U + 1i*omega );
    a0 = kx^4/Re + 1i*kx^3*U + 1i*kx*Uyy + 1i*omega*kx^2;
    A.op = @(y,phi) (diff(phi,4)/Re + a2*diff(phi,2) + a0*phi);

    % Specify boundary conditions
    A.lbc = @(phi) [phi;diff(phi)];
    A.rbc = @(phi) [phi;diff(phi)];

    % B operator
    B = chebop([-1 1]);
    B.op = @(y,d1,d2) (-diff(d1) + 1i*kx*d2);

    % C operator
    C = chebop([-1 1]);
    C.op = @(y,phi) [diff(phi); -1i*kx*phi];

    % Solving for the principal singular value
    [PhiAndPsi,sval] = svdfr(A,B,C,1); % change if '1' to '2' if you want the second resolvent mode
    uandv = C(PhiAndPsi.blocks{1,1});   % First variable is the regular variable, phi.
                                        % Note that C(Phi) gives the output, [u;v]
    u = uandv{1};                       % streamwise velocity
    v = uandv{2};                       % wall-normal velocity

    % discretized values for plotting
    uvec(:,n) = u(yd,1);
    vvec(:,n) = v(yd,1);

end


% Getting physical fields of u and v
kx = abs(kxval(1));
xval = linspace(0, 4*pi/kx, 100); % streamwise coordinate

Up = zeros(N,length(xval));   % physical value of u
Vp = zeros(N,length(xval));   % physical value of v

for indx = 1:length(xval)

    x = xval(indx);

    for n = 1:2

        kx = kxval(n);

        Up(:,indx) = Up(:,indx) + uvec(:,n)*exp(1i*kx*x);
        Vp(:,indx) = Vp(:,indx) + vvec(:,n)*exp(1i*kx*x);

    end

end

Up = real(Up); Vp = real(Vp); % extract real parts

xval=0.125*xval*Re_tau;
yd = flip(yd)*Re_tau+Re_tau;

Up = Up./max(Up,[],'all').*3;

% Plotting the most amplified streamwise velocity structures
pcolor(xval,yd,Up); shading interp;
cb = colorbar('vert');
xlabel('z^+');
ylabel('y^+');
ax = gca;
cb.Ticks = -3:1:3;
colormap turbo;
xlim([0 150])
ylim([0 150])



