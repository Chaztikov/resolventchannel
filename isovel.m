%isovel.m
% Set parameters and base velocity
ii = sqrt(-1);      % Imaginary unit
kx = 1.13; kz = 1.13;     % Spatial wavenumbers
k2 = kx*kx + kz*kz;
k4 = k2*k2;
Re = 3250; Re_tau = 186;                 % Reynolds number
kxval = 1.13*[1 1 -1 -1];        % streamwise wavenumber
kzval = 1.13*[1 -1 1 -1];        % spanwise wavenumber
omval  = 0.113*[-1 -1 1 1]; % temporal frequency
y = chebfun('y');

U = 1 - y^2;                % base flow (Poiseuille; U = y for Couette)
Uy = diff(U);               % derivative of U
Uyy = diff(U,2);            % second derivative of U

% Looping over kx, kz, and om
N = 200;                    % Grid size for plotting in y-direction
uvec = zeros(N,length(kxval));
dzvec = zeros(N,length(kxval));

for n = 1:length(kxval)

    kx = kxval(n); kz = kzval(n); omega = omval(n);
    kx2 = kx*kx; kz2 = kz*kz;
    k2 = kx2 + kz2;

    A = chebop([-1,1]);     % Operator A
    B = chebop([-1,1]);     % Operator B
    C = chebop([-1,1]);     % Operator C

    A.op = @(x,v,eta)([1i*omega*(diff(v,2) - k2*v) - (diff(v,4)-2*k2*diff(v,2) ...
            + k4*v)/Re - 1i*kx*Uyy*v  + 1i*kx*U*(diff(v,2) - k2*v);...
            1i*omega*eta + 1i*kz*Uy*v + 1i*kx*U*eta - (diff(eta,2) - k2*eta)/Re]);

    A.lbc = @(v,eta)[diff(v);v;eta];
    A.rbc = @(v,eta)[diff(v);v;eta];

    B.op = @(x,dx,dy,dz) ([-1i*kx*diff(dx) - k2*dy - 1i*kz*diff(dz);...
                            1i*kz*dx - 1i*kx*dz]);
    C.op = @(x,v,eta)([1i*kx*diff(v)/k2 - 1i*kz*eta/k2;...
                    v ; ...
                    1i*kz*diff(v)/k2 + 1i*kx*eta/k2]);

    % Compute the singular function
    [PhiAndPsi,sval] = svdfr(A,B,C,1);

    % velocities
    uvw = C(PhiAndPsi(1:2,:));  % First two variables are the regular variables,
                                % v and eta, so that Phi = [v;eta].
                                % Note C(Phi) gives the output, [u;v;w]
    u = uvw.blocks{1};          % streamwise velocity
    v = uvw.blocks{2};          % wall-normal velocity
    w = uvw.blocks{3};          % spanwise velocity

    % Body forces:
    Bad = adjointNS(B);             % The body forces are computed as Bad(Psi);
    dxdydz = Bad(PhiAndPsi(3:4,:)); % The second two arguments are Psi, the
                                    % auxiliary variables for the adjoint system
    dx = dxdydz.blocks{1};
    dy = dxdydz.blocks{2};
    dz = dxdydz.blocks{3};

    uvec(:,n) = u(chebpts(N));
    dzvec(:,n) = dz(chebpts(N));

end

kx = abs(kxval(1)); kz = abs(kzval(1));
zval = linspace(-7.8, 7.8, 100);    % spanwise coordinate
xval = linspace(0, 12.7, 100);      % streamwise coordinate

Up = zeros(length(zval),length(xval),N);    % physical value of u
Dp = zeros(length(zval),length(xval),N);    % physical value of dz

for indx = 1:length(xval)
    x = xval(indx);
    for indz = 1:length(zval)
        z = zval(indz);
        for n = 1:4

            kx = kxval(n); kz = kzval(n);

            Up(indz,indx,:) =  squeeze(Up(indz,indx,:)) + ...
                uvec(:,n)*exp(1i*kx*x + 1i*kz*z);

            Dp(indz,indx,:) =  squeeze(Dp(indz,indx,:)) + ...
                dzvec(:,n)*exp(1i*kx*x + 1i*kz*z);
        end
    end
end
Up = real(Up); Dp = real(Dp); % only real part exist

%% Plot isosurfaces of optimal flow structures
clf
Upmax = max(max(max(Up)));
val = 0.2*Upmax;
ypt = chebpts(N);

ypt=ypt*Re_tau+180;
xval=xval*Re_tau;
zval=zval*Re_tau+120;



[X,Z,Y] = meshgrid(xval,zval,ypt);
p = patch(isosurface(X,Z,Y,Up,val));
hh = isonormals(xval,zval,ypt,Up,p);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect('auto')
% view(3);
view([-15 30])
ax = gca;
xlabel('x^+');
ylabel('z^+');
% ax.YTick = -10:5:10;
zlabel('y^+');
% ax.ZTick = -1:0.5:1;
zlim([0 200])
xlim([0 2000])
ylim([-500 500])
camlight

hold on
p = patch(isosurface(X,Z,Y,Up,-val));
isonormals(xval,zval,ypt,Up,p);
p.FaceColor = 'blue';
p.EdgeColor = 'none';
disp('done');
hold off
