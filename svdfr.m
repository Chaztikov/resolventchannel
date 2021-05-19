function varargout = svdfr(A,varargin)
%% SVDFR
% Computes singular values corresponding to frequency responses
% For linear block matrix operator of the form,
%
% Needs adjointFormal.m of Chebfun as a separate function (provided with this code).
%%
%
% $$\partial_t M u(y,t) = L u(y,t) +  B d(y)e^{1i \omega t}$$
% $$\eta(y) = C u(y,t),$$
%
% where, L,M,B and C are block matrix operators in $y$, L and M are square.
% The solution is given by,
%%
%
% $$\eta(y,t) = (B*inv(1i \omega M - L)*C)d(y)e^{1i \omega t},$$
% hence,
%%
%
% $$\eta(y) = (B*inv(1i \omega M - L)*C)d(y),$$
%
% The singular values of (B*inv(1i \omega M - L)*C) is computed by this
% function.
% Physically this represents the maximum amplitude of the output for a
% given sinusoidal input, e^{1i \omega t}. The right singular functions
% \eta(y) left singular function d(y) give shape of the output and
% input corresponding to the maximum amplification.
%
% Usage:
%
% A: (i \omega M) - L, A is of type chebop. Boundary conditions have to be
% specified on this operator as A.lbc and A.rbc. Other types of Chebfun
% boundary conditions like A.bc will not work (will consider in a future
% release).
%
% B: B is of type chebop
%
% C: C is of type chebop
%
% Usage is similar to chebop/eigs.m, most of these comments on usage are
% taken from chebop/eigs.m.
%
% D = SVDFR(A, B, C) returns a vector of 6 largest singular values 
% of the system represented by A, B and C.
%
% [V, D] = SVDFR(A,B,C) returns a diagonal 6x6 matrix D of A's largest
% singular values, and a quasimatrix V of the corresponding both left and 
% right singular functions.
%
% SVDFR(A, B, C, K) for an integer K > 0 finds the K largest
% singular values.
%
% SVDFR(A, B, C, K, SIGMA) finds K singular values. If SIGMA is a
% scalar, the singular values found are the ones closest to SIGMA. Other
% possibilities are 'LR' and 'SR' for the eigenvalues of largest and smallest
% real part, and 'LM' (or Inf) and 'SM' for largest and smallest magnitude.
% SIGMA must be chosen appropriately for the given operator; for example, 'LM'
% for an unbounded operator will fail to converge! Use flag 'inverse' to
% solve for the inverse of singular values using the 'SM' option.
%
% SVDFR(..., PREFS) accepts a CHEBOPPREF to control the behavior of the
% algorithm. If empty, defaults are used.
%
% Written by Gokul, harih020@umn.edu


%cheboppref.setDefaults('discretization',@ultraS)
[lbc2,rbc2] = adjointBcEndpoints(A);

B = varargin{1};
C = varargin{2};
%lbc2 = varargin{3};
%rbc2 = varargin{4};

% Determine size of lbcs and rbcs
Atemp = A;
Atemp.rbc = [];
Atempop = linearize(Atemp,0*Atemp.init);
slbc = size(Atempop.constraint.values);
Atemp = A;
Atemp.lbc = [];
Atempop = linearize(Atemp,0*Atemp.init);
srbc = size(Atempop.constraint.values);
Aop = linearize(A,0*A.init);
Bop = linearize(B,0*B.init);
Cop = linearize(C,0*C.init);
[nrows,ncols]=size(Aop);
for ii = 1:nrows
    for jj = 1:ncols
        if ( isa(Aop.blocks{ii,jj}, 'chebfun') )
            b = Aop.blocks{ii,jj};
            Aop.blocks{ii,jj} = operatorBlock.mult(b, b.domain);
        end
    end
end
[nrows,ncols]=size(Bop);
for ii = 1:nrows
    for jj = 1:ncols
        if ( isa(Bop.blocks{ii,jj}, 'chebfun') )
            b = Bop.blocks{ii,jj};
            Bop.blocks{ii,jj} = operatorBlock.mult(b, b.domain);
        end
    end
end
[nrows,ncols]=size(Cop);
for ii = 1:nrows
    for jj = 1:ncols
        if ( isa(Cop.blocks{ii,jj}, 'chebfun') )
            b = Cop.blocks{ii,jj};
            Cop.blocks{ii,jj} = operatorBlock.mult(b, b.domain);
        end
    end
end


pref = chebfunpref();
Astarop = adjointFormal(Aop,pref);
Bstarop = adjointFormal(Bop,pref);
Cstarop = adjointFormal(Cop,pref);

BBstar = linop(Bop*Bstarop);
CstarC = linop(Cstarop*Cop);

Zl = chebop(Aop.domain);
Zl.op = @(x,u) (0*u);
Zl = linearize(Zl,0*Zl.init);


Zlinop = linop(zeros(size(Aop)));
[nrows, ncols] = size(Aop);
for i = 0:nrows-1
    for j = 0:ncols-1
        Zlinop.blocks{i+1,j+1} = Zl.blocks{1};
    end
end
Zlinop.domain = Aop.domain;
M = [Aop, Zlinop;
    Zlinop, Astarop];
L = [Zlinop,BBstar;CstarC,Zlinop];

% Need to manipulate the strings so that the
% adjoint and regular boundary conditions have different names.

% new names for adjoint operator:
varlist = cell(ncols,1);
varAdlist = cell(ncols,1);
for i = 1:ncols
    varlist{i} = ['v' int2str(i)];
    varAdlist{i} = ['w' int2str(i)];
end

rbc1 = func2str(A.rbc);
[rbc1,varlist] = replaceNames(rbc1,varlist);
lbc1 = func2str(A.lbc);
[lbc1,varlist] = replaceNames(lbc1,varlist);
rbc2 = func2str(rbc2);
[rbc2,varAdlist] = replaceNames(rbc2,varAdlist);
lbc2 = func2str(lbc2);
[lbc2,varAdlist] = replaceNames(lbc2,varAdlist);

if (srbc(1)>1)
    rbcnew = strtok(rbc1,')');
    rbcnew = [rbcnew ',' strtok(rbc2(3:end),')') ')' ];
    [~,tempchar] = strtok(rbc1,'[');
    rbcnew = [rbcnew tempchar];
    rbcnew(end) = [];
    rbcnew = [rbcnew ';'];
    [~,tempchar] = strtok(rbc2,'[');
    tempchar(1) = [];
    rbcnew = [rbcnew tempchar ';'];
else
   rbcnew = strtok(rbc1,')');
    rbcnew = [rbcnew ',' strtok(rbc2(3:end),')') ')' ];
    [~,tempchar] = strtok(rbc1,')');
    rbcnew = [rbcnew '[' tempchar(2:end)];
    %rbcnew(end) = [];
    rbcnew = [rbcnew ';'];
    [~,tempchar] = strtok(rbc2,'[');
    tempchar(1) = [];
    rbcnew = [rbcnew tempchar ';'];
end

if (slbc(1)>1)
    lbcnew = strtok(lbc1,')');
    lbcnew = [lbcnew ',' strtok(lbc2(3:end),')') ')' ];
    [~,tempchar] = strtok(lbc1,'[');
    lbcnew = [lbcnew tempchar];
    lbcnew(end) = [];
    lbcnew = [lbcnew ';'];
    [~,tempchar] = strtok(lbc2,'[');
    tempchar(1) = [];
    lbcnew = [lbcnew tempchar ';'];
else
    lbcnew = strtok(lbc1,')');
    lbcnew = [lbcnew ',' strtok(lbc2(3:end),')') ')' ];
    [~,tempchar] = strtok(lbc1,')');
    lbcnew = [lbcnew '[' tempchar(2:end)];
    %rbcnew(end) = [];
    lbcnew = [lbcnew ';'];
    [~,tempchar] = strtok(lbc2,'[');
    tempchar(1) = [];
    lbcnew = [lbcnew tempchar ';'];
end
eval(['rbcnew = ', rbcnew]);
eval(['lbcnew = ', lbcnew]);

orderOfA = Aop.diffOrder;
orderOfAstar = orderOfA';
dummy = chebop(Aop.domain);
strdummy = '@(x,';
for i = 1:length(varlist)
    strdummy = [strdummy varlist{i} ','];
end
for i = 1:length(varAdlist)
    strdummy = [strdummy varAdlist{i} ','];
end
strdummy(end) = [];
strdummy = [strdummy ')['];

for i = 1:length(varlist)
for j = 1:length(varlist)
   strdummy =  [strdummy 'diff(' varlist{j} ',orderOfA(' int2str(i)...
                ',' int2str(j) '))+'];
end
strdummy(end) = [];
strdummy = [strdummy ';' newline];
end

for i = 1:length(varAdlist)
for j = 1:length(varAdlist)
   strdummy =  [strdummy 'diff(' varAdlist{j} ',orderOfAstar(' int2str(i)...
                ',' int2str(j) '))+'];
end
strdummy(end) = [];
strdummy = [strdummy ';' newline];
end
strdummy(end) = [];
strdummy = [strdummy '];'];

eval(['dummy.op = ', strdummy]);
dummy.lbc = lbcnew;
dummy.rbc = rbcnew;
linCheck = true;
dummy = linearize(dummy, dummy.init, [], linCheck);



sigflag1 = false;

if (numel(varargin) == 2)
    numval = 6;
    sigflag1 = true;
    pref = cheboppref();
    sigma = 'LR';
    pref.dicretization = @ultraS;
    pref = determineDiscretization(A, length(A.domain), pref);
elseif (numel(varargin) == 3)
    numval = varargin{3};
    sigflag = true;
    sigma = 'LR';
    pref = cheboppref();
    pref = determineDiscretization(A, length(A.domain), pref);
    pref.dicretization = @ultraS;
elseif (numel(varargin)>3 && numel(varargin)<5)
    if (ischar(varargin{4}) || isnumeric(varargin{4}))
        numval = varargin{3};
        sigma = varargin{4};
        pref = cheboppref();
        pref = determineDiscretization(A, length(A.domain), pref);
        pref.dicretization = @ultraS;
    else
        sigma = 'LR';
        numval = varargin{3};
        pref = varargin{4};
        pref = determineDiscretization(A, length(A.domain), pref);
    end

elseif (numel(varargin)>3 && numel(varargin)<6)
    numval = varargin{3};
    sigflag = false;
    sigma = varargin{4};
    pref = varargin{5};
    pref = determineDiscretization(A, length(A.domain), pref);
else
    error('Wrong number of arguments for svdfr');
end
if strcmp(sigma,'inverse')
    M.constraint = dummy.constraint;
else
    L.constraint = dummy.constraint;
end


%     pref = varargin{6};
%     pref = determineDiscretization(A, length(A.domain), pref);
% else
%     pref = cheboppref();
%     pref = determineDiscretization(A, length(A.domain), pref);
% end


if strcmp(sigma,'inverse')
    [V,lam] = eigs(M,L,numval,'SM',pref);
else
    [V,lam] = eigs(L,M,numval,sigma, pref);
end

if ( nargout < 2 )  % Return the eigenvalues only
    varargout = { diag(lam) };
else
    varargout = {V,lam};
end

end

function [out,newerNames] = replaceNames(in,newName)
%   Gives new names to variables in function handles.

% Make a cell list of names
ne = split(in,'(');
ne1 = split(ne{2},')');
ne2 = split(ne1{1},',');

out = in;
for i = 1:length(newName)
    out = replace(out,ne2{i},newName{i});
end
ne = split(out,'(');
ne1 = split(ne{2},')');
newerNames = split(ne1{1},',');
end

function [lbcopH,rbcopH] = adjointBcEndpoints(A)
%   Computes adjoint boundary conditions when boundary conditions are
%   specified at the end points. The input A must be of type chebop with
%   boundary conditions specified as A.lbc and A.rbc for this to work. For
%   periodic boundary conditions, use adjointBcPeriodic(A).

% Linearize A:
L = linearize(A,0*A.init,true);
[nrows,ncols] = size(L);

if (nrows~=ncols)
    error('This routine is only for square matrices')
end

orderA = L.diffOrder;
n = max(max(orderA));
orderAstar = orderA';
hh = max(orderA');

% Order of boundary conditions have to one less.
orderAbc = orderA - ones(size(orderA));
orderAbc(orderAbc<0) = 0;
orderAstarbc = orderAstar - ones(size(orderAstar));
orderAstarbc(orderAstarbc<0) = 0;
% Make n + 1 chebmatrices for alphais.
alphais = cell(n+1,1);
for i = 1:n+1
    alphais{i} = chebmatrix(linop(zeros(size(L))));
end

% Load coeffs, in the form alpha{0+1} + alpha{1 +1} D1 + alpha{1 + 2} D2
% etc.


for i = 1:nrows
    for j = 1:ncols
        tempMat = L.blocks{i,j}.toCoeff;
        lentempMat = length(tempMat);
        for k = 1:n+1
            if (k <= lentempMat)
                alphais{k}(i,j) = tempMat.blocks{lentempMat - k + 1};
            else
                alphais{k}(i,j) = chebfun(0,L.domain);
            end
        end
    end
end
% make a chebmatrix with zero chebfuns:
zchebmat = chebmatrix(linop(zeros(size(L))));
for i = 1:nrows
    for j = 1:ncols
        zchebmat(i,j) = chebfun(0,L.domain);
    end
end
% Compute Amat based on alphais, based on Cij in theorem.
Amat = cell(n,n);
for i = 1:n
    for j = 1:n
        Amat{i,j} = zchebmat;
    end
end

for i = 0:n-1
    for j = 0:n-i-1
        for k = i:n-j-1
           Amat{i + 1,j + 1} =  Amat{i + 1,j + 1} + (-1)^k * nchoosek(k, i) .* ...
               diff(alphais{k + j + 2},k-i);
        end
    end
end

% Evaluate Amat at the endpoints to get Aplus and Amin:
Aplus = zeros(size(Amat));
Amin = zeros(size(Amat));

for i = 0:n-1
    for j = 0:n-1
        for l = 0:nrows-1
            for m=0:ncols-1
                Aplus(i*nrows + l + 1, j*ncols +m+1) = Amat{i+1,j+1}.blocks{l+1,m+1}(L.domain(1));
                Amin(i*nrows + l + 1, j*ncols +m+1) = Amat{i+1,j+1}.blocks{l+1,m+1}(L.domain(2));
            end
        end
    end
end

% Parse boundary conditions.
% Gokul: This part can probably be done more efficiently. I convert A.lbc
% or A.rbc to a
% chebop, then linearize it to get an operator that gives the coefficients
% as chebfun and then evaluate at the endpoints.


% variables in the lbc:
ne = split(func2str(A.lbc),'(');
ne1 = split(ne{2},')');
ne2l = split(ne1{1},',');

slbc = nargin(A.lbc);
if (slbc(1) > 1)
me = split(func2str(A.lbc),'[');
me1 = split(me{2},']');
else
    catexp = '@\(\w*\)';
    me = regexp(func2str(A.lbc),catexp,'split');
    me1{1} = me{2};
end
lbcopH = ['@(' ne1{1} ') [' ];

lbcstr = func2str(A.lbc);
%lbcstr = replace(lbcstr,ne1{1},['x,' ne1{1}]);
lbcstr = [ '@(x,' ne1{1} ')[' me1{1} ']'];
lbcop = chebop(A.domain);
lbcstr = [lbcstr ';'];
eval(['lbcop.op = ',lbcstr]);
lbcop = linearize(lbcop,0*lbcop.init);
[lbcr,lbcc] = size(lbcop);
blbc = zeros(lbcr,n*lbcc);
for ii = 1:lbcr
    for jj = 1:lbcc
        if ( isa(lbcop.blocks{ii,jj}, 'chebfun') )
            b = lbcop.blocks{ii,jj};
            lbcop.blocks{ii,jj} = operatorBlock.mult(b, b.domain);
        end
    end
end

for i = 1:lbcr
    for j = 1:lbcc
        tempMat = lbcop.blocks{i,j}.toCoeff;
        lentempMat = length(tempMat);
        for k = 1:n
            if (k <= lentempMat)
                blbc(i,(j-1)*n + k) = tempMat.blocks{lentempMat - k + 1}(A.domain(1));
            end
        end
    end
end


ne = split(func2str(A.rbc),'(');
ne1 = split(ne{2},')');
ne2r = split(ne1{1},',');

srbc = nargin(A.rbc);
if (srbc(1) > 1)
me = split(func2str(A.rbc),'[');
me1 = split(me{2},']');
else
    catexp = '@\(\w*\)';
    me = regexp(func2str(A.rbc),catexp,'split');
    me1{1} = me{2};
end
rbcopH = ['@(' ne1{1} ') [' ];

% me = split(func2str(A.rbc),'[');
% me1 = split(me{2},']');
% rbcopH = ['@(' ne1{1} ') [' ];

rbcstr = func2str(A.rbc);
%lbcstr = replace(lbcstr,ne1{1},['x,' ne1{1}]);
rbcstr = [ '@(x,' ne1{1} ')[' me1{1} ']'];;
rbcop = chebop(A.domain);
rbcstr = [rbcstr ';'];
eval(['rbcop.op = ',rbcstr]);
rbcop = linearize(rbcop,0*rbcop.init);
[rbcr,rbcc] = size(rbcop);
brbc = zeros(rbcr,n*rbcc);
for ii = 1:rbcr
    for jj = 1:rbcc
        if ( isa(rbcop.blocks{ii,jj}, 'chebfun') )
            b = rbcop.blocks{ii,jj};
            rbcop.blocks{ii,jj} = operatorBlock.mult(b, b.domain);
        end
    end
end
for i = 1:rbcr
    for j = 1:rbcc
        tempMat = rbcop.blocks{i,j}.toCoeff;
        lentempMat = length(tempMat);
        for k = 1:n
            if (k <= lentempMat)
                brbc(i,(j-1)*n + k) = tempMat.blocks{lentempMat - k + 1}(A.domain(2));
            end
        end
    end
end



if (rank(blbc)~= lbcr || rank(brbc)~= rbcr)
    error('Boundary conditions are not linearly independent');
end
solutionL = null(blbc);
solutionR = null(brbc);
blbcad = (Amin*solutionL)';
brbcad = (Aplus*solutionR)';
blbcad = blbcad/norm(blbcad);
brbcad = brbcad/norm(brbcad);

% Remove linearly dependent bcs:
[~,~,P] = qr(blbcad.');
blbcad = P*blbcad;
blbcad = blbcad(1:rank(blbcad.'),:);

[~,~,P] = qr(brbcad.');
brbcad = P*brbcad;
brbcad = brbcad(1:rank(brbcad.'),:);

bclim = max(orderAstarbc);

lbcadr = rank(blbcad);
for k = 1:lbcadr
    for l = 1:n
        for i = 1:length(ne2l)
            if (bclim(i) >= l-1)
                lbcopH =  [lbcopH ['(' num2str(blbcad(k,(i-1)*n + l)) ')'] '*diff(' ne2l{i} ',' int2str(l-1) ')+ '];
            end
        end
    end
    lbcopH(end-1:end) = [];
    lbcopH = [lbcopH '; ...' newline];
end
lbcopH(end-5:end) = [];
lbcopH = [lbcopH '];'];
lbcopH = eval(lbcopH);

rbcadr = rank(brbcad);
for k = 1:rbcadr
    for l = 1:n
        for i = 1:length(ne2l)
            if (bclim(i) >= l-1)
                rbcopH =  [rbcopH ['(' num2str(brbcad(k,(i-1)*n + l)) ')'] '*diff(' ne2l{i} ',' int2str(l-1) ')+ '];
            end
        end
    end
    rbcopH(end-1:end) = [];
    rbcopH = [rbcopH '; ...' newline];
end
rbcopH(end-5:end) = [];
rbcopH = [rbcopH '];'];
rbcopH = eval(rbcopH);

end
