function MaskPer = BraggOnly(SpectrumName,Wavelength,RL,resol,material,x,thermoOpticCoeff,temperature,count)
% This program calculates the first few modes at the Bragg wavelength to
% find out the mask period to be used in CCoeff
%%
disp('calculating the core mode effective index at the Bragg wavelength ')
global nL r_nf nf
global nLB 
M = 1;                                   % Only azimuthal order number used (only M=1 needed for fund core mode)
Neig = 1;						 	     % Number of eigenvalues demanded from the program


%======================================================================================


k0 = 2*pi/Wavelength;          				% [MR: free space wave vector, constant]
N_layer=numel(RL)-1;                        % Calculates number of layers from RL input
for i=1:N_layer
    nL(i) = Ref_indexV17(material{i},x(i),Wavelength, thermoOpticCoeff,temperature,count)
end
nLB=nL(1);
          

%============================= Generation of the index profile =============================

% Generation of the radial position and index of refraction vectors. Points
% of occurence of an index discontinuity are repeated in order to give both
% values of the refractive index on either side of the jump. This is used
% later to process index discontinuities. 

N_inc = ceil(diff(RL)./resol);              % vector with one row and the number of radial points in each layer in each column
numel(N_inc);                              % Number of elements in the cell here equal to dim(diff(RL)) = dim(RL) - 1
resolc = diff(RL)./N_inc;                   % Recalculate resolution in microns for each layer (vector with N_layer columns)
r_nf = cell(1,numel(N_inc));                % Initialization of cell for radial positions
nf = cell(1,numel(N_inc));                  % Initialization of cell for refractive index profile

for i = 1:numel(N_inc)
    r_nf{i} = RL(i) + resolc(i)*(0:N_inc(i)).';  %is a cell element with the radial locations for this element
    nf{i} = nL(i)*ones(N_inc(i) + 1,1);           % is a cell element with the refractive index values for this element
end


%================== Initializations ===================

N_layer = numel(r_nf);                      %Number of radial layers including the surrounding medium
nf2 = cell(1,N_layer);                      %Index of refraction squared
dnf2 = cell(1,N_layer);                     %Derivative of the index of refraction squared
Dr = cell(1,N_layer);                       %Radial component of the Electric Displacement field
Er = cell(1,N_layer);                       %Radial component of the Electric field
Ephi = cell(1,N_layer);                     %Azimutal componentof the electric field
Ez = cell(1,N_layer);                       %Longitudinal component of the electric field
dDr = cell(1,N_layer);                      %Derivative of the displacement field radial component
dEphi = cell(1,N_layer);                    %Derivative of the electric field azimutal component

%================== Processing of the index profile ===================

inc = zeros(1,numel(r_nf));

for i = 1:N_layer
    t_layer = r_nf{i}(end) - r_nf{i}(1);    % [MM : radial layer thickness]
    inc(i) = t_layer/(numel(r_nf{i}) - 1);  % [MM : radial increment on each layer]
end

% The index value at r = 0 is adjusted to ensure a zero index gradient at
% this point.

nf{1}(1) = (4*nf{1}(2) - nf{1}(3))/3; 

for i = 1:N_layer
    nf2{i} = nf{i}.^2;
    dnf2{i} = deriv_c(nf2{i},inc(i));
end

% Creation of vectors covering the whole fiber with doubled points at
% discontinuities

r_d = cat(1,r_nf{:});                       % [MM : position vector with doubled points at index jumps]
ni_d = cat(1,nf{:});                        % [MM : index of refraction with doubled points at index jumps]
ni2_d = cat(1,nf2{:});
dni2_d = cat(1,dnf2{:});

% Elimination of the point at r = 0 and of the end point assumed to be that
% at which the field vanishes (guidance condition).

r = r_d(2:(end - 1));                        
ni = ni_d(2:(end - 1));
ni2 = ni2_d(2:(end - 1));
dni2 = dni2_d(2:(end - 1));

% Refractive index discontinuities

jump_index = find(eq(diff(r),0));           % [MM : localization of doubled points (two values at same radius)]
r_jump = r(jump_index);                     % [MM : position of discontinuities] 
ni2_in = ni2(jump_index);                   % [MM : square of the index of refraction inside of the discontinuities]
ni2_out = ni2(jump_index + 1);              % [MM : square of the index of refraction outside of the discontinuities]

% Radial increment vector

h = inc(1)*ones(1,numel(ni2));
N = 1;
while lt(N,N_layer)
    h(gt(r,r_jump(N))) = inc(N + 1);
    N = N + 1;
end

% Elimination of discontinuities

r([jump_index (jump_index + 1)]) = [];
r2 = r.^2;
ni2([jump_index (jump_index + 1)]) = [];
dni2([jump_index (jump_index + 1)]) = [];
h([(jump_index - 1) jump_index]) = [];

%============= CALCULATION OF MODES ==========

%================== Construct d_x Operator ==================================

N = numel(r);
B = zeros(N,3); 						% [MR: Setup a Nx3 Zero Matrix. It is merely used as a step towards D]
B(:,1) = -(1/2)./h;						% [MR: Set the first column to -1's]
B(:,3) = (1/2)./circshift(h,[1,0]);			% [MR: Set the last column to 1's]
D_rho = spdiags(B,[-1,0,1],N,N); 		% [MR: FDM approx for first derivative. Useful in implementing L operator in Eq. (2.36), NxN matrix]
D_phi = D_rho;                          % [MR: FDM approx for first derivative. Useful in implementing L operator in Eq. (2.36), NxN matrix]
D_11 = D_rho(1,1);
D_12 = D_rho(1,2);

%================== Construct d_x^2 Operator ================================
% [MR: O(h^2) FDM approx for 2nd derivative]

h2 = h.^2;
C = zeros(N,3);								% [MR: Used as a step towards D2, Nx3 Ones Matrix]
C(:,1) = 1./h2;
C(:,2)= -2./h2;									% [MR: Set the second column to -2's]
C(:,3) = 1./circshift(h2,[1,0]);
D2_rho = spdiags(C,[-1,0,1],N,N);		% [MR: FDM approx for second derivative. Useful in implementing the L operator in Eq. (2.36), NxN matrix]
D2_phi = D2_rho;
D2_11 = D2_rho(1,1);
D2_12 = D2_rho(1,2);

%=====Contribution of discontinuities to D an D2================================

jump = jump_index + 1 - (1:numel(jump_index))'*2;
N_jump = numel(jump);
h_m = h(jump);              %[MM: Radial increment inside of the discontinutity]
h_p = h(jump + 1);          %[MM: Radial increment outside of the discontinuity]
h_mp = h_m'./ni2_out;
h_pp = h_p'./ni2_in;

for i = 1:N_jump

    fac = 3*r_jump(i)*(h_mp(i) + h_pp(i)) + 2*(ni2_out(i) - ni2_in(i))*h_mp(i)*h_pp(i);
    fac2 = 2*h_m(i)*fac; 
    D_rho(jump(i),jump(i) - 1) = -(fac + r_jump(i)*h_pp(i))/fac2 ;
    D_rho(jump(i),jump(i)) = 4*r_jump(i)*h_pp(i)/fac2;
    D_rho(jump(i),jump(i) + 1) = 4*r_jump(i)*h_mp(i)/fac2;
    D_rho(jump(i),jump(i) + 2) = -r_jump(i)*h_mp(i)/fac2;
    fac2 = 2*h_p(i)*fac;
    D_rho(jump(i) + 1,jump(i) - 1) = r_jump(i)*h_pp(i)/fac2 ;
    D_rho(jump(i) + 1,jump(i)) = -4*r_jump(i)*h_pp(i)/fac2;
    D_rho(jump(i) + 1,jump(i) + 1) = -4*r_jump(i)*h_mp(i)/fac2;
    D_rho(jump(i) + 1,jump(i) + 2) = (r_jump(i)*h_mp(i) + fac)/fac2;
    fac2 = h_m(i)^2*fac; 
    D2_rho(jump(i),jump(i) - 1) = (fac - r_jump(i)*h_pp(i))/fac2;
    D2_rho(jump(i),jump(i)) = -2*(fac - 2*r_jump(i)*h_pp(i))/fac2;
    D2_rho(jump(i),jump(i) + 1) = 4*r_jump(i)*h_mp(i)/fac2;
    D2_rho(jump(i),jump(i) + 2) = -r_jump(i)*h_mp(i)/fac2;
    fac2 = h_p(i)^2*fac;
    D2_rho(jump(i) + 1,jump(i) - 1) = -r_jump(i)*h_pp(i)/fac2 ;
    D2_rho(jump(i) + 1,jump(i)) = 4*r_jump(i)*h_pp(i)/fac2;
    D2_rho(jump(i) + 1,jump(i) + 1) = -2*(fac - 2*r_jump(i)*h_mp(i))/fac2;
    D2_rho(jump(i) + 1,jump(i) + 2) = (fac - r_jump(i)*h_mp(i))/fac2;
    fac = 6*h_m(i)*(h_m(i) + h_p(i)); 
    D_phi(jump(i),jump(i) - 1) = -(3*h_m(i) + 4*h_p(i))/fac;
    D_phi(jump(i),jump(i)) = 4*h_p(i)/fac;
    D_phi(jump(i),jump(i) + 1) = 4*h_m(i)/fac;
    D_phi(jump(i),jump(i) + 2) = -h_m(i)/fac;
    fac = 6*h_p(i)*(h_m(i) + h_p(i));
    D_phi(jump(i) + 1,jump(i) - 1) = h_p(i)/fac ;
    D_phi(jump(i) + 1,jump(i)) = -4*h_p(i)/fac;
    D_phi(jump(i) + 1,jump(i) + 1) = -4*h_m(i)/fac;
    D_phi(jump(i) + 1,jump(i) + 2) = (4*h_m(i) + 3*h_p(i))/fac;
    fac = 3*h_m(i)^2*(h_m(i) + h_p(i)); 
    D2_phi(jump(i),jump(i) - 1) = (3*h_m(i) + 2*h_p(i))/fac;
    D2_phi(jump(i),jump(i)) = -(6*h_m(i) + 2*h_p(i))/fac;
    D2_phi(jump(i),jump(i) + 1) = 4*h_m(i)/fac;
    D2_phi(jump(i),jump(i) + 2) = -h_m(i)/fac;
    fac = 3*h_p(i)^2*(h_m(i) + h_p(i));
    D2_phi(jump(i) + 1,jump(i) - 1) = -h_p(i)/fac ;
    D2_phi(jump(i) + 1,jump(i)) = 4*h_p(i)/fac;
    D2_phi(jump(i) + 1,jump(i) + 1) = -(2*h_m(i) + 6*h_p(i))/fac;
    D2_phi(jump(i) + 1,jump(i) + 2) = (2*h_m(i) + 3*h_p(i))/fac;
    
end

%================== Useful Intermediate Operators ================================

p = (1./r) - dni2./ni2;                 % [MR: N component column vector, useful for Eq. (2.21) for E_phi (and E_rho)]
P = spdiags(p,0,N,N);                   % [MR: NxN Operator, useful for Eq. (2.21)]
q = ni2*k0^2 - dni2./ni2./r;            % [MM]
a_ = 2*ni2./r2;                         % [MM: Upper off-diagonal of (2.56) N component column vector]
a = spdiags(a_,0,N,N);					% [MR: NxN diagonal matrix (operator) with u1 along the diagonal]
b_ = (2./r2 + dni2./ni2./r)./ni2;       % [MM]
b = spdiags(b_,0,N,N);                  % [MM] 
Ri = spdiags((1./r),0,N,N);             % [MR: NxN Matrix/Operator form of 1/rho -- to be used in L-hat (2.36)]

%================== Contribution of discontinuities to b ================================ 

for i = 1:N_jump
    
    fac = 3*r_jump(i)*(h_mp(i) + h_pp(i)) + 2*(ni2_out(i) - ni2_in(i))*h_mp(i)*h_pp(i);
    fac2 = h_mp(i)*h_pp(i)*(ni2_out(i) - ni2_in(i))/(3*h_m(i)*(h_m(i) + h_p(i))*fac);
    fac3 = (2/h_m(i) + 1/r(jump(i)))*fac2; 
    b(jump(i),jump(i) - 1) = -h_pp(i)*fac3 ;
    b(jump(i),jump(i)) = 4*h_pp(i)*fac3;
    b(jump(i),jump(i) + 1) = 4*h_mp(i)*fac3;
    b(jump(i),jump(i) + 2) = -h_mp(i)*fac3;
    fac2 = h_mp(i)*h_pp(i)*(ni2_out(i) - ni2_in(i))/(3*h_p(i)*(h_m(i) + h_p(i))*fac);
    fac3 = (2/h_p(i) - 1/r(jump(i) + 1))*fac2; 
    b(jump(i) + 1,jump(i) - 1) = -h_pp(i)*fac3;
    b(jump(i) + 1,jump(i)) = 4*h_pp(i)*fac3;
    b(jump(i) + 1,jump(i) + 1) = 4*h_mp(i)*fac3;
    b(jump(i) + 1,jump(i) + 2) = -h_mp(i)*fac3;

end



   
    % Corrections aux opérateurs dérivée lorsque M = 1
    

        D_rho(1,1) = -2/3/h(1);             % Corrections to the d_x operator at r=h that takes into account boundary condition at r=0
        D_rho(1,2) = 2/3/h(1);
        D2_rho(1,1) = -2/3/h(1)^2;          % Similar BC correction to d_x^2 at origin
        D2_rho(1,2) = 2/3/h(1)^2;
        D_phi(1,1) = -2/3/h(1);             % Corrections to the d_x operator at r=h that takes into account boundary condition at r=0
        D_phi(1,2) = 2/3/h(1);
        D2_phi(1,1) = -2/3/h(1)^2;          % Similar BC correction to d_x^2 at origin
        D2_phi(1,2) = 2/3/h(1)^2;
    
    %================== Useful Intermediate Operators ================================
    
    v_tmp = (M^2 + 1)./r2;
    u_rho = q - v_tmp;                      % [MM: N-component column vector, Eq. (2.35) --  definition of potential well]
    U_rho = spdiags(u_rho,0,N,N); 			% [MR: NxN Matrix/Operator Form of Eq. (2.35)]
    u_phi = ni2*k0^2 - v_tmp;               % [MM: N-component column vector, Eq. (2.35) --  definition of potential well]
    U_phi = spdiags(u_phi,0,N,N); 			% [MR: NxN Matrix/Operator Form of Eq. (2.35)]

    %================== Construct L2hat Operator - Eq. (2.57)=======================

    L2 = D2_rho + P*D_rho + U_rho;                       	% [JA: Direct from (2.55) and (2.58)]

    %================== Construct Lhat Operator - Eq. (2.57) =======================

    L = D2_phi + Ri*D_phi + U_phi;                              	% [JA: Direct from (2.55)]

    %================ Final assembly ================%
    
    Ls = [L2, M*a; M*b, L];                     % [MR: Exact solution operator, a 2Nx2N matrix implementing Eq (2.56)]

    %====================== Solve for Eigenvalues =============================%
    
    Umax = real((nL(1)+nL(2))/2)^2*k0^2;        % [will find the eigenvalue  closest to the average between core and clad index
    [V,Eig] = eigs(Ls, Neig, Umax);		% [MR: This line finds Neig eigenvectors of Ls that are close to Umax, stored in V (2NxNeig matrix), and diagonal
                                        % NeigxNeig matrix of eigenvalues, Eig. This is solving Eq. (2.59)]
    Eig = diag(Eig);                    % [MR: creates a Neig-component vector from L's  diagonal elements:  these are the eigenvalues]
    gamma = sqrt(Eig);
    gamma = sign(real(gamma)).*gamma;               % [MM : on assure que la partie réelle de la constante de propagation est positive]
         
    %=============== Calculate Effective indices and orders ============
    
    Neff = gamma./k0;
%     [v,ind] = sort(real(Neff),'descend'); %Sometimes the modes (eigenvalues) do not come out in decreasing order  (v is a dummy variable)
%     Neff = Neff(ind)                     %Ordering them (as done here) seems to help getting better spectra
%     CoreIndex=find((nf{1,1}(1) > real(Neff)) & (real(Neff) > nf{1,2}(1)),1);    %Finds the index location of the first effective index that are between n_core and n_clad
    NeffcoreB = real(Neff)
    MaskPer = Wavelength/NeffcoreB;
       
    
    
end

function deriv_ = deriv_c(matrix,inc)

% Derivative along columns of a matrix

deriv_ = zeros(size(matrix));
deriv_(1,:) = 4*matrix(2,:) - 3*matrix(1,:) - matrix(3,:);
deriv_(end,:) = matrix(end - 2,:) + 3*matrix(end,:) - 4*matrix(end - 1,:);
deriv_(2:(end - 1),:) = matrix(3:end,:) - matrix(1:(end - 2),:);
deriv_ = deriv_/2/inc;
end

