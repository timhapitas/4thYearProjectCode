function [dummy, Neff] = Modes(SpectrumName,Wavelength,Mmax,Neig,Neffmax)
% This program calculates "Neig" modes for each M=0 to "Mmax" around "Neffmax" at
% "Wavelength"
% No need to recalculate index profile here since nf(i) was updated in the
% "Coreonly" function
global nL r_nf nf
k0 = 2*pi/Wavelength;          				% k0 at calc wavelength
 
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

%Note: for M=0, the first mode is TE01 and the second one TM01, and they
%continue to alternate (low order TE and TM of the same order have the same Neff but not high order ones)
%%
for M = 0 : Mmax;
tic;  %start the timer for each mode family

    fprintf('Calculating mode family # %d\n', M);
    
    % Corrections aux opérateurs dérivée lorsque M = 1
    
    if eq(M,1)
        D_rho(1,1) = -2/3/h(1);             % Corrections to the d_x operator at r=h that takes into account boundary condition at r=0
        D_rho(1,2) = 2/3/h(1);
        D2_rho(1,1) = -2/3/h(1)^2;          % Similar BC correction to d_x^2 at origin
        D2_rho(1,2) = 2/3/h(1)^2;
        D_phi(1,1) = -2/3/h(1);             % Corrections to the d_x operator at r=h that takes into account boundary condition at r=0
        D_phi(1,2) = 2/3/h(1);
        D2_phi(1,1) = -2/3/h(1)^2;          % Similar BC correction to d_x^2 at origin
        D2_phi(1,2) = 2/3/h(1)^2;
    end

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
    
    Umax = Neffmax^2*k0^2;        % [MR: Constant. This sets an upper 'bound' for B_max^2 (effective wave vector squared) equal to n1^2*k^2]
    [V,Eig] = eigs(Ls, Neig, Umax);		% [MR: This line finds Neig eigenvectors of Ls that are close to Umax, stored in V (2NxNeig matrix), and diagonal
                                        % NeigxNeig matrix of eigenvalues, Eig. This is solving Eq. (2.59)]
    Eig = diag(Eig);                    % [MR: creates a Neig-component vector from L's  diagonal elements:  these are the eigenvalues]
    
    %================ Fields at r = 0 and r = rho_max ====================%
    
    if eq(M,1)
        Dr_0 = (4*V(1,:) - V(2,:))/3;
        Ephi_0 = (4*V(1 + end/2,:) - V(2 + end/2,:))/3;
    else
        Dr_0 = zeros(1,Neig);
        Ephi_0 = zeros(1,Neig);
    end

    Dr_rho_max = zeros(1,Neig);
    Ephi_rho_max = zeros(1,Neig);
        
    %================ Fields at index discontinuities ====================%
    
    Dr_jump = zeros(numel(jump),Neig);
    Ephi_jump = zeros(numel(jump),Neig);

    for i = 1:N_jump
        fac = 3*r_jump(i)*(h_mp(i) + h_pp(i)) + 2*(ni2_out(i) - ni2_in(i))*h_mp(i)*h_pp(i); 
        Dr_jump(i,:) = r_jump(i)/fac*(-h_pp(i)*V(jump(i) - 1,:) + 4*h_pp(i)*V(jump(i),:) + 4*h_mp(i)*V(jump(i) + 1,:) -h_mp(i)*V(jump(i) + 2,:));
        Ephi_jump(i,:) = (-h_p(i)*V(jump(i) - 1 + end/2,:) + 4*h_p(i)*V(jump(i) + end/2,:) + 4*h_m(i)*V(jump(i) + 1 + end/2,:) -h_m(i)*V(jump(i) + 2 + end/2,:))/3/(h_p(i) + h_m(i));
        Ephi_jump(i,:) = Ephi_jump(i,:) + 2*M*h_mp(i)*h_pp(i)*(ni2_out(i) - ni2_in(i))*Dr_jump(i,:)/r_jump(i)/3/(h_p(i) + h_m(i)); 
    end

    %================ Complete Radial and Azimutal Fields ====================%
    
    Dr{1} = [Dr_0;V(1:jump(1),:);Dr_jump(1,:)];
    Ephi{1} = [Ephi_0;V((1 + end/2):(jump(1) + end/2),:);Ephi_jump(1,:)];
    for i = 1:(N_layer - 2)
        Dr{i + 1} = [Dr_jump(i,:);V((jump(i) + 1):jump(i + 1),:);Dr_jump(i + 1,:)];
        Ephi{i + 1} = [Ephi_jump(i,:);V((jump(i) + 1 + end/2):(jump(i + 1) + end/2),:);Ephi_jump(i + 1,:)];
    end
    Dr{N_layer} = [Dr_jump(end,:);V((jump(end) + 1):(end/2),:);Dr_rho_max];
    Ephi{N_layer} = [Ephi_jump(end,:);V((jump(end) + 1 + end/2):end,:);Ephi_rho_max];

    %================ Field derivatives ====================%
    
    for i = 1:N_layer
        dDr{i} = deriv_c(Dr{i},inc(i));
        dEphi{i} = deriv_c(Ephi{i},inc(i));
    end
    
    %================ Ez component ====================%

    gamma = sqrt(Eig);
    gamma = sign(real(gamma)).*gamma;               % [MM : on assure que la partie réelle de la constante de propagation est positive]
    
    for i = 1:N_layer
        Ez{i} = -j./(ones(numel(r_nf{i}),1)*gamma.').*(dDr{i}./(nf2{i}*ones(1,Neig)) + Dr{i}./((r_nf{i}.*nf2{i})*ones(1,Neig)) - M*Ephi{i}./(r_nf{i}*ones(1,Neig)));
    end

    if eq(M,0)
        Ez{1}(1,:) = (4*Ez{1}(2,:) - Ez{1}(3,:))/3;
    else
        Ez{1}(1,:) = 0;
    end

    %====================== Normalize Eigenvectors =============================%
    
    if eq(M,0)
        
        norme1 = zeros(1,Neig);
        norme2 = zeros(1,Neig);
        
        for i = 1:N_layer
            norme1 = norme1 + Integral((r_nf{i}./nf2{i}).',Dr{i}.^2,inc(i));
            norme2 = norme2 + Integral(r_nf{i}.',Ephi{i}.^2,inc(i));
        end

        norme = 2*pi*(k0*norme1./gamma.' + gamma.'/k0.*norme2);
    
    else
            
        norme1 = zeros(1,Neig);
        norme2 = zeros(1,Neig);
        norme3 = zeros(1,Neig);
        norme4 = zeros(1,Neig);
        
        for i = 1:N_layer
            norme1 = norme1 + Integral(r_nf{i}.',Ephi{i}.^2,inc(i));
            norme2 = norme2 + Integral(M^2./r_nf{i}.',Ephi{i}.^2,inc(i));
            norme3 = norme3 + Integral((k0^2*r_nf{i}./nf2{i} - M^2./r_nf{i}./nf2{i}.^2).',Dr{i}.^2,inc(i));
            norme4 = norme4 + Integral(M./nf2{i}.',Dr{i}.*dEphi{i} - Ephi{i}.*dDr{i},inc(i));
        end
        norme = pi/k0./gamma.'.*(Eig.'.*norme1 + norme2 + norme3 + norme4);
    end

    for i = 1:N_layer
        mat_norme = sqrt(ones(numel(r_nf{i}),1)*norme);
        Dr{i} = Dr{i}./mat_norme;
        Er{i} = Dr{i}./(nL(i)^2);      % these are the final field with N_layer cell elements each containing a matrix with Neig columns
                                       % and a number of lines equal to the number of grid points       
        Ephi{i} = Ephi{i}./mat_norme;
        Ez{i} = Ez{i}./mat_norme;      % but they are not sorted yet (this is done in the next step
    end

        
    %=============== Calculate Effective indices and orders ============
    
    Neff = gamma./k0;
    [v,ind] = sort(real(Neff),'descend'); %Sometimes the modes (eigenvalues) do not come out in decreasing order  (v is a dummy variable)
    Neff = Neff(ind)                     %Ordering them (as done here) seems to help getting better spectra
    
    for i = 1:N_layer
        Er{i} = Er{i}(:,ind);           %this orders the fields to match the reordered Neff
        Ephi{i} = Ephi{i}(:,ind);
        Ez{i} = Ez{i}(:,ind);
    end

    if eq(M,0)
        GuidedIndex=((nf{1,1}(1) > real(Neff)) & (real(Neff) > nf{1,N_layer}(1))); %non leaky modes only
        Er{1} = Er{1}(:,GuidedIndex);           %this orders the fields to match the reordered Neff
        Ephi{1} = Ephi{1}(:,GuidedIndex);       %this only ordered the fields in the core
        Ez{1} = Ez{1}(:,GuidedIndex);           %therefore the fields in the other layers are no longer "indexed" the same way and must not be used
        Neff = Neff(GuidedIndex);
        Neig = numel(Neff);                    %this will ensure that following mode sets will have the same number of modes  
    else
        GuidedIndex= find(real(Neff),Neig);
        Er{1} = Er{1}(:,GuidedIndex);           %this orders the fields to match the reordered Neff
        Ephi{1} = Ephi{1}(:,GuidedIndex);       %this only ordered the fields in the core
        Ez{1} = Ez{1}(:,GuidedIndex);           %therefore the fields in the other layers are no longer "indexed" the same way and must not be used
        Neff = Neff(GuidedIndex);
    end
    %=============== Save results ============
    %NOTE: when the index profile has absorbing layers, the fields are
    %complex
    %NOTE: It would be possible to store the TE and TM modes separately
    %because they come out in alternation (or we could check for the Ez
    %component. But then, earlier in the program we would have to double
    %Neig when M=0
    
    save([SpectrumName,'/Neff_', int2str(M)], 'Neff');	% [MR: Save the effective indices for each propagation constants (as MATLAb variables
                                            %      for the given m and documents the m value. The file to be saved is 'Neff_M']
    save([SpectrumName,'/Er_', int2str(M)], 'Er');       % [MR: Save the normalized eigenvectors into file named Er_M, and the variable name of the cell array is Er
    save([SpectrumName,'/Ephi_', int2str(M)], 'Ephi');       % [MR: Save the normalized eigenvectors into file named Ephi_M, and the variable name of the cell array is Ephi
    save([SpectrumName,'/Ez_', int2str(M)], 'Ez');       % [MR: Save the normalized eigenvectors into file named Ez_M, and the variable name of the cell array is Ez
    
    dummy = 'All modes saved';
    
    % Élimination des corrections aux opérateurs dérivée après calcul pour M = 1
    
    if eq(M,1)
        D_rho(1,1) = D_11;                        % Corrections to the d_x operator at r=h that takes into account boundary condition at r=0
        D_rho(1,2) = D_12;
        D2_rho(1,1) = D2_11;                     % Similar BC correction to d_x^2 at origin
        D2_rho(1,2) = D2_12;
        D_phi(1,1) = D_11;                        % Corrections to the d_x operator at r=h that takes into account boundary condition at r=0
        D_phi(1,2) = D_12;
        D2_phi(1,1) = D2_11;                     % Similar BC correction to d_x^2 at origin
        D2_phi(1,2) = D2_12;
    end
%====================================================
%       FACULATIVE BLOCK TO PLOT FIELDS
%====================================================
%figure;hold on;grid on
%     for i = 1:N_layer
%         plot(r_nf{i},Dr{i});
%     end
%     xlabel('r');ylabel('Dr');title(['M = ',num2str(M)]);
%     
%     figure;hold on;grid on
%     for i = 1:N_layer
%         plot(r_nf{i},Ephi{i});
%     end
%     xlabel('r');ylabel('E\phi');title(['M = ',num2str(M)]);
%     
%     figure;hold on;grid on
%     for i = 1:N_layer
%         plot(r_nf{i},imag(Ez{i}));
%     end
%     xlabel('r');ylabel('Ez');title(['M = ',num2str(M)]);
toc
end

function deriv_ = deriv_c(matrix,inc)

% Derivative along columns of a matrix

deriv_ = zeros(size(matrix));
deriv_(1,:) = 4*matrix(2,:) - 3*matrix(1,:) - matrix(3,:);
deriv_(end,:) = matrix(end - 2,:) + 3*matrix(end,:) - 4*matrix(end - 1,:);
deriv_(2:(end - 1),:) = matrix(3:end,:) - matrix(1:(end - 2),:);
deriv_ = deriv_/2/inc;
return

function integral_ = Integral(vector,matrix,dr)

% This program performs multiple trapezoidal integrations through a matrix
% multiplication. When filter = 0, the first element of the integrand is
% set to zero (to avoid indeterminations such as 0/0).
% vector : line vector with N elements (1 x N)
% matrix : matrix with N lines and M columns (N x M)
% dr : scalar
% integral_ : line vector with M elements

if or(isnan(vector(1)),isinf(vector(1)))
    vector(1) = 0;
end
integral_ = vector*matrix;
integral_ = integral_ - 0.5*(vector(1)*matrix(1,:) + vector(end)*matrix(end,:));
integral_ = dr*integral_;
return
