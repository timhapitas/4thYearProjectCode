function dummy = CCoeff(SpectrumName,MaskPer,theta_int,Mmax)

global nL   							    % [MR: Refractive index profile, 3-component row vector] and radii of fiber layer boundaries in microns

global r_nf                                 % Cell array  of the radial grid points (starting at h, not at zero)
global nf                                   % Cell array of the refractive index profile at each radial point

global Mpm CpmP CpmS Npm 

%=========================================
%Loading input core mode fields in the core only
load([SpectrumName,'/E1r']);
load([SpectrumName,'/E1p']);
load([SpectrumName,'/E1z']);

%Loading and calculating for M=1 family first............
M=1;

load([SpectrumName,'/Er_1']);             					%  Load radial E component of all modes stored in matrix Er_1.mat, under the variable name "Er"] 
load([SpectrumName,'/Ephi_1']);                             %  Load azimuthal E component of all modes stored in matrix Ephi_1.mat, under the variable name "Ephi"]
load([SpectrumName,'/Ez_1']);             					%  Load axial E component of all modes stored in matrix Ez_1.mat, under the variable name "Ez"] 
load([SpectrumName,'/Neff_1']);                 %  Load Effective index of all modes stored in matrix Neff_1.mat, under the variable name Neff


Vnsize = size(Er{1});                           % finds the size of the matrix in the first cell of Er 
Nmax = Vnsize(2);                               % first dim is number of radial points in cell, second number is number of modes found
                                                % for this value of M
                                                                                                
                                                %Initialize matrices
Mpm=zeros(Mmax+1,Nmax);                         % M values of the modes to be coupled to by the incident core mode
Npm=zeros(Mmax+1,Nmax);                         % Effective indices of the modes
CpmP=zeros(Mmax+1,Nmax);                        % Coupling coefficients (C+- factors from Michel)for even-even (HE11x input, of 'P') coupling
CpmS=zeros(Mmax+1,Nmax);                        % Coupling coefficients (C+- factors from Michel)for odd-odd (HE11y input, or 'S')
                                                % Note that for the chosen
                                                % axis orientation of
                                                % Michel these are the only
                                                % allowed couplings
                                                
%=========== Restrict variables to core region only for integration of C+- factors =========%

r_c = r_nf{1};										% [MR: coreIndex x 1 column vector. Restrict spatial extent to core radius]
nc = nL(1);                                        % the core index for use in the integrals

%  Variable naming convention: 
%		E = Electric field  ;  D = Electric displacement (= epsilon * E)

%		r/p/z = spatial component
%		c = in the core]

E1_rc = Er{1};		% The first cell of Er is a matrix with rows corresponding to radial positions in the core and columns representing modes of the M=1 family
E1_pc = Ephi{1};	
E1_zc = Ez{1};		


% Calculate the coupling integrals for the M=1 modes, including the self coupling between the "HE11" mode and itself
%(i.e. the first mode of the M=1 group

for j=1 : Nmax        %Loop through all the M=1 modes
    
    E2r=E1_rc(:,j);
    E2p=E1_pc(:,j);
    E2z=E1_zc(:,j);        

    Ir0 = f_Iterm(nc,r_c,MaskPer,theta_int,E1r,E2r,0);   %the overlap intergrals "I" are defined in function Iterm
    Ip0 = f_Iterm(nc,r_c,MaskPer,theta_int,E1p,E2p,0);
    Iz0 = f_Iterm(nc,r_c,MaskPer,theta_int,E1z,E2z,0);
    
    Ir2 = f_Iterm(nc,r_c,MaskPer,theta_int,E1r,E2r,2);
    Ip2 = f_Iterm(nc,r_c,MaskPer,theta_int,E1p,E2p,2);
    Iz2 = f_Iterm(nc,r_c,MaskPer,theta_int,E1z,E2z,2);
    
    %row index for Cpm is M+1 in order for the M=0 mode family will correspond to first line of the coefficient matrix.
     
    CpmP(M+1,j) = Ir0 + Ip0 - Ir2 + Ip2 + Iz0 - Iz2;    %  Formula for even-even coupling (P-polarized input) (i.e. HE11x input when tilt is in x-z)
     
    CpmS(M+1,j) = Ir0 + Ip0 + Ir2 - Ip2 + Iz0 + Iz2;    %  Formula for odd-odd coupling (S-polarized input) (i.e. HE11y input when tilt is in x-z)
    
                   
end
Npm(M+1,:) = Neff(:);                      % building the corresponding matrices of effective indices and M values
Mpm(M+1,:) = M;

%Now do the group of modes with M=0, for which only one I integral is
%needed. 

M=0;     

load([SpectrumName,'/Er_0']);             					% variable name is again Er but that is OK since we saved the incident mode in another variable 
load([SpectrumName,'/Ephi_0']);             				
load([SpectrumName,'/Ez_0']);             					
load([SpectrumName,'/Neff_0']);
                                                                                                             
%There is a small problem with M=0 because the solutions returned by the mode solver alternate between TE and TM (with TM being the even mode and TE the odd mode)
%BUT this problems solves itself since
%it turns out that the Cpm values alternate between near zero and finite
%values, i.e. CpmS has zeros where CpmP has finite values and the opposite
                                                            
for j=1 : Nmax
    E2r = Er{1}(:,j);			% Defining new temporary field vectors for the fields in the core
    E2p = Ephi{1}(:,j);			
    E2z = Ez{1}(:,j);	
    
    Ir1 = f_Iterm(nc,r_c,MaskPer,theta_int,E1r,E2r,1);
    Ip1 = f_Iterm(nc,r_c,MaskPer,theta_int,E1p,E2p,1);
    Iz1 = f_Iterm(nc,r_c,MaskPer,theta_int,E1z,E2z,1);
    
     CpmP(M+1,j) = 2*(Ir1+Iz1);           %even-even (P) HE11x this is coupling to TM0
     CpmS(M+1,j) = 2*(Ip1);               %odd-odd (S) HE11y  this is coupling to TE0
  end
Npm(M+1,:) = Neff(:);
Mpm(M+1,:) = M;
%Finally, we can loop through the other M values in the same manner

for M = 2:Mmax;

load([SpectrumName,'/Er_',int2str(M)]);             					
load([SpectrumName,'/Ephi_',int2str(M)]);             				 
load([SpectrumName,'/Ez_',int2str(M)]);             					 
load([SpectrumName,'/Neff_',int2str(M)]);
    for j=1 : Nmax
        E2r = Er{1}(:,j);									% [MR: coreExtent by Nmax matrix of radial electric field component in the core]
        E2p = Ephi{1}(:,j);				% [MR: coreExtent by Nmax matrix of azimuthal electric field component in the core]
        E2z = Ez{1}(:,j);	% [MR: coreExtent by Nmax matrix of longitudinal electric field component in the core]
        kA=M-1;
        IrA = f_Iterm(nc,r_c,MaskPer,theta_int,E1r,E2r,kA);
        IpA = f_Iterm(nc,r_c,MaskPer,theta_int,E1p,E2p,kA);
        IzA = f_Iterm(nc,r_c,MaskPer,theta_int,E1z,E2z,kA);
        kB=M+1;
        IrB = f_Iterm(nc,r_c,MaskPer,theta_int,E1r,E2r,kB);
        IpB = f_Iterm(nc,r_c,MaskPer,theta_int,E1p,E2p,kB);
        IzB = f_Iterm(nc,r_c,MaskPer,theta_int,E1z,E2z,kB);
        
        if rem(M-1,2) == 0      %case when M-1 is even (formula 119), where 1 is the M number of the incident mode
            
        CpmP(M+1,j) = (-1)^(kA/2)*(IrA+IpA+IzA)+(-1)^(kB/2)*(IrB-IpB+IzB); %even-even (P) HE11x 
        CpmS(M+1,j) = (-1)^(kA/2)*(IrA+IpA+IzA)-(-1)^(kB/2)*(IrB-IpB+IzB); %odd-odd (S) HE11y        
        else                    %case when M-1 is odd
        CpmP(M+1,j) = (-1)^((kA-1)/2)*(IrA+IpA+IzA)+(-1)^((kB-1)/2)*(IrB-IpB+IzB);  %even-even (P) HE11x 
        CpmS(M+1,j) = (-1)^((kA-1)/2)*(IrA+IpA+IzA)-(-1)^((kB-1)/2)*(IrB-IpB+IzB);  %odd-odd (S) HE11y  
        end
    
    end
    Npm(M+1,:) = Neff(:);     % these two output matrices are the same for S and P polarized input.
    Mpm(M+1,:) = M;           % They are used to transfer the matrices of effective indices and M values to the next program
end

%This block saves the raw matrices of M, Neff, Cp and Cs just for info
CoeffP=CpmP;
CoeffS=CpmS;
Nefflist=Npm;
Mlist=Mpm;
save([SpectrumName,'/CoeffP'],'CoeffP');
save([SpectrumName,'/CoeffS'],'CoeffS');
save([SpectrumName,'/Nefflist'],'Nefflist');
save([SpectrumName,'/Mlist'],'Mlist');


%This block only transfers calculation matrices for the fiber modes (not
%the SPP one on the outer metal boundary
%and it transforms the 2D matrices into 1D arrays to facilitate their use
%in the following programs

FibermodeIndex = nf{1,1}(1) > real(Npm);  %Eliminates modes for which the Neff is larger than the core index nf{1,1}(1)
Npm=Npm(FibermodeIndex);
Mpm=Mpm(FibermodeIndex);
CpmS=CpmS(FibermodeIndex);
CpmP=CpmP(FibermodeIndex);

dummy = 'Coefficients complete';
end