
%{

This function calculates the right hand side of ODE system (128)
it is used as a function handle in ODE45 (in Spectrum1_Ja.m)
For the special case of M=1 and P=1 coupling (same azim family)
%}

function DE = f_MEcalc_P(z,E)



global wavLoc                               % Transver of variables from Step3.......

global MaskPer

global DnAC
global Neff_core modeMpm modeCpm modeNpm NewNmax


Neff1=Neff_core;
Neffcl=modeNpm;

kvec=2*pi/wavLoc;
locKmask = 4 * pi /MaskPer;                %is the local mask period (if we could have chirp, see below)

%in the chirp calculation we correct for the fact that z is zero at the beginning.
%locKmask = 2 * pi * (2-Chirp*(z-Grat_length/2)/MaskPer)/MaskPer;    %is the local (at this value of z when chirp is present) value of the mask period
                                                    %NOTE: the formula (128) assumes z=0 at the center of the grating (mask), here
                                                      
exparg = (kvec*(Neff1+Neffcl)-locKmask);     %will need to be multiplied by z in the exponential
                                                                             
                              
%NOW READY TO CALCULATE Matrix elements for ODE
        
Mat=zeros(NewNmax,NewNmax);
        
   for k_index = 2 : NewNmax    
        
        if ne(rem(modeMpm(k_index-1)-1,2),0)           % Checking for the parity of abs(M-P), where P=1 the azimuthal index of the core mode. This group is ODD
            Mat(1,k_index) = kvec*DnAC*modeCpm(k_index-1)*exp(1i*exparg(k_index-1)*z)/2; %First line of matrix for M-P odd  (same as even, except divided by -i)
            Mat(k_index,1) = kvec*DnAC*modeCpm(k_index-1)*exp(-1i*exparg(k_index-1)*z)/2; %First column of matrix for M-P odd (same as even except divided by i)
        else
            Mat(1,k_index) = -1i*kvec*DnAC*modeCpm(k_index-1)*exp(1i*exparg(k_index-1)*z)/2; %First line of matrix for M-P even
            Mat(k_index,1) = 1i*kvec*DnAC*modeCpm(k_index-1)*exp(-1i*exparg(k_index-1)*z)/2; % First column of matrix for M-P even
        end
   end
                                               
 DE = Mat*E;



