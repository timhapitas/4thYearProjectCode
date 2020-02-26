%Program TFBGsimVM
% In this version the refractive index profile is updated between the
% calculation of the core mode at the Bragg wavelength (to get the mask period) and
% the calculation of the core mode at the calculation wavelength. It is not
% recalculated afterwards when calculating cladding mode indices so material dispersion is neglected for the range
% of wavelengths of the spectrum
% For the core mode calculations, only 1 eigenmode is calculated, nearest
% to the average between the core and cladding refractive index
% this avoids problems when guidance is possible in one of the outer layers
% but it may cause problems if the core is not single mode (to be verified)
%  A few parameters are preset for the case of a MIP layer on top of SMF

global nLB nLC

currentOutsideIndex = 1.315;

troughs = cell(1, 1);
wavelengths = cell(1, 1);

thermExpCoeff = 5.5e-7; %K^-1
thermoOpticCoeff = 12.9e-6; %K^-1

temperature = [0 5 10 15 20 25 30]; %degree increase from room temperature (Mask Period obtained at room temperature from BraggOnly subroutine)

for count = 1:length(temperature)

%First call the solver at the Bragg wavelength to calculate the mask period
%of the experimental grating needed

mkdir SMFlowRef1550b     % a name for the output spectrum files (change it to save several spectra one after the other)
%NEED TO CHANGE IT BELOW AS WELL%%%%%%%%%%%%
SpectrumName = 'SMFlowRef1550b'  
 
BraggWav = 1.610       % the Bragg wavelength of grating (will be used to calculate the effective mask period

theta_int = 10        % the internal angle of the tilt
Length = 10000          % length of grating in um
Dn_AC = 0.004           % index modulation amplitude of the grating

Mmax = 1;               % the maximum M value (modes will be calculated for M = 0 to Mmax

%The next input parameters can be changed and then only the "Spectrum"
%function called to recalculate the spectrum with the same modes involved
InputPol='P'          % P=HE11x ; S=HE11y  ;  U=unpolarized input, i.e. 1/2*(P^2+S^2)

WavStart = 1.535       %Starting wavelength (in microns)
WavEnd = 1.620 %End wavelength

Npoints=1+round((WavEnd-WavStart)/0.0000001)   %For spectrum points every 20 pm
CalcWav = (WavStart+WavEnd)/2;                   % The median value of wavelength for the spectrum studied 

maxdetuning = .003   %For TFBGs with metal layers 0.005 is quite good (0.01 better but very slow)
                     %Otherwise, .001 is sufficient and much faster

Neig = round(1000*(WavEnd-WavStart)*6)              % number of modes found below NeffMax for each value of M (the azimuthal mode order) NEED TO COVER THE WHOLE SPECTRUM CALCULATED                    
                                                    %this particular
                                                    %estimate is based on
                                                    %about 2 modes per nm
                     
%FWG definition

rho_max = 80;                               % [MM : radius at which the field is assumed to vanish to ensure guidance of the mode]
RL = [0 4.1 62.5 rho_max];         % [MM: RL = [0,outside radii of fiber radial layers,rho_max]]         
numlayers = numel(RL)-1;
resol = ones(1,numlayers)*.001;
       % [MM : radial sampling resolution in nm over each layer]
index = lt(diff(RL),.01);   % find if any layer is less than 10 nm thick
resol(index) = .0001;       % sets resol to 0.1 nm in layers less than 10 nm thick
%NOTE: loss corresponds to negative imaginary part of the index
resol
material{1}='SiO2-GeO2';
material{2}='SiO2-GeO2';
material{3}='Water1_315';
%material{4}='Air';
x=zeros(1,numlayers);                                     
GeConc = 0.03485;                   % [JA : refractive indices of the fiber layers (x(1) is the Ge mole fraction in the core )
x(1)=GeConc;                     % the other elements of x could be used to pass on more variables for the  materials

if 1                                    
%Calculation of the mask period corresponding to the Bragg wavelength (based on Neff of the core at that wavelength)
MaskPer = BraggOnly(SpectrumName,BraggWav,RL,resol,material,x,thermoOpticCoeff,temperature,count);
MaskPer = MaskPer + MaskPer*thermExpCoeff*temperature(count);

end
ncoreB=nLB  %to display the core refractive index (at the Bragg) and save it in the workspace in full resolution

if 1
%Now finding the core mode index and core fields at the new calculation wavelength
Neffcore = CoreOnly(SpectrumName,CalcWav,RL,resol,material,x,thermoOpticCoeff,temperature,count);
end 
ncoreC=nLC  %to display the core refractive index (at the Calc) and save it in the workspace in full resolution
if 1
%Now calculate a limited number of modes near a certain value of effective index
%Store the values as before and call the function for the coupling
%coefficients
Neffmax = (1+.00)*((2*CalcWav/MaskPer)-Neffcore)   %will be used to determine Umax for the eigenvalue solver. 
                                                    % turns out Umax is the "center" value of the search....
[dummy0, Neff] = Modes(SpectrumName,CalcWav,Mmax,Neig,Neffmax);
end

% ----- DETERMINING CUTOFF WAVELENGTH ----- %

disp("Modes computed.");

Neff = sort(Neff);

[NeffDiff, NeffCladCutoffIndex] = min(abs(Neff - currentOutsideIndex));
NeffCladCutoff = Neff(NeffCladCutoffIndex);

cutoffWavelength = (Neffcore + NeffCladCutoff)*MaskPer/2;

disp(cutoffWavelength);
disp("-----------------------------------------------------");

% ---------------------------------------- %

if 1
%Calculation of the coupling coefficients between the core mode found
%earlier and the modes just calculated (if the core mode is not part of the
%modes just calculated, it will not show up as a resonance in the spectrum
dummy1 = CCoeff(SpectrumName,MaskPer,theta_int,Mmax)
end

%Finally calculation of the spectrum (the wavelength span must correspond
%to the effective indices chosen in "Modes", i.e. around Neffmax)

[spectrumResponse, trueWavelengths, trueTroughs] = Spectrum(Neffcore,MaskPer,Length,Dn_AC,SpectrumName,InputPol,Npoints,WavStart,WavEnd,maxdetuning,cutoffWavelength);

troughs{1, count} = trueTroughs;
wavelengths{1, count} = trueWavelengths;

end