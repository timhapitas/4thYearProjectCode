function dummy = Spectrum(Neffcore,Mask_Per,Length,Dn_AC,SpectrumName,InputPol,Npoints,WavStart,WavEnd,maxdetuning)
%{

This program calculates the RUNGE-KUTTA matrix for the solution of Eqn. 128 (Michel)
The input mode is always the same, i.e. the HE11 (EVEN for now)!!
This is done in two steps: first define the vector of the unknown mode amplitudes
This vector is name "E" (for Eigenmodes) and give it an "initial" value at
z=L the at the end of the grating. E is made up of a "1" (the amplitude of the incoming
mode A), and (B1,B2,....,BNeig) the amplitudes of all the reflected modes
(including the fundamental, B1).  Therefore E has a dimension Neig+1.
The value of E at z=L is zero for all components, except 1 for the first.
This value (the amplitude of the forward mode at the end of the grating
should be an unknown, but it can be renormalized once we find its calculated value at
the beginning of the grating (where it should be 1).
This is because ALL initial conditions must be at the
same z.  The Runge-Kutta solution will then be run from the end towards the
beginning of the grating (see Alex B. Thesis page 84-ff).

This uses function handle @ME defined by:
            function DE = ME(z,E)

%}

global DnAC                                 
global MaskPer
global wavLoc                               % Transver of variables to f_MEcalc_V2.......
global Grat_length                           
global Mpm CpmS CpmP Npm modeNpm modeCpm modeMpm NewNmax Neff_core

Grat_length = Length;  %Grating length in microns
DnAC = Dn_AC;           %Index modulation amplitude
MaskPer = Mask_Per;
Neff_core = Neffcore;  %to be able to transfer to f_MEcalc_P through "global"
%=====================================================================                      

clear TdB wav transS transP TSpectrum


tic;   %just a timer to calculate how long it takes to get a grating.

DWav = (WavEnd-WavStart)/(Npoints-1);
zspan = [Grat_length,2500,0];                    %integration window in z from L to zero return solutions at these three z positions


figure
%Looping through wavelength spectrum........
if or(InputPol=='P',InputPol=='U')
    for Iwav = 1 : 1 : Npoints
        
    wavLoc = WavStart + (Iwav-1)*DWav;      %wavelength of current iteration, passed to the function f_MEcalc.m
    wav(Iwav,1) = wavLoc;                      %building a wavelength array for plotting spectrum with TdB.
    wavLoc;

%Now we select which modes will participate in the RK integration by checking phase matching.    
                                                  
    Detuning = abs((real(Neffcore+Npm))-2*wavLoc/MaskPer); %is a matrix containing the detunings of all modes at that wavelength
    IB= (Detuning<maxdetuning);                   %Keeping only modes within +/-maxdetuning in effective index from perfect phase matching
                                                  %need to increase detuning allowed to 0.01 to do lossy spectra (slows down the program by some...
    
    modeNpm=Npm(IB);  %Vectors of dimension IB keeping the resonant modes only
    modeCpm=CpmP(IB);
    modeMpm=Mpm(IB);
    NewNmax = length(modeCpm)+1;    %New number of modes in the matrix (A,B1,...Bmax)
    
    E_init =  zeros(NewNmax,1);                %gives an amplitude of 0 to all modes at grating output
    E_init(1,1)=1;                          %gives an amplitude of 1.0 for the core mode at the grating output
           
     [z,E0] = ode45(@f_MEcalc_P,zspan,E_init);
     % Now get the transmitted amplitude (absolute value) for P input
     transP(Iwav,1)=1/(abs(E0(3,1)^2));         %The first column of matrix E0 has the amplitude squared of the forward core mode at each calculated z position
                                            %The last line corresponds to the end of the integrationrange i.e. the input of the
                                            %grating. Since we had fixed the amplitude of this mode at 1 at the end of
                                            %the grating, renormalizing amplitudes so that the INPUTamplitude is 1, means that the
                                            %real output amplitude is 1/A1 and that we have to divide
                                            %all amplitudes by A1
    end
    TdBP=10*log10(transP);
    TSpectrumP = [wav,TdBP];
    save([SpectrumName,'/',SpectrumName,'_P'],'TSpectrumP','-v6');
    plot(wav,TdBP,'b');
    title('Transmission Spectrum of TFBG (Bragg Wavelength of 1610 nm, 0 Degree Tilt Angle)');
    xlabel('Wavelength (microns)');
    ylabel('Transmitted Power (dB)');
end

if or(InputPol=='S',InputPol=='U');
    figure
    for Iwav = 1 : 1 : Npoints;
        
     wavLoc = WavStart + (Iwav-1)*DWav;      %wavelength of current iteration, passed to the function f_MEcalc.m
     wav(Iwav,1) = wavLoc ;                     %building a wavelength array for plotting spectrum with Tmin.
     wavLoc;                                    %displays the iteration  wavelength

%Now we select which modes will participate in the RK integration by checking phase matching.    
                                                  
     Detuning = abs((real(Neffcore+Npm))-2*wavLoc/MaskPer); %is a matrix containing the detunings of all modes at that wavelength;   I am taking the real part of the detunings for this!
     IB= (Detuning<maxdetuning);                   %Keeping only modes within +/-0.001 in effective index from perfect phase matching
    
     modeNpm=Npm(IB);  %Vectors of dimension IB keeping the resonant modes only
     modeCpm=CpmS(IB);
     modeMpm=Mpm(IB);
     NewNmax = length(modeCpm)+1;    %New number of modes in the matrix (A,B1,...Bmax)
    
     E_init =  zeros(NewNmax,1);                %gives an amplitude of 0 to all modes at grating output
     E_init(1,1)=1;                          %gives an amplitude of 1.0 for the core mode at the grating output
           
     [z,E0] = ode45(@f_MEcalc_P,zspan,E_init);
     %Now get the transmitted amplitude for S input (absolute value)
     transS(Iwav,1)=1/(abs(E0(3,1)^2));         %The first column of matrix E0 has the amplitude of the forward core mode at each calculated z position
                                            %The last line corresponds to the end of the integrationrange i.e. the input of the
                                            %grating. Since we had fixed the amplitude of this mode at 1 at the end of
                                            %the grating, renormalizing amplitudes so that the INPUTamplitude is 1, means that the
                                            %real output amplitude is 1/A1 and that we have to divide
                                            %all amplitudes by A1
     
    end
    TdBS=10*log10(transS);
    TSpectrumS = [wav,TdBS];
    save([SpectrumName,'/',SpectrumName,'_S'],'TSpectrumS','-v6');
    plot(wav,TdBS,'b');
    title('S transmission');
end
toc
 
if InputPol=='U';
    figure
    TdBU=10*log10((transS+transP)/2);
    TSpectrumU = [wav,TdBU];
    save([SpectrumName,'/',SpectrumName,'_U'],'TSpectrumU','-v6');
    plot(wav,TdBU,'b');
    title('Unpolarized transmission');
end
dummy = 'Spectrum calc finished';
end





