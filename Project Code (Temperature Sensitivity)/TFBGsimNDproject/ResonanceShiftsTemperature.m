close all;
clear;
clc;

wavelengthData = load('wavelengths2.mat');
troughData = load('troughs2.mat');

wavelengths = wavelengthData.wavelengths;
troughs = troughData.troughs;

wavelengthShifts = cell(1, (max(size(wavelengths)) - 1));
temperature = [0 5 10 15 20];

figure;
plot(wavelengths{1,1}, troughs{1,1}, 'b.');
grid on;

% clean up data then compute resonance shifts for each step in temperature
for i = 1:max(size(wavelengths))

    tempIndices = find(wavelengths{1,i} < 1.59);
    deleteIndices = find(troughs{1,i}(tempIndices) > -8.0);
    troughs{1,i}(deleteIndices) = [];
    wavelengths{1,i}(deleteIndices) = [];
        
    if i == 1
        
        figure;
        plot(wavelengths{1,i}, troughs{1,i}, 'b.');
        title('Peaks Tracked in Simulation (Including Bragg Peak)');
        xlabel('Wavelength (\mu m)');
        ylabel('Transmission (dB)');
        grid on;
        
    else
        wavelengthShifts{1,i-1} = vpa(abs(wavelengths{1,i} - wavelengths{1,1}));
        wavelengthShifts{1,i-1} = double(wavelengthShifts{1,i-1});
          
    end
      
end

%Reshape shift data to make slope computation simpler
wavelengthShiftsReshape = zeros(length(wavelengthShifts{1,1}), max(size(wavelengthShifts)));
slopes = zeros(1, length(wavelengthShifts{1,1}));
slopeDistance = zeros(1, length(wavelengthShifts{1,1}));

for i = 1:length(wavelengthShifts{1,1})

    for j = 1:max(size(wavelengthShifts))
        
        wavelengthShiftsReshape(i,j) = wavelengthShifts{1,j}(i);
        
    end

end

% obtain slopes of shift vs temperature for all tracked modes. Plot slope
% vs distance from Bragg wavelength
temperature(1) = [];
for i = 1:length(slopes)
    
    coeffs = polyfit(temperature, wavelengthShiftsReshape(i,:), 1);
    slopes(i) = coeffs(1);
    slopeDistance(i) = (wavelengths{1,1}(end) - wavelengths{1,1}(i));
    
    if i == 1
        
        linFunc = @(x) slopes(i)*x + coeffs(2);
    
        figure;
        plot(temperature, wavelengthShiftsReshape(i,:), 'r.');
        hold on;
        fplot(linFunc);
        title('Resonance Shift Vs Temperature for \lambda = ' + string(wavelengths{1,1}(i)) + ' \mu m (Degrees C Above Room Temp)');
        xlabel('Temperature - Degrees Above Room Temp (C^\circ)');
        ylabel('Wavelength Shift (\mu m)');
        xlim([0 20])
        grid on;
    
    end
    
end

slopes = flip(slopes);
slopeDistance = flip(slopeDistance);

linFitResults = polyfit(slopeDistance, slopes, 1);
fittedFunc = @(x) linFitResults(1)*x + linFitResults(2);

figure;
plot(slopeDistance, slopes, 'r.');
hold on;
fplot(fittedFunc);
hold off;
title('Slopes of Linear Fits of Resonance Shifts Vs Temperature (\Delta\lambda / \DeltaT) (Unfinished)');
xlabel('Wavelength Distance from Bragg (\mum)');
ylabel('\Delta\lambda / \DeltaT (\mum / Degree)');
grid on;
xlim([slopeDistance(1), slopeDistance(end)]);
ylim([0.75e-6, 1.90e-6]);


