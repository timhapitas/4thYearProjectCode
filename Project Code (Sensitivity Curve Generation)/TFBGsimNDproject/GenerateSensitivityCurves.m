WavelengthData = load('Wavelengths');
TroughData = load('Troughs');

Wavelengths = WavelengthData.wavelengths;
Troughs = TroughData.troughs;
refractiveIndices = [1.315 1.320 1.325 1.330 1.335 1.340];

wavelengthsForPlot = cell(1, (max(size(Wavelengths))));
fittedCurves = cell(1, max(size(Wavelengths)));
sensCurves = cell(1, (max(size(Wavelengths))));
indexDifferences = zeros(1, (max(size(Wavelengths))));
% indexDiff = 0.005;


for i = 1:max(size(Wavelengths))
    
    [dummy, indices] = find(diff(Wavelengths{1,i}) == min(diff(Wavelengths{1,i})));
    
    if (abs(Wavelengths{1,i}(indices) - Wavelengths{1,i}(indices + 1)) <= 0.0001)
        Wavelengths{1,i}(indices) = [];
        Troughs{1,i}(indices) = [];
    end
        
    if i > 1
            
%         lengthDiff = abs(length(Wavelengths{1, i}) - length(Wavelengths{1, i-1}));
        lengthDiffPureWater = abs(length(Wavelengths{1, i}) - length(Wavelengths{1, 1}));
        pureWaterWavelengths = Wavelengths{1, 1}(lengthDiffPureWater+1:end);
        wavelengthsForPlot{1, i} = pureWaterWavelengths;
        pureWaterTroughs = Troughs{1, 1}(lengthDiffPureWater+1:end);
%         Wavelengths{1,i-1} = Wavelengths{1, i-1}(lengthDiff+1:end);
        wavelengthDiff = abs(Wavelengths{1, i} - pureWaterWavelengths);
        indexDifferences(i) = abs(refractiveIndices(i) - refractiveIndices(1)); 
        sensCurves{1, i} = wavelengthDiff/indexDifferences(i);
        
        currentSensCurve = sensCurves{1,i};
        
        % Set up fittype and options.
        [xData, yData] = prepareCurveData(pureWaterWavelengths, sensCurves{1,i});
        ft = fittype('power1');
        opts = fitoptions('Method', 'NonlinearLeastSquares');
        opts.Display = 'Off';
        opts.StartPoint = [2.31449424555992e+35 -191.333217731478];

        % Fit model to data.
        [fitresult, gof] = fit(xData, yData, ft, opts);
        
        fitParams = coeffvalues(fitresult);
        fittedSensCurve = @(x) fitParams(1).*x.^(fitParams(2));
        fittedCurves{1,i} = fittedSensCurve;
        
        figure('Renderer', 'painters', 'Position', [100 100 900 600]);
        plot(wavelengthsForPlot{1, i}, sensCurves{1, i}, 'r.');
        hold on;
        fplot(fittedSensCurve, [(min(pureWaterWavelengths)-0.0025) (max(pureWaterWavelengths)+0.0025)], 'b');
        xlim([(min(pureWaterWavelengths)-0.0025) (max(pureWaterWavelengths)+0.0025)]);
        title('Sensitivity of Guided Resonances to Change in External Index of ' + string(abs(refractiveIndices(i) - refractiveIndices(1))) + '(Relative to Pure Water)');
        xlabel('Guided Mode Wavelengths in Pure Water (pm)');
        ylabel('Sensitivity of Guided Modes (\Delta \lambda / \Delta SRI)');
        legend('Modes Tracked in Simulation', 'Fitted Function of the Form ax^b (a = ' + string(fitParams(1)) + ', b = ' + string(fitParams(2)) + ')');
        hold off;
        grid on
        
    end
    
end

% generate plots in 3D space

figure('Renderer', 'painters', 'Position', [100 100 900 600]);
h = axes;
z = transpose(indexDifferences)*ones(1, 1000);
ylim([min(Wavelengths{1,1}) max(Wavelengths{1,1})]);

for j = 1:max(size(sensCurves))
    
    if j == 1
        continue;
        
    else
        x = linspace(min(wavelengthsForPlot{1,j}), max(wavelengthsForPlot{1,j}), 1000);
        y = fittedCurves{1,j}(x);
        ztemp = z(j,:);
        
        plot3(ztemp, x, y);
        set(h,'ydir','reverse');
        grid on;
        hold on;
        
    end
    
end

hold off;
title('3D View of Fitted Sensitivity Curves for Increasing \Delta SRI (Not Yet Interpolated)');
xlabel('Increasing \Delta SRI');
ylabel('Guided Mode Resonance Wavelengths (pm)');
zlabel('Sensitivity (\Delta \lambda/\Delta SRI)');


