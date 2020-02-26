data = importdata('0.txt');
wavelengths = data(:,1);
intensity = data(:,2);

figure;
plot(wavelengths, intensity, 'r');
grid on;