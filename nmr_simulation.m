%Pauli Matrices & identity
sigma_x = [0 , 1/2 ; 1/2 , 0];
sigma_y = [0, -1i/2; 1i/2, 0];
sigma_z =[1/2, 0; 0, -1/2];
identity = [1 0; 0 1];


%Defining the 4x4 matrices for each nucleus
I_Ax = kron(sigma_x,identity);
I_Ay = kron(sigma_y, identity);
I_Az = kron(sigma_z, identity);
I_Bx = kron(identity,sigma_x);
I_By = kron(identity,sigma_y);
I_Bz = kron(identity,sigma_z);


%Frequencies in Hamiltonian
omega_0 = 100;
delta = 20;
omega_A = omega_0 + delta/2;
omega_B = omega_0 - delta/2;
omega_J = 15;

%Expression for 4x4 Hamiltonian matrix
H = omega_A*I_Az + omega_B*I_Bz + omega_J*(I_Ax*I_Bx + I_Ay*I_By + I_Az*I_Bz);

%Initial density matrix & observable operator
rho = I_Ax + I_Bx;
I_plus = (I_Ax + I_Bx) + 1i*(I_Ay + I_By);

%Time step at which rho evolves & sampling frequency (1/delta_T)
time_step = 0.25 / norm(H,2);
samp_freq = 1 / time_step;

%Time-evolved density matrix
p = expm(-1i * H * time_step);

%Observable signal
fid = zeros(1,2048);

%Evaluate observable signal by evolving density matrix
for n=1:2048
    fid(n) = trace(rho*I_plus);
    rho = p*rho*p';
end

%Signal is infinite, cut off the end using exponential
window_function = exp(-5*linspace(0,1,2048));
fid = fid .* window_function;

%Get Fourier transform of signal and centre at f=0 Hz
spectrum = fftshift(fft(fid,8192));

%Get frequency axis
freq = (-4096:4095)*(2*pi*samp_freq/8192);

%Plot FT of FID signal against frequency axis
plot(freq,real(spectrum),"LineWidth",2);
title("Frequency spectrum for two coupled spins", "with $\omega_0=100$, $\delta = 20$, and $J=15$", "FontWeight","bold","interpreter","latex")
xlabel("Frequency / $rad\:s^{-1}$","interpreter","latex", "FontWeight","bold")
ylabel("Signal amplitude","interpreter","latex", "FontWeight","bold")
set(gca, 'FontSize', 14)
axis([50 150 0 200])