f_samp = 425e3;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-7.95) / (2.285*0.02353*pi));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 13;

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) - ideal_lp(0.826*pi,n) + ideal_lp(0.614*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
plot(f,abs(H))
grid

FIR_BandStop

function hd = ideal_lp(wc,M);

alpha = (M-1)/2;
n = [0:1:(M-1)];
m = n - alpha + eps;
hd = sin(wc*m) ./ (pi*m);
end