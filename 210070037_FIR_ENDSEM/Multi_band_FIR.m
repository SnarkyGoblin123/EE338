f_samp = 600e3;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-7.95) / (2.285*0.01667*pi));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 23;

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) - ideal_lp(0.725*pi,n) + ideal_lp(0.2583*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
% plot(f,abs(H))
grid


f_samp = 600e3;

%Band Edge speifications


%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-7.95) / (2.285*0.01667*pi));           %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 23;

%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(0.8417*pi,n) - ideal_lp(0.1417*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, f_samp);
% plot(f,abs(H))
% plot(f,angle(H))
grid

MULTI_BAND_FIR = conv(FIR_BandStop, FIR_BandPass);
MULTI_BAND_FIR
fvtool(MULTI_BAND_FIR); 
[H,f] = freqz(MULTI_BAND_FIR,1,1024, f_samp);
plot(f,abs(H))
xlabel("frequency")
ylabel("Magnitude")
grid
% impz(MULTI_BAND_FIR)
function hd = ideal_lp(wc,M);

alpha = (M-1)/2;
n = [0:1:(M-1)];
m = n - alpha + eps;
hd = sin(wc*m) ./ (pi*m);

end