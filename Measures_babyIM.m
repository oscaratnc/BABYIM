clear all 
clc

%Read data from csv File Measurements
filename ='C:\Users\danos\Dropbox\ITH\BABYIM\DataFile_Read.csv';
M = csvread(filename);
RedOriginal = M(1,:);
IRoriginal = M(2,:);

SampleRate = 1000;
NyqRate = SampleRate/2;

%Extract DC Component of both signals for SpO2 calculation
DCRed = mean(RedOriginal);
DCIR = mean(IRoriginal);

%Eliminate DC Component and Normalization of the signal
 RedAC = RedOriginal-DCRed;
 IRAC = IRoriginal - DCIR;
 
RedACN = RedAC/max(abs(RedAC));
IRACN = IRAC/max(abs(IRAC));

MMfiltRed = movmedian(RedACN,250);
MMfiltIR = movmedian(IRACN,250);

fs = SampleRate;
i=0;
while  i <= 2
if i == 1
    x1= MMfiltRed;
else
    x1= MMfiltIR;
end

n = (200 * fs)/1000;
ki = 1;
for k=1:length(x1)
    linea0 = zeros(1,n);
    lim1 = k - (n/2);
    lim2 = k + (n/2) - 1;
    if lim2 > length(x1)
        lim2 = length(x1);
    end
    if lim1 <= 0
        for k2 = (n/2) + lim1: lim2
            linea0(k2) = x1(k);
            k=k+1;
        end
    else
        kl = 1;
        for k2 = lim1: lim2
            linea0(kl) = x1(k2);
            kl=kl+1;
        end 
    end
    linea0 = sort(linea0);
    if lim1 <= 0 
        prom = (linea0(n+lim1-1)+linea0(n+lim1))/2;
    else
        prom = (linea0(n/2)+linea0(n/2 +1))/2;
    end
    linea1(ki) = prom;
    ki=ki+1;
end
x1 = x1-linea1;
% ventana 600mseg
% n=216;
% Cálculo de número de muestras para ventana de 600ms
n = (600 * fs)/1000;
ki = 1;
for k=1:length(x1)
    linea0 = zeros(1,n);
    lim1 = k - (n/2);
    lim2 = k + (n/2) - 1;
    if lim2 > length(x1)
        lim2 = length(x1);
    end
    if lim1 <= 0
        for k2 = (n/2) + lim1: lim2
            linea0(k2) = linea1(k);
            k=k+1;
        end
    else
        kl = 1;
        for k2 = lim1: lim2
            linea0(kl) = linea1(k2);
            kl=kl+1;
        end 
    end
    linea0 = sort(linea0);
    if lim1 <= 0 
        prom = (linea0(n+lim1-1)+linea0(n+lim1))/2;
    else
        prom = (linea0(n/2)+linea0(n/2+1))/2;
    end
    linea2(ki) = prom;
    ki=ki+1;
end
 
 x1=x1-linea2;
 
 if i == 1
    RedBLD= x1;
else
    IRBLD= x1;
end
 i = i+1;
end 


%medFiltRed = medfilt1(RedOriginal);
%medFiltIR = medfilt1(IRoriginal);


% bdr = 30;
% opr = 40;
fc= 6;
% Wp = fc/NyqRate;
% Ws=((.1*fc)+fc)/NyqRate;
% [N,Wn] = buttord(Wp,Ws, bdr,opr);
% [b,a] = butter(N,Wn,'low');

N=50;
fc2= 60;
Wn = fc/NyqRate;
WnL  = (fc2-.5)/NyqRate;
WnF = (fc2+.5)/NyqRate;
bandwidth = [WnL WnF];


A = fir1(100,[.5,5]/NyqRate);
B = fir1(N,bandwidth,'stop');

Filt60Red = filtfilt(B,1,RedBLD);
Filt60IR = filtfilt(B,1,IRBLD);


% FiltACRed = Filt60Red-DCRed;
% FiltACIR = Filt60IR - DCIR;

FilteredRed = filtfilt(A,1,Filt60Red);
FilteredIR  = filtfilt(A,1,Filt60IR);



Ts = 1/SampleRate;

FFTFilteredRed = fft(FilteredRed);
LRed = length(FilteredRed);
tRed = (0:LRed-1)*Ts;
P2Red = abs(FFTFilteredRed/LRed);
P1Red = P2Red(1:LRed/2+1);
P1Red(2:end-1) = 2*P1Red(2:end-1);
fRed = (SampleRate*(0:(LRed/2)))/LRed;

FFTFilteredIR = fft(FilteredIR);
LIR = length(FilteredIR);
tIR = (0:LIR-1)*Ts;
P2IR = abs(FFTFilteredIR/LIR);
P1IR = P2IR(1:LIR/2+1);
P1IR(2:end-1) = 2*P1IR(2:end-1);
fIR = (SampleRate*(0:(LIR/2)))/LIR;

SPO2 = round(100-(25*(FilteredRed/DCRed)/(FilteredIR/DCIR)));
x= ['SPO2: ', num2str(SPO2), '%'];

disp(x)


%Plot Original
subplot(3,2,1),plot(RedOriginal,'r'), title('Red Led Measurements');
subplot(3,2,2),plot(IRoriginal), title('IR Led Measurements');
subplot(3,2,3),plot(RedACN,'r'), title('Normalized Red Led ');
subplot(3,2,4),plot(IRACN), title('Normalized IR Led');
subplot(3,2,5),plot(MMfiltRed,'r'), title('MM Filter Red Led ');
subplot(3,2,6),plot(MMfiltIR), title('MM Filter IR Led');

figure
subplot(2,1,1),plot(RedBLD,'r'), title('BLD Red Led ');
subplot(2,1,2),plot(IRBLD), title('BLD IR Led');

figure
subplot(2,1,1),plot(Filt60Red,'r'), title('60Hz Stop Band Filter Red Led ');
subplot(2,1,2),plot(Filt60IR), title('60Hz Stop Band Filter IR Led');

figure
subplot(3,1,1),plot(FilteredRed,'r'), title('.5-5Hz BP Filter Red Led ');
subplot(3,1,2),plot(FilteredIR), title('.5-5Hz BP Filter Red Led ');

subplot(3,1,3),plot(FilteredRed,'r'), title('.5-5Hz BP Filter Red Led ');
%subplot(2,1,2),plot(fRed(1:200),P1Red(1:200)), title ('FFT Filtered Red Signal ');
%subplot(2,1,1),plot(RedOriginal,'r'), title('Red Led Measurements')
hold on 
plot(FilteredIR), title('.5-5Hz BP Filter IR Led');
%subplot(5,1,5),plot(IRoriginal), title('IR Led Measurements');
%subplot(2,1,2),plot(fIR(1:200),P1IR(1:200)), title ('FFT Filtered IR Signal ');



