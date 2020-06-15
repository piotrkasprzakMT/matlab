% Przedmiot: Techniki Obliczeniowe 
% Kierunek studiów: Mechatronika 
% Semestr: 2
% Rok akademicki: 2019/2020
% Data (dzień-miesiąc-rok): <<05-06-2020>>
%
% Imię:             <<Piotr>>
% Nazwisko:         <<Kasprzak>>
% Numer albumu ZUT: <<46748>>

% Obliczanie PSD z pomocą FFT

% Czytanie wartości y, 
% rekonstrukcja wartości x.
%
%[y, fs] = audioread('data/labrador-barking-daniel_simon.wav');

fsignal = 350; % częstotliwość sygnału, Hz
fs = 8.23*pi * fsignal;

t1 = 0.0; % początek, czas w sekundach
t2 = 2.0; % koniec, czas w sekundach
N_SAMPLES = (t2 - t1) * fs;
t = linspace(t1, t2, N_SAMPLES);

y = sin(2 * pi * fsignal .* t);
y = y.';

% Tylko 1 kanał
%
y = y(:,1);

y = y(1:2:end)

%sound(y, fs);

N = length(y);
Delta = 1 ./ fs; 
x = (0:(N-1))' .* Delta;

% Rysowanie danych wejściowych
%
figure(1);
clf;
subplot(2, 1, 1); 
plot(x,y); 
title('dane'); 
xlabel('t [sekundy]')



% Jawne utworzenie skali częstotliwości.
%
f = (-N/2:N/2)' ./ (N .* Delta);

% Szybka transformacja Fouriera,
% mnożenie przez Delta jest konieczne
% jeżeli chcemy mieć dobre jednoski fizyczne.
%
F = fft(y) .* Delta;

% Przetasowanie wyników tak, aby przebiegały
% od najmniejszej wartości f, a nie od zera.
% Inaczej trochę dla parzystych/nieparzystych N.
%
if mod(N,2) == 0
  F = fftshift(F); % parzyste N
  F = [F; F(1) ]; 
else
  F = [F; F(1) ];  % nieparzyste N
  F = fftshift(F);
end

% PSD liczymy według wzoru:
%
%   p  = abs(F).^2 + abs(flipud(F)).^2;
%
% Poniżej jest nieco szybszy sposób obliczenia
% dający te same wyniki.
%
p  = 2.0*abs(F).^2;

% Zostawiamy tylko PSD dla wartości dodatnich f.
% To samo robimy ze skalą częstotliwości f.
%
p  = p(f >= 0);
pf = f(f >= 0);

% Rysujemy PSD po transformacji.
%
subplot(2, 1, 2); 
semilogx(pf, p, '-ro'); 
title('PSD'); 
xlabel('f [Hz]'); 
xlim([1, 10000]);
grid on;
grid minor;

