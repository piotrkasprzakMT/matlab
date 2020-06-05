% Przedmiot: Techniki Obliczeniowe 
% Kierunek studiów: Mechatronika 
% Semestr: 2
% Rok akademicki: 2019/2020
% Data (dzień-miesiąc-rok): <<05-06-2020>>
%
% Imię:             <<Piotr>>
% Nazwisko:         <<Kasprzak>>
% Numer albumu ZUT: <<46748>>

f1 = @(x)   2 + 1 ./ (1 + x.^2);
f2 = @(x)   2 + sin(x);

a = 0;
b = 5;
N = 1000;

x = linspace(a, b, N);

figure(1);
clf;
plot(x, f1(x));
grid on;
grid minor;
hold on;
plot(x, f2(x));

f21 = @(x) f2(x) - f1(x);

for i = 1:(N-1)
    try
        xL = fzero(f21, [x(i), x(i+1)]);
        break
    catch
    end
end

for i = N:-1:2
    try
        xR = fzero(f21, [x(i-1), x(i)]);
        break
    catch
    end
end

xLR = linspace(xL, xR, N);
xRL = flip(xLR);

yLR = f2(xLR);
yRL = f1(xRL);

xi = [xLR, xRL];
yi = [yLR, yRL];

fill(xi, yi, 'y');

text(xL,2.4,'\uparrow x=0.70; y=2.66 ')
text(xR,f1(xR),'\leftarrow x=3.04; y=2.10')

format long;
quadgk(f21, xL, xR)
