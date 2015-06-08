% Autor:        Richard Gut
%***************************************************************
%Dieses M-File berechnet die Schrittantwort der mitgegebenen
%***************************************************************
%Übertragungsfunktion
%Param:
%       B   = Zählerpolynom
%       A   = Nennerpolynom
%       ws  = Abtastkreisfrequenz
%       N   = Blocklänge
%Rückgabe:
%       y   = Amplitude der Schrittantwort
%       t   = Zeitintervall der Schrittantwort

function [y, t] = schrittIfft(B,A,ws,N)

T = 1/ws;
w = linspace(0, pi*ws, N/2);
H = freqs(B, A, w); 
H = [H(1:N/2) 0 conj(H(N/2:-1:2))];
h = ifft(H); 
y = conv(h, ones(1,N+1));
y = y(1:length(y)/2)';
t = linspace(0, (length(y)-1)*T, length(y))';
