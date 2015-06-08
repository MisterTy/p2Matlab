% Autor:        Pascal Ackermann, Michael Bos
%***************************************************************
%Dieses M-File berechnet und plottet die Schrittantwort des PI- sowie des 
%PID-Reglers. Zum Vergleich wird die Berechnung einmal mit der
%Step-Funktion durchgeführt und einmal mit der FFT.
%Dieses M-File wird vom File "phasengangmethode_sa_Prozent_Ueberschw.m" so
%oft ausgeführt bis sich das gewünschte Überschwingen einstellt.
%***************************************************************
%Übertragungsfunktion
%Param:
%       kR      = Regelverstärkung
%       Tm      = Nachstellzeit
%       Tv      = Vorstellzeit
%       Tp      = parasitäre Zeitkonstante
%       kS      = Streckenbeiwert
%       T       = Zeitkonstanten der Strecke
%       w       = Kreisfrequenzspektrum
%       typ     = "PI"/"PID" als String
%Param:
%       kR      = Regelverstärkung
%       Tm      = Nachstellzeit
%       Tv      = Vorstellzeit
%       Tp      = parasitäre Zeitkonstante
%       kS      = Streckenbeiwert
%       T       = Zeitkonstanten der Strecke
%       w       = Kreisfrequenzspektrum
%       typ     = "PI"/"PID" als String
%Rückgabe:
%       prozUenSchw      = Überschwingen in Prozent
%       t                = Zeitintervall
%       y1               = Amplitude der Schrittantwort

function [prozUebSchw, t,y1] = schrittantwort_Prozent_Ueberschw( kR, Tn, Tv, Tp, kS, T, w, typ)
    if (strcmp(typ,'Pi'))
        %--------------------------Mit Step--------------------------------
        %Zähler- und Nennerpolynom der Strecke
        ZaehlerPolynomB_Strecke=kS;        
        N=1;                               %Initialisiere N
        s= tf('s');
        for y=1:1:length(T)                 %ohne kS
            N = N*(1+s*T(y));
        end
        NennerPolynomA_Strecke = N;
        %Zähler- und Nennerpolynom des Reglers
        ZaehlerPolynomB_Regler = kR*(s*Tn+1);   
        NennerPolynomA_Regler = s*Tn;
        %Uebertragungsfunktion geschl. Regelkreis
        UET= (ZaehlerPolynomB_Strecke*ZaehlerPolynomB_Regler)/((NennerPolynomA_Strecke*NennerPolynomA_Regler)+(ZaehlerPolynomB_Strecke*ZaehlerPolynomB_Regler));
        clf('reset');           %Damit nur letzte Iteration angezeigt wird
        figure(1);
        step(UET);hold on;
        %--------------------------Mit FFT--------------------------------
        %Zähler- und Nennerpolynom der Strecke
        ZaehlerPolynomB_Strecke=kS;   
        N = [T(1),1];               % Initialisiere N = (sT1 +1) f?r conv
        for y=2:1:length(T)

              if(T(y)==0)           
                  break
              end
              N = conv(N,[T(y) 1]);  % Multipliziert die Terme des Nenners
        end
        NennerPolynomA_Strecke = N; 
        %Zähler- und Nennerpolynom des Reglers
        ZaehlerPolynomB_Regler=kR*[Tn 1];  
        NennerPolynomA_Regler= [Tn 0];      
        %Uebertragungsfunktion geschl. Regelkreis
        B = conv(ZaehlerPolynomB_Strecke, ZaehlerPolynomB_Regler);
        A = conv(NennerPolynomA_Strecke, NennerPolynomA_Regler);
        A(length(A)-length(B)+1:length(A)) = A(length(A)-length(B)+1:length(A)) + B;
        %FFT durchführen
        fs= max(w);                         %Abtastfrequent
        N= 2^floor(log2(length(w)));        %Anzahl Punkte im Kreisfrequenz w-Array
        H = freqs(B, A, w); 
        %Schrittantwort berechnen
        [y1, t] = schrittIfft(B,A,fs,N);    % Methode via ifft().
        %Plot
        plot(t, y1,'red'),grid on
        %Überschwingen berechnen
        maximum = max(y1);
        endwert = y1(N);
        prozUebSchw = ((100/endwert)*maximum)-100;
    end

    if (strcmp(typ,'Pid'))
        %--------------------------Mit Step--------------------------------   
        %Zähler- und Nennerpolynom der Strecke
        ZaehlerPolynomB_Strecke=kS;        
        N=1;                               %Initialisiere N
        s= tf('s');
        for y=1:1:length(T)                 %ohne kS
            N = N*(1+s*T(y));
        end
        NennerPolynomA_Strecke = N;
        %Zähler- und Nennerpolynom des Reglers
        ZaehlerPolynomB_Regler = kR*(s*Tn*(1+s*Tp)+(1+s*Tp)+(s^2)*Tn*Tv);   
        NennerPolynomA_Regler = s*Tn+(s^2)*Tn*Tp;
        % Uebertragungsfunktion geschl. Regelkreis
        UET= (ZaehlerPolynomB_Strecke*ZaehlerPolynomB_Regler)/((NennerPolynomA_Strecke*NennerPolynomA_Regler)+(ZaehlerPolynomB_Strecke*ZaehlerPolynomB_Regler));
        clf('reset')%Damit nur letzte Iteration angezeigt wird
        figure(1);
        step(UET);hold on;
    
        %--------------------------Mit FFT--------------------------------
        %Zähler- und Nennerpolynom der Strecke
        ZaehlerPolynomB_Strecke=[kS];        
        N = [T(1),1];               % Initialisiere N = (sT1 +1) f?r conv
        for y=2:1:length(T)

              if(T(y)==0)           
                  break
              end
              N = conv(N,[T(y) 1]);  % Multipliziert die Terme des Nenners
        end
        NennerPolynomA_Strecke = N; 
        %Zähler- und Nennerpolynom des Reglers
        ZaehlerPolynomB_Regler=[kR*(Tn*Tp+Tn*Tv) kR*(Tn+Tp) kR];  
        NennerPolynomA_Regler= [Tn*Tp Tn 0];  % Keine Konstante daher 0    
        % Uebertragungsfunktion Geschlossener Regelkreis 
        B = conv(ZaehlerPolynomB_Strecke, ZaehlerPolynomB_Regler);
        A = conv(NennerPolynomA_Strecke, NennerPolynomA_Regler);
        A(length(A)-length(B)+1:length(A)) = A(length(A)-length(B)+1:length(A)) + B;

        %FFT durchführen
        fs= max(w);                         %Grenzfrequenz
        N= 2^floor(log2(length(w)));        % Anzahl Punkte im Kreisfrequenz w-Array
        % Schrittantwort berechnen
        [y1, t] = schrittIfft(B,A,fs,N);    % Methode via ifft().
        %Plot
        plot(t, y1,'red'),grid on
        %Überschwingen berechnen
        maximum = max(y1);
        endwert = y1(N);
        prozUebSchw = ((100/endwert)*maximum)-100;
    end
end
