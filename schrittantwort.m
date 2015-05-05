function [] = schrittantwort( kR, Tn, Tv, Tp, kS, T, w, typ)

%Dieses M-File berechnet und plottet die Schrittantwort des PI- sowie des 
%PID-Reglers

if (strcmp(typ,'Pi'))
    
    %%%Mit Step
    
    %Strecke
    ZaehlerPolynomB_Strecke=kS;        
    N=1;                               %Initialisiere N
    s= tf('s');
    for y=1:1:length(T)                 %ohne kS
        N = N*(1+s*T(y));
    end
    NennerPolynomA_Strecke = N;
    %Regler
    ZaehlerPolynomB_Regler = kR*(s*Tn+1);   
    NennerPolynomA_Regler = s*Tn;
    % Uebertragungsfunktion geschl. Regelkreis
    UET= (ZaehlerPolynomB_Strecke*ZaehlerPolynomB_Regler)/((NennerPolynomA_Strecke*NennerPolynomA_Regler)+(ZaehlerPolynomB_Strecke*ZaehlerPolynomB_Regler));
    step(UET);
    
    %%Mit FFT
    %Strecke
    ZaehlerPolynomB_Strecke=kS;        
%     N=1;                               %Initialisiere N
%     syms s;
%     for y=1:1:length(T)                 %ohne kS
%         N = N*(1+s*T(y));
%     end

% Berechnung des Nennerpolynoms 
    N = [T(1),1];               % Initialisiere N = (sT1 +1) für conv
    for y=2:1:length(T)
          
          if(T(y)==0)           
              break
          end
          N = conv(N,[T(y) 1]);  % Multipliziert die Terme des Nenners
    end
    
    
    NennerPolynomA_Strecke = N; %sym2poly(N);
    %Regler
    ZaehlerPolynomB_Regler=kR*[Tn 1];  
    NennerPolynomA_Regler= [Tn 0];      
    % Übertragungsfunktion Geschlossener Regelkreis ***************************

    B = conv(ZaehlerPolynomB_Strecke, ZaehlerPolynomB_Regler);
    A = conv(NennerPolynomA_Strecke, NennerPolynomA_Regler);
    A(length(A)-length(B)+1:length(A)) = A(length(A)-length(B)+1:length(A)) + B; %Keine Ahnung was er da macht

    % Plot Schrittantwort Geschlossener Regelkreis ****************************
    
    %Berechnung fs , N
    fs= max(w);                         %Grenzfrequenz
    N= 2^floor(log2(length(w)));        % Anzahl Punkte im Kreisfrequenz w-Array
    
    % Schrittantwort berechnen
    [y1, t] = schrittIfft(B,A,fs,N);    % Methode via ifft().
    SchrittantwortStrecke=figure(1);
    set(SchrittantwortStrecke,  'name', 'Schrittantwort Strecke',...        % Setzt Titel des Fensters neu
                                'numbertitle', 'off');                      % Lässt 'Figure 1' verschwinden

    hold on;                 
    plot(t, y1,'red'),grid on
end

if (strcmp(typ,'Pid'))
       
    %%%Mit Step

    %Strecke
    ZaehlerPolynomB_Strecke=kS;        
    N=1;                               %Initialisiere N
    s= tf('s');
    for y=1:1:length(T)                 %ohne kS
        N = N*(1+s*T(y));
    end
    NennerPolynomA_Strecke = N;
    %Regler
    ZaehlerPolynomB_Regler = kR*(s*Tn*(1+s*Tp)+(1+s*Tp)+(s^2)*Tn*Tv);   
    NennerPolynomA_Regler = s*Tn+(s^2)*Tn*Tp;
    % Uebertragungsfunktion geschl. Regelkreis
    UET= (ZaehlerPolynomB_Strecke*ZaehlerPolynomB_Regler)/((NennerPolynomA_Strecke*NennerPolynomA_Regler)+(ZaehlerPolynomB_Strecke*ZaehlerPolynomB_Regler));
    step(UET);
    
    %%Mit FFT
    %Strecke
    ZaehlerPolynomB_Strecke=[kS];        

% Berechnung des Nennerpolynoms 
    N = [T(1),1];               % Initialisiere N = (sT1 +1) für conv
    for y=2:1:length(T)
          
          if(T(y)==0)           
              break
          end
          N = conv(N,[T(y) 1]);  % Multipliziert die Terme des Nenners
    end
    
    NennerPolynomA_Strecke = N; 
    %Regler
    ZaehlerPolynomB_Regler=[kR*(Tn*Tp+Tn*Tv) kR*(Tn+Tp) kR];  
    NennerPolynomA_Regler= [Tn*Tp Tn 0];  % Keine Konstante daher 0    
    % Übertragungsfunktion Geschlossener Regelkreis ***************************

    B = conv(ZaehlerPolynomB_Strecke, ZaehlerPolynomB_Regler);
    A = conv(NennerPolynomA_Strecke, NennerPolynomA_Regler);
    A(length(A)-length(B)+1:length(A)) = A(length(A)-length(B)+1:length(A)) + B; %Keine Ahnung was er da macht

    % Plot Schrittantwort Geschlossener Regelkreis ****************************
    
    %Berechnung fs , N

    fs= max(w);                         %Grenzfrequenz
    N= 2^floor(log2(length(w)));        % Anzahl Punkte im Kreisfrequenz w-Array
    
    % Schrittantwort berechnen
    [y1, t] = schrittIfft(B,A,fs,N);    % Methode via ifft().
    SchrittantwortStrecke=figure(1);
    set(SchrittantwortStrecke,  'name', 'Schrittantwort Strecke',...        % Setzt Titel des Fensters neu
                                'numbertitle', 'off');                      % Lässt 'Figure 1' verschwinden

    hold on;                 
    plot(t, y1,'red'),grid on
    
    
end
end

