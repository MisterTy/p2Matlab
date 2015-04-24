function [Gr] = PiRegler( Gs,w,phir,k,T )
%PIREGLER Summary of this function goes here
%   Detailed explanation goes here

    %Berechnung von wpi/Tn
    phi_s=angle(Gs);                 %Phasengang Strecke
    for y=1:1:length(phi_s)             % Entfernt Sprung bei -pi
        if phi_s(y) > 0
            phi_s(y) = -2*pi+phi_s(y); 
        end
    end
    
    [ind_left,ind_right] = int_ver(phi_s,-pi/2);    %Berechnet die beiden Indizien die ges. wpi einschliessen
    wpi = (w(ind_left)+w(ind_right))/2;             %Weil Wpi irgendwo zwischen W[Left] und W[Right] liegt nehmen wir den arithm Mittelwert
    Tn=1/wpi
    Grp = (1 + (1./(Tn.*j*w)));
    
    %Durchtrittspunkt wD bestimmen
    phi_rprov= angle(Grp);              %Phasengang des prov. Pi-Reglers mit Kr=1
    phi_O = phi_s+phi_rprov;            %Phasengang des offenen Regelkreises            
    for y=1:1:length(phi_O)             % Entfernt Sprung bei -pi
        if phi_O(y) > 0
            phi_O(y) = -2*pi+phi_O(y); 
        end
    end
    
    [ind_left,ind_right] = int_ver(phi_O,-pi+phir);    %Berechnet die beiden Indizien die ges. wD einschliessen
    wD = (w(ind_left)+w(ind_right))/2;                 %Weil WD irgendwo zwischen W[Left] und W[Right] liegt nehmen wir den arithm Mittelwert
    
% Provisorische Übertragungsfunktion mit wD
    Grp_wd = 1*(1+1/(j*wD*Tn));
    
    % KR bestimmen 
    ampl_s = abs(Gs);                       %Amplitudengang Strecke
    ampl_rprov = abs(Grp_wd);               %Amplitudengang des prov. Pi-Reglers mit Kr=1
    ampl_O = ampl_s.*ampl_rprov;            %Amplitudengangs des offenen Regelkreises
    
    amplOwd = 20*log10((ampl_O(ind_left)+ampl_O(ind_right))/2);    %Amplitude bei wD
    KrdB=-amplOwd;                          %Reglerverstärkung in DB
    Kr=10^(KrdB/20)                         %Reglerverstärkung in DB
    Gr=Kr*(1+1./(1j.*w.*Tn));               %Übertragungsfunktion Regler

end

% Ende File ---------------------------------------------------------------

