% Beschreibung:     Dieses File berechnet die grundlegenden Parameter zur
%                   Dimensionerung eines PI-Reglers. Diese Dimensionierung
%                   erfolgt aufgrund einer Regelstrecke, deren
%                   Übertragungsfunktion (Verstärkung und Zeitkonstanten)
%                   bekannt ist.
% Eingabeparameter: Gs,w,phir,kS,T 
%                   - Gs :  Übertragungsfunktion der Strecke als Array
%                   - w:    Kreisfrequenzspektrum in Listenform
%                   - phir: gewünschter Phasenrand (->Überschwingen)
%                   - kS:   Verstärkung der Strecke
%                   - T:    Zeitkonstanten der Strecke als Array
%
% Ausgabeparameter: - kR:   Reglerverstärkung
%                   - Tn:   Nachstellzeit des Reglers
%                   - Tv:   Vorstellzeit des Reglers
%                   - Tp:   Parasitäre Zeitkonst. des Reglers
%
% Autor:            Michael Bos,Pascal Ackermann
% Datum:            28.04.2015

function [kR, Tn, Tv, Tp] = PiRegler(Gs,w,phir,kS,T)

    %Berechnung von Tn-----------------------------------------------------
    phi_s=angle(Gs);                    % Phasengang Strecke
    
    for y=1:1:length(phi_s)             % Entfernt Sprung bei -pi
        if phi_s(y) > 0                 %--> Punktsuche funktioniert sonst nicht
            for z=y:1:length(phi_s)
                phi_s(z) = -2*pi+phi_s(z); 
            end
            break
        end
    end
    
    [ind_left,ind_right] = int_ver(phi_s,-pi/2);% berechnet Indizien, die ges. wpi einschliessen
    wpi = (w(ind_left)+w(ind_right))/2;         %Weil Wpi irgendwo zwischen W[Left] und W[Right] liegt nehmen wir den arithm Mittelwert
    Tn=1/wpi;               
    Grp = (1 + (1./(Tn.*1j*w)));                %Übertragungfk. :Pi-Regler mit kR=1
    
    %Durchtrittspunkt wD bestimmen
    phi_rprov= angle(Grp);              % Phasengang: Pi-Regler mit kR=1
    phi_O = phi_s+phi_rprov;            % Phasengang des offenen Regelkreises 
    
    for y=1:1:length(phi_O)             % Entfernt Sprung bei -pi
        if phi_s(y) > 0                 %--> Punktsuche funktioniert sonst nicht
            for z=y:1:length(phi_O)
                phi_O(z) = -2*pi+phi_O(z); 
            end
            break
        end
    end
    
    [ind_left,ind_right] = int_ver(phi_O,-pi+phir);    %Berechnet die beiden Indizien die ges. wD einschliessen
    wD = (w(ind_left)+w(ind_right))/2;                 %Weil WD irgendwo zwischen W[Left] und W[Right] liegt nehmen wir den arithm Mittelwert
    
    % Provisorische Übertragungsfunktion bei wD
    Grp_wd = 1*(1+1/(1j*wD*Tn));
    
    % G(s) Strecke bei Wd
    G_str_wd=1;                         %Initialisiere Gs
    for y=1:1:length(T)                 %Berechne Übertragungsfunktion    
        G_str_wd = G_str_wd.*(1./(1+1j.*wD.*T(y)));
    end
    G_str_wd=kS*G_str_wd;               %Vervollständige mit Verstärkung 
    
    % G(s) offener Regelkreis bei Wd
    GOffwd=Grp_wd*G_str_wd; 
    
    %Reglerverstärkung in dB
    amplOffwd=20*log10(abs(GOffwd));    
    kRdB=-amplOffwd;  
    %Reglerverstärkung linear
    kR=10^(kRdB/20);                   
    
    %Parameter die beim PI-Regler nicht vorkommen
    Tv=0;
    Tp=0;
    
end

% Ende File ---------------------------------------------------------------

