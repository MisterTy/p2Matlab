% Beschreibung:     Dieses File berechnet die grundlegenden Parameter zur
%                   Dimensionerung eines PID-Reglers. Diese Dimensionierung
%                   erfolgt aufgrund einer Regelstrecke, deren
%                   Übertragungsfunktion (Verstärkung und Zeitkonstanten)
%                   bekannt ist.
% Eingabeparameter: - Gs :  Übertragungsfunktion der Strecke als Array
%                   - w:    Kreisfrequenzspektrum in Listenform
%                   - phir: gewünschter Phasenrand (->Überschwingen)
%                   - k:    Verstärkung der Strecke
%                   - T:    Zeitkonstanten der Strecke als Array
%
% Ausgabeparameter: - kR:   Reglerverstärkung
%                   - Tn:   Nachstellzeit des Reglers
%                   - Tv:   Vorstellzeit des Reglers
%                   - Tp:   Parasitäre Zeitkonst. des Reglers
% Autor:            Michael Bos, Pascal Ackermann
% Datum:            28.04.2015

function [kR, Tn, Tv, Tp] = PidRegler( Gs,w,phir,kS,T )


% Berechnung von wpid -----------------------------------------------------
phi_s=angle(Gs);                        % Phase der Strecke
   
for y=1:1:length(phi_s)             % Entfernt Sprung bei -pi
    if phi_s(y) > 0                 %--> Punktsuche funktioniert sonst nicht
        for z=y:1:length(phi_s)
            phi_s(z) = -2*pi+phi_s(z); 
        end
        break
    end
end
    
% Liefert die Indices der Liste, wo sich phi= -135 befindet
[ind_left,ind_right] = int_ver(phi_s,-3*pi/4);          %Berechnet die beiden Indizien die ges. wpid einschliessen
wpid=(w(ind_left)+w(ind_right))/2;                      %Weil WpiD irgendwo zwischen W[Left] und W[Right] liegt nehmen wir den arithm Mittelwert

%Steigung = dphis/dw
phi_s_m = (phi_s(ind_right)-phi_s(ind_left))/(w(ind_right)-w(ind_left)); 

% Berechnung Beta und Tnk/Tvk ---------------------------------------------
     
 Ko = -0.5 - (wpid*phi_s_m); 
        
        %Fallunterscheidung Ko
        if(Ko>1)                     
              beta = 1;
        else  % Auflösen der quad. Gleichung nach b (0.5+ wpid*phi_s_m +(2*b/(1+b^2)=0)
              %syms b
              %beta = solve(0.5+ wpid*phi_s_m +(2*b/(1+b^2)),b);
              diskr =((4-4*(0.5+wpid*phi_s_m)*(0.5+wpid*phi_s_m)))^(1/2);
              beta1= ( (-2)+ diskr)/(2*(0.5+wpid*phi_s_m));
              beta2= ( (-2)- diskr)/(2*(0.5+wpid*phi_s_m));
              beta= min(beta1,beta2);       
        end

 Tnk = 1/(wpid*beta);
 Tvk = 1/(wpid/beta);
 
 %Übertragungsfunktion des prov Reglers(K=1) ------------------------------
 Krk=1;                  % provisorische Verstärkung
 Tp = Tvk/10;            % Faustformel für parasitäre Zeitkonstante
        
% Grp = Krk.*(1+( 1./ (1j.*w.*Tnk) ) ).*( 1+1j.*w.*Tvk ).*( 1./ ( 1+1j.*w.*Tp ) ); % Bodekonf 
  Grp = Krk.*((1+1j.*w.*Tnk)./(1j.*w.*Tnk)).*((1+1j.*w.*Tvk)./(1+1j.*w.*Tp));
 % Übertragungsfunktion des offenen Regelkreises
 GOff = Grp.*Gs; 
 
 %Durchtrittspunkt wD bestimmen -------------------------------------------
  
 phi_Go = angle(GOff);              % Phasengang des offenen Regelkreises
 
 for y=1:1:length(phi_Go)             % Entfernt Sprung bei -pi
    if phi_Go(y) > 0                 %--> Punktsuche funktioniert sonst nicht
        for z=y:1:length(phi_Go)
            phi_Go(z) = -2*pi+phi_Go(z); 
        end
        break
    end
 end
 
 indWd = int_ver(phi_Go,-pi+phir);  % Listenindex, wo sich wD befindet
 wD    = w(indWd);
 
 % Krk bestimmen-----------------------------------------------------------
 % G(s) Strecke mit Wd
     % Je nach Ordnung wird die Übertragungsfunktion der Strecke gebildet
    Gs=1;                               %Initialisiere Gs
    for y=1:1:length(T)                 %Berechne Übertragungsfunktion    
        Gs = Gs.*(1./(1+1j.*wD.*T(y)));
    end
    G_str_wd=kS*Gs;                      %Vervollständige mit Verstärkung 

%G(s) des Regler bei wD
G_reg_wd = Krk * ((1+1j*wD*Tnk)/(1j*wD*Tnk))*((1+1j*wD*Tvk)/(1+1j*wD*Tp));

 % Offener Regelkreis bei wD
 GOffwd=G_reg_wd*G_str_wd;        

 % Reglerverstärkung in DB
 amplOffwd=20*log10(abs(GOffwd));
 KrkdB=-amplOffwd;
 % Reglerverstärkung linear
 Krk=10^(KrkdB/20);                    
 
 % Umrechnung in Reglerkonform
 Tv = (Tnk*Tvk)/(Tnk+Tvk-Tp)-Tp;
 Tn = (Tnk+Tvk-Tp);
 kR = (Krk*(Tnk+Tvk-Tp))/Tnk;
 
 
end

% Ende File ---------------------------------------------------------------