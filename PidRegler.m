% Beschreibung:     Dieses File berechnet die grundlegenden Parameter zur
%                   Dimensionerung eines PID-Reglers. Diese Dimensionierung
%                   erfolgt aufgrund einer Regelstrecke, deren
%                   Übertragungsfunktion (Verstärkung und Zeitkonstanten)
%                   bekannt ist.
%                   Die Dimensionierung wird zusätzlich beeinflusst durch
%                   den gewünschten Phasenrand.
% Eingabeparameter: - Gs :  Übertragungsfunktion der Strecke in Listenform
%                   - w:    Kreisfrequenzspektrum in Listenform
%                   - phir: gewünschter Phasenrand (->Überschwingen)
%                   - k:    Verstärkung der Strecke
%                   - T:    Zeitkonstanten der Strecke in Arrayform
% Ausgabeparameter:
%
% Autor:            Michael Bos
% Datum:            22.04.2015

function [Gr] = PidRegler( Gs,w,phir,k,T )


% Berechnung von wpid -----------------------------------------------------
phi_s=angle(Gs);                        % Phase der Strecke
   
for m=1:1:length(phi_s)                 % Entfernt Sprung bei -pi
        if phi_s(m) > 0
            phi_s(m) = -2*pi+phi_s(m);  
        end
end
    
% Liefert die Indices der Liste, wo sich phi befindet
[ind_left,ind_right] = int_ver(phi_s,-3*pi/4);          %Berechnet die beiden Indizien die ges. wpid einschliessen
wpid=(w(ind_left)+w(ind_right))/2;                      %Weil WpiD irgendwo zwischen W[Left] und W[Right] liegt nehmen wir den arithm Mittelwert

% % Überprüfen, dass die beiden Omegas nicht gleich sind
% while(w_next_ind == wpid)
%     copy_ind= copy_ind+1;
%     w_next_ind = w(copy_ind);
% end
    
phi_s_m = (phi_s(ind_right)-phi_s(ind_left))/(w(ind_right)-w(ind_left)); %Steigung = dphis/dw

% Berechnung Beta und Tnk/Tvk ---------------------------------------------
     
 Ko = -0.5 - (wpid*phi_s_m);
         
        if(Ko>1)                     
              beta = 1;
        else    
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
        
 Grp = Krk.*(1.+( 1./ (1j.*w.*Tnk) ) ).*( 1.+1j.*w.*Tvk ).*( 1./ ( 1.+1j.*w.*Tp ) ); % Bodekonf 
 
 % Übertragungsfunktion des provisorischen offenen Regelkreises
 GOff = Grp.*Gs; 
 
 %Durchtrittspunkt wD bestimmen -------------------------------------------
  
 phi_Go = angle(GOff);              % Phasengang des offenen Regelkreises
 
 indWd = int_ver(phi_Go,-pi+phir);  % Listenindex, wo sich wD befindet
 wD    = w(indWd);
 
 % Krk bestimmen-----------------------------------------------------------
 % G(s) Strecke mit Wd
     % Je nach Ordnung wird die Übertragungsfunktion der Strecke gebildet
    Gs=1;                                %Initialisiere Gs
    for y=1:1:length(T)                 %Berechne Übertragungsfunktion    
        Gs = Gs.*(1./(1+1j.*wD.*T(y)));
    end
    G_str_wd=k*Gs;                            %Vervollständige mit Verstärkung 
 
 %G_str_wd = k/(((1+1j*wD*T(1))*(1+1j*wD*T(2))*(1+1j*wD*T(3))*(1+1j*wD*T(4))*(1+1j*wD*T(5))*(1+1j*wD*T(6))*(1+1j*wD*T(7))*(1+1j*wD*T(8)))); 
 
 % G(s) Regler mit Wd
 G_reg_wd = Krk*(1+(1/1j*wD*Tnk))*(1+1j*wD*Tvk)*(1/(1+1j*wD*Tp));    
 
 % Offener Regelkreis 
 GOffwd=G_reg_wd*G_str_wd;        

 % Reglerverstärkung in DB
 amplOffwd=20*log10(abs(GOffwd));
 KrkdB=-amplOffwd;                                             
 Krk=10^(KrkdB/20);                    
 Gr= Krk.*(1+(1./1j*w.*Tnk)).*(1+1j*w.*Tvk).*(1./(1+1j*w.*Tp));
 
 % Umrechnung in Reglerkonform
 Tv = (Tnk*Tvk)/(Tnk+Tvk-Tp)-Tp;
 Tn = (Tnk+Tvk-Tp);
 Kr = (Krk*(Tnk+Tvk-Tp))/Tnk;
 
 
end

% Ende File ---------------------------------------------------------------