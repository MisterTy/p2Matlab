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
ind = int_ver(phi_s,-3*pi/4);           
wpid = w(ind);

% Nächster Punkt zur Bildung der Steigungstangente
copy_ind = ind+1;
w_next_ind = w(copy_ind);

% Überprüfen, dass die beiden Omegas nicht gleich sind
while(w_next_ind == wpid)
    copy_ind= copy_ind+1;
    w_next_ind = w(copy_ind);
end
    
phi_s_m = (phi_s(copy_ind)-phi_s(ind))/(w_next_ind-wpid);

% -------------------------------------------------------------------------
% Alte Lösungsvariante zur Berechnung von wpid

% differ=abs(phi_s+3*pi/4);                %Suchen Index bei min Abweichung
% diffmin= min(differ);
%              
% indices=find(differ == diffmin); 
% 
%     if length(indices)>1             %Falls zwei Punkte gleichen Abstand
%         indices= indices(2);         % Vermerk für Pascal!!!!!!!!
%                                      % indices(2) nehmen weil unten
%                                      % indices+1
%     end
%   
%  wpid= w(indices);                   % wpid bei -135° 
%  
%  % Nächsten Punkt zur Bildung der Tangente 
%  diffmin_next= differ(indices+1);    
%  w_next_ind = w(indices+1);
%  
%  %Steigung der Tangente bei Wpid
%  phi_s_m = -(diffmin_next-diffmin)/(w_next_ind-wpid);
%--------------------------------------------------------------------------

% Berechnung Beta und Tnk/Tvk ---------------------------------------------
     
 Ko = -0.5 - phi_s_m;
         
        if(Ko>1)                     
              beta = 1;
        else    
              % syms b
              %beta = solve(0.5+ wpid*phi_s_m +(2*b/(1+b^2)),b);
              diskr =((4-4*-(0.5+wpid*phi_s_m*(0.5+wpid*phi_s_m)))^(1/2));
              beta1= ( -(-2)+ diskr)/(2*-(0.5+wpid*phi_s_m));
              beta2= ( -(-2)- diskr)/(2*-(0.5+wpid*phi_s_m));
              beta= max(beta1,beta2);
               
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
 
% -------------------------------------------------------------------------
% Alte Lösungsvariante zur Berechnung von wD
%  % Bei (-180° + Phasenrand) befindet sich Wd, weshalb Punkt von Phase
%  % offener Regelkreis suchen, der diesem Wert am nächsten ist.
%  differWd = abs(phi_Go -(-pi + phir));            
%  diffminWd = min(differWd);    
%  
%  indWd=find(differWd == diffminWd); 
% 
%  %Falls zwei Punkte gleichen Abstand
%     if length(indWd)>1           
%         indWd= indWd(2);         
%     end
%   
%  wD= w(indWd);
%--------------------------------------------------------------------------
 
 % Krk bestimmen-----------------------------------------------------------
 % G(s) Strecke mit Wd
 G_str_wd = k/(((1+1j*wD.*T(1)).*(1+1j*wD.*T(2)).*(1+1j*wD.*T(3)).*(1+1j*wD.*T(4)).*(1+1j*wD.*T(5)).*(1+1j*wD.*T(6)).*(1+1j*wD.*T(7)).*(1+1j*wD.*T(8)))); 
 
 % G(s) Regler mit Wd
 G_reg_wd = Krk*(1+(1/1j*wD*Tnk))*(1+1j*wD*Tvk)*(1/(1+1j*wD*Tp));    
 
 % Offener Regelkreis 
 GOffwd=G_reg_wd*G_str_wd;        

 % Reglerverstärkung in DB
 amplOffwd=20*log10(abs(GOffwd));
 KrkdB=-amplOffwd;                                             
 Krk=10^(KrkdB/20);                    
 Gr= Krk.*(1+(1./1j*w.*Tnk)).*(1+1j*w.*Tvk).*(1./(1+1j*w.*Tp));
 
 
end

% Ende File ---------------------------------------------------------------