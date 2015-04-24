%**************************************************************************
% Projekt:      P2 - Java Applikation - Reglerdimensionierung
%               (Phasengangmethode)
% Autor:        Ackermann Pascal
% Beginndatum:  08.04.2015
% Enddatum:     **.**.2015
% Version:      6.0
%**************************************************************************
% Das folgende m-File führt eine komplette Reglerdimensonierung nach 
% Professor Zellweger's Phasengangmethode durch. Als Eingangsparameter
% werden streckenspezifische Zeitkonstanten mitgegeben, mit denen
% schliesslich die Übertragungsfunktion des zur Strecke passenden,
% dimensionierten Reglers berechnet wird.
% (Diese Berechnung findet Anwendung in der Java-Applikation für das 
% Projekt 2 (FS 2015), sodass dieses m-File schliesslich in Java 
% implementiert wird) 
%**************************************************************************
%
% Benutzte Funktionen:  p2_sani.m
%                       
% Eigene Funktionen:    PiRegler.m, ampl_S.m, phase_S.m
%
% Eingabeparameter:     Verzugszeit Tu, Anstiegszeit Tg, Verstärkung k
%                       Reglertyp typ, Phasenrand phir
% Rückgabeparameter:    Übertragungsfunktion des Reglers Gr(s)      
%
%History:               Ver. 01: Dimensionierung Pi-Regler bis max. n=8
%                       Ver. 02: Falsche Formel für Pi-Regler angepasst
%                       Ver. 03: Andere Benennung diff=abw, diffmin=abwmin
%                       diffminkontrolle = abwminkontrolle, k=y
%                       Ver. 04: p2_sani.m auf Stand von Herrn Niklaus
%                       rückgängig machen, Übertragungsfunktion der
%                       Strecken dynamisch für n-te Ordnung machen 
%                       Ver. 05:angle() und abs() Funktionen werden durch
%                       eigene ersetzt. Für Strecke: phase_S(), ampl_S()
%                       Für PI-Regler: phase_Pi, ampl_Pi 
%                       Ver. 06: Eigene Funktionen phase_S(), ampl_S(),
%                       phase_Pi(),ampl_Pi() werden durch angle() und abs()
%                       ersetzt. Einbindung des PID-Reglers, 
%                       Punktsuche mittels Intervallschachtelung
%
%
%Test Aufruf:

 %Präsentation: phasengangmethode(2.8,20.56,0.5,1,pi/4) ;
                         
 %Str1(n=2):phasengangmethode(3.08,30.8,0.5,1,pi/4) 
 %Str2(n=3):phasengangmethode(0.95,9,2,1,pi/4) 
 %Str3(n=3):phasengangmethode(1.4e-3,7.7e-3,1,1,pi/4) 
 %Str4(n=4):phasengangmethode(3.45e-6,15.5e-6,5,1,pi/4) 
 %Str5(n=4):phasengangmethode(64.8e-3,245.2e-3,0.5,1,pi/4) 
 %Str6(n=4):phasengangmethode(8.6,29.6,1,1,pi/4)                                                
 %Str7(n=5):phasengangmethode(16.6,41.7,1,1,pi/4)
 
%**************************************************************************


function[Gr] = phasengangmethode(Tu,Tg,k,typ,phir)

    % Identifikation ------------------------------------------------------
    
    [n,T]=p2_sani(Tu,Tg);                    %Ordnung und Zeitkonstante  
    Tmax=max(T(1:n));                        %ermittle max von T1 bis Tn
    Tmin=min(T(1:n));                        %ermittle min von T1 bis Tn
    wmin= 1/(Tmax*10);
    wmax= 1/(Tmin/10);
    w = logspace(log10(wmin),log10(wmax),1000);
    
      
    % Je nach Ordnung wird die Übertragungsfunktion der Strecke gebildet
    Gs=1;                                %Initialisiere Gs
    for y=1:1:length(T)                 %Berechne Übertragungsfunktion    
        Gs = Gs.*(1./(1+1j.*w.*T(y)));
    end
    Gs=k*Gs;                            %Vervollständige mit Verstärkung 
    
    % PI Regler ************************************************************
   if (typ==1)
       [Gr] = PiRegler( Gs,w,phir,k,T );
   end
    % PID Regler ***********************************************************
   if (typ==2)
       [Gr] = PidRegler(Gs,w,phir,k,T);
   end
    
end



