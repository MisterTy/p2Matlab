%**************************************************************************
% Projekt:      P2 - Java Applikation - Reglerdimensionierung
%               (Phasengangmethode)
% Autor:        Ackermann Pascal, Michael Bos
% Beginndatum:  08.04.2015
% Enddatum:     28.04.2015
% Version:      8.0
%**************************************************************************
% Das folgende m-File führt eine komplette Reglerdimensonierung nach 
% Professor Zellweger's Phasengangmethode durch. Als Eingangsparameter
% werden streckenspezifische Zeitkonstanten(tu,tg) und die Verstärkung resp. 
% Dämpfung mitgegeben(Ks).Mit diesen wird schliesslich die 
% Übertragungsfunktion des zur Strecke passenden Reglers berechnet.

%**************************************************************************
%
% Benutzte Funktionen:  p2_sani.m
%
% Benutzte Ressourcen:  p2_sani_tu_tg.mat(enthält Sanikurven bis n=8) 
%                       
% Eigene Funktionen:    PiRegler.m, PidRegler.m, (int_ver.m)
%
% Eingabeparameter:     Verzugszeit Tu, Anstiegszeit Tg, Verstärkung kS
%                       Reglertyp typ(String), Phasenrand phir
% Rückgabeparameter:    kR, Tn, Tv, Tp   
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
%                       Ver. 07: Anpassung Intervallverschachtelung mit 
%                       Abbruchbedingung Rechterrand - Linkerrand = 1,
%                       Kommentare vervollständigt, 
%                       Sprung bei Phasengang richtig behoben(war zuvor buggy)
%                       Rückgabeparameter ist nicht mehr Gr(s) sondern
%                       Tn,Tv,kR,Tp
%                       Ver. 08: Einbinden der Schrittantwort. 
%                       Berechnung mit step() sowie FFT



%Test Aufruf:
                         
 %Str1(n=2):[kR, Tn, Tv, Tp]=phasengangmethode_sa(3.08,30.8,0.5,'Pi',pi/4) 
 %Str2(n=3):[kR, Tn, Tv, Tp]=phasengangmethode_sa(0.95,9,2,'Pi',pi/4) 
 %Str3(n=3):[kR, Tn, Tv, Tp]=phasengangmethode_sa(1.4e-3,7.7e-3,1,'Pi',pi/4) 
 %Str4(n=4):[kR, Tn, Tv, Tp]=phasengangmethode_sa(3.45e-6,15.5e-6,5,'Pi',pi/4) 
 %Str5(n=4):[kR, Tn, Tv, Tp]=phasengangmethode_sa(64.8e-3,245.2e-3,0.5,'Pi',pi/4) 
 %Str6(n=4):[kR, Tn, Tv, Tp]=phasengangmethode_sa(8.6,29.6,1,'Pi',pi/4)                                                
 %Str7(n=5):[kR, Tn, Tv, Tp]=phasengangmethode_sa(16.6,41.7,1,'Pi',pi/4)
 
 %Str1(n=2):[kR, Tn, Tv, Tp]=phasengangmethode_sa(3.08,30.8,0.5,'Pid',pi/4) 
 %Str2(n=3):[kR, Tn, Tv, Tp]=phasengangmethode_sa(0.95,9,2,'Pid',pi/4) 
 %Str3(n=3):[kR, Tn, Tv, Tp]=phasengangmethode_sa(1.4e-3,7.7e-3,1,'Pid',pi/4) 
 %Str4(n=4):[kR, Tn, Tv, Tp]=phasengangmethode_sa(3.45e-6,15.5e-6,5,'Pid',pi/4) 
 %Str5(n=4):[kR, Tn, Tv, Tp]=phasengangmethode_sa(64.8e-3,245.2e-3,0.5,'Pid',pi/4) 
 %Str6(n=4):[kR, Tn, Tv, Tp]=phasengangmethode_sa(8.6,29.6,1,'Pid',pi/4)                                                
 %Str7(n=5):[kR, Tn, Tv, Tp]=phasengangmethode_sa(16.6,41.7,1,'Pid',pi/4)
%**************************************************************************


function[kR, Tn, Tv, Tp, kS, T] = phasengangmethode_sa(Tu,Tg,kS,typ,phir)

    % Identifikation der Strecke-------------------------------------------
    [n,T]=p2_sani(Tu,Tg);           %n = Ordnung, T=Zeitkonstanten 
    Tmax=max(T(1:n));               %ermittle max von T1 bis Tn
    Tmin=min(T(1:n));               %ermittle min von T1 bis Tn
    wmin= 1/(Tmax*10);              
    wmax= 1/(Tmin/10);
    w = logspace(log10(wmin),log10(wmax),2^12);   %w im relevanten Bereich
    
      
    % Übertragungsfunktion der Strecke berechnen
    Gs=1;                               %Initialisiere Gs
    for y=1:1:length(T)                 %ohne kS
        Gs = Gs.*(1./(1+1j.*w.*T(y)));
    end
    Gs=kS*Gs;                           %Vervollständige mit kS
    
    % PI-Regler ***********************************************************
   if (strcmp(typ,'Pi'))
       [kR, Tn, Tv, Tp] = PiRegler( Gs,w,phir,kS,T );
       schrittantwort( kR, Tn, Tv, Tp, kS, T, w, 'Pi');
   end
    % PID-Regler **********************************************************
   if (strcmp(typ,'Pid'))
       [kR, Tn, Tv, Tp] = PidRegler(Gs,w,phir,kS,T);
       schrittantwort( kR, Tn, Tv, Tp, kS, T, w, 'Pid');
   end
end



