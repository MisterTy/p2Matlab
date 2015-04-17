function [Gr] = PiRegler( Gs,w,phir,k,T )
%PIREGLER Summary of this function goes here
%   Detailed explanation goes here

    %Berechnung von wpi/Tn
    phi_s=phase_S(w,T);                 %Phasengang Strecke
    for y=1:1:length(phi_s)             % Entfernt Sprung bei -pi
        if phi_s(y) > 0
            phi_s(y) = -2*pi+phi_s(y); 
        end
    end                                     
    abw=abs(phi_s+pi/2);                %Suchen Index bei min Abweichung
    abwmin= min(abw);
    indices=find(abw == abwmin); 
    if length(indices)>1                 %Falls zwei Punkte gleichen Abstand
        indices= indices(1);
    end
    
    wpi= w(indices);
    Tn=1/wpi
%     %% Kontrolle mit berechnetem wpi und angle()
%     Gskontrolle=k/((1+1j*wpi*T(1))*(1+1j*wpi*T(2))*(1+1j*wpi*T(3))*(1+1j*wpi*T(4))*(1+1j*wpi*T(5))*(1+1j*wpi*T(6))*(1+1j*wpi*T(7))*(1+1j*wpi*T(8)));
%     phikontrolle= angle(Gskontrolle)
%     abwminkontrolle = phikontrolle+pi/2
%     abwmin
         
    %Durchtrittspunkt wD bestimmen
    phi_rprov= phase_Pi(w,Tn);          %Phasengang des prov. Pi-Reglers mit Kr=1
    phi_O = phi_s+phi_rprov;            %Phasengang des offenen Regelkreises            
    for y=1:1:length(phi_O)             % Entfernt Sprung bei -pi
        if phi_O(y) > 0
            phi_O(y) = -2*pi+phi_O(y); 
        end
    end                                     
    abw=abs(phi_O+(+pi-phir));          %Suchen Index bei min Abweichung
    abwmin= min(abw);
    indices=find(abw == abwmin); 
    if length(indices)>1                 %Falls zwei Punkte gleichen Abstand
        indices= indices(1);
    end
    wD = w(indices);
%     %% Kontrolle mit berechnetem wD und angle()
%     Grpkontrolle = 1*(1+1/(j*wD*Tn));
%     Gskontrolle=k/((1+1j*wD*T(1))*(1+1j*wD*T(2))*(1+1j*wD*T(3))*(1+1j*wD*T(4))*(1+1j*wD*T(5))*(1+1j*wD*T(6))*(1+1j*wD*T(7))*(1+1j*wD*T(8)));
%     GOkontrolle = Grpkontrolle*Gskontrolle;
%     phiOkontrolle= angle(GOkontrolle)
%     abwminkontrolle = phiOkontrolle+(+pi-phir)
%     abwmin
    
    % KR bestimmen
    ampl_s=ampl_S(w,T,k);               %Amplitudengang Strecke
    ampl_rprov= ampl_Pi(w,Tn,1);        %Amplitudengang des prov. Pi-Reglers mit Kr=1
    ampl_O = ampl_s.*ampl_rprov;        %Amplitudengangs des offenen Regelkreises
    amplOwd = 20*log10(ampl_O(indices));% Amplitude bei wD
    KrdB=-amplOwd;                       %Reglerverstärkung in DB
    Kr=10^(KrdB/20)                      %Reglerverstärkung in DB
    Gr=Kr*(1+1./(1j.*w.*Tn));            %Übertragungsfunktion Regler
    
     %% Plot Frequenzgang des offenen Regelkreises
    phi_Op=phi_O/pi*180;            %Umrechnung in Grad für Plot
    ampl_Op=20*log10(ampl_O);       %Umrechnung in dB für Plot
    figure(4)
    subplot(211),semilogx(w,ampl_Op),xlabel('w/log'),ylabel('Ampl[dB]');
    title('Frequenzgang des offenen Regelkreises')
    subplot(212),semilogx(w,phi_Op),xlabel('w/log'),ylabel('Phase[°]');
    
end

