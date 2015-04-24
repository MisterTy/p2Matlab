% Beschreibung:      Dieses File erhält eine Liste mit komplexen Werten,
%                    welche den Phasengang einer Strecke bilden. Aus diesem
%                    Phasengang wird der gewünschte Winkel "phi" mit der
%                    Intervallverschachtelung gesucht. Der Index, bei
%                    welchem sich "phi" befindet, entspricht auch dem Index
%                    der Omega-Liste, wo sich wpi/wpid befindet.
%
% Übergabeparameter: - Phasengang in Form einer Liste
%                    - Gesuchte Phase 
% Rückgabeparameter: - Linke und Rechter Rand , welche den gesuchten Index
%                      eingrenzen
%
% Autor:             Michael Bos
% Datum:             22.04.2015
%**************************************************************************
function [left, right] =  int_ver(phase,phi)

right = length(phase);         % Länge der Liste ist der rechte Rand  
left = 1;                      % Linke Rand beginnt bei Indices 1

while(right-left ~=1)           % Iteration bis Linker und Rechter Rand gesuchten Wert eingrenzen
                                % Bemerkung: gesuchter Index liegt zwischen
                                % beiden Ränder
    
    ind=floor((right+left)/2);  % Hälfte der Anzahl Indices
    curr_phase = phase(ind);    % Phase bei der Hälfte
    
    if(curr_phase>phi)
        left = floor(right+left)/2;   % Linker Rand auf gemittelten Wert setzen
        right = right;                % Rechter Rand bleibt
    else
        right = floor(left+right)/2;  % Rechter Rand auf gemittelten Wert setzen  
        left = left;                  % Linker Rand bleibt
    end
end
end

% File ENDE ---------------------------------------------------------------