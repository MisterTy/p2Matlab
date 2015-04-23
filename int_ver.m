% Beschreibung:      Dieses File erhält eine Liste mit komplexen Werten,
%                    welche den Phasengang einer Strecke bilden. Aus diesem
%                    Phasengang wird der gewünschte Winkel "phi" mit der
%                    Intervallverschachtelung gesucht. Der Index, bei
%                    welchem sich "phi" befindet, entspricht auch dem Index
%                    der Omega-Liste, wo sich wpi/wpid befindet.
%
% Übergabeparameter: - Phasengang in Form einer Liste
%                    - Gesuchte Phase 
% Rückgabeparameter: - Index der Liste
%
% Autor:             Michael Bos
% Datum:             22.04.2015
%**************************************************************************
function [ind] =  int_ver(phase,phi)

right = length(phase);         % Länge der Liste ist der rechte Rand  
left = 1;                      % Linke Rand beginnt bei Indices 1

for n=1:1:15
   
    ind=floor((right+left)/2);  % Hälfte der Anzahl Indices
    curr_phase = phase(ind);    % Phase bei der Hälfte
    
    if(curr_phase>phi)
        left = (right+left)/2;   % Linker Rand auf gemittelten Wert setzen
        right = right;
    else
        right = (left+right)/2;  % Rechter Rand auf gemittelten Wert setzen  
        left = left;
    end
end
end

% File ENDE ---------------------------------------------------------------