function d = proc_naninterpolation(s, typ)
% d = proc_naninterpolation(s, [typ])
%
% The function interpolates the input signal s over the NaN value. It returns a signal
% d with the same size of s. Signal are in the format [samples x channels].
% Optional: type of interpolation ['linear']. See interp1 for other
% interpolation types.

    if nargin < 2
        typ = 'linear';
    end

    NumChannels = size(s, 2);
    d = zeros(size(s));
    
    for chId = 1:NumChannels
        cs = s(:, chId);
        
        % Interpolate NaN values
        nanId     = isnan(cs);
        
        if isequal(length(unique(nanId)), 1) == false
            tval      = 1:numel(cs);
            cs(nanId)  = interp1(tval(~nanId), cs(~nanId), tval(nanId), typ);
        end
        
        d(:, chId) = cs;
        
    end
    
end