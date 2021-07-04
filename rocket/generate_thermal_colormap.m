function cmap = generate_thermal_colormap(param, min_param, max_param, cyclic)
% Written by Caitlin Cash '18
% Generates a thermal colormap for coloring particles by some parameter.
%   It colors everything less than the min_param as the "minimum" color and
%   everything more than the max_param as the "maximum" color. 
%
% INPUTS:
%   param: A list of the parameter values that you are making the cmap from
%   min_param: The floor, for colormapping purposes. Having this too low
%       can make seeing differences on a picture difficult. This is often
%       set to the minimum recorded parameter value.
%   max_param: The ceiling, for colormapping purposes. Having this too high
%       can make seeing differences on a picture difficult. This is often
%       set to mean + 2 standard deviations.
%   cyclic: a boolean that says whether or not we want the colormap to be
%       cyclic, meaning the maximum value is the same color as the minimum
%       value. This is useful for phase measurements. 
%
%


%cyclic says whether or not to make a cyclic colormap
if nargin < 4
    cyclic = false;
end

if cyclic
    spectrum = hsv(100);
else
    spectrum = hsv(140);
    spectrum = spectrum(1:100,:);
end

cmap = zeros(length(param),3);
for index = 1:length(param)
    if param(index) > max_param
        spectrum_entry = 100;
    elseif param(index) < min_param
        spectrum_entry = 0;
    else
        spectrum_entry = ceil(100*(abs(param(index) - min_param)/(max_param - min_param)));
    end
    spectrum_entry = 100 - spectrum_entry;
    if (spectrum_entry == 0) || (spectrum_entry == -1)
        spectrum_entry = 1;
    end
    
    if ~isnan(spectrum_entry)
        val = spectrum(spectrum_entry,:);
    else
        val = 0;
    end
    cmap(index,:) = val;
end


end

