%an efficent way of smoothing the low frequnecy microtremor signals
%this code is compiled  from konno_ohmachi.py by Jian Shi 
%erman sentürk&hliva

%  Original paper:
%         K. Konno & T. Ohmachi (1998) "Ground-motion characteristics estimated
%         from spectral ratio between horizontal and vertical components of
%         microtremor." Bulletin of the Seismological Society of America.
%         Vol.88, No.1, 228-241.
%[Inputs]
%         signal: Signal to be smoothed in frequency domain.
%         freq_array: Frequency array corresponding to the signal.
%                     It must have the same length as "raw_signal".
%         smooth_coeff: A parameter determining the degree of smoothing.
%                       The lower this parameter, the more the signal
%                       is smoothed.
%         
% 
%           (Note: "raw_signal" and "freq_array" can be Python lists, 1D numpy
%                  arrays, or 2D 1-column/1-row arrays. The data type affects the
%                  running time. For optimum speed, use 1D numpy array as input.)
%     [Output]
%         y: Smoothed signal (1D numpy array).

% The formula of Konno-Ohmachi smoothing window is here:
%         http://www.geopsy.org/wiki/index.php/Smoothing_details

function y = kohmachi(signal,freq_array,smooth_coeff)
    
%     x = signal;
%     f = freq_array;
%     %f_shifted = f/(1+1e-4);
%     y = zeros(length(x),1);
%     
%     for i = 1 : length(x)
%      %   if (i ~= 1) && (i ~= L)
%             w(:,i) = (sin(log10(f/f(i))*smooth_coeff)./(log10(f/f(i))*smooth_coeff)).^4;
%             %z = f_shifted / f(i);
%             %w = ((sin(smooth_coeff * log10(z)) / smooth_coeff) ./ log10(z)) .^ 4;
%             w(isnan(w)) = 0;
%             y(i) = (w(:,i)' * x) / sum(w(:,i));
%       %  end
%       %  clearvars fc
%     end
%     %y(1) = y(2);
%     %y(L) = y(L-1);
    
    ii = (1:length(freq_array))';
%     [~, b] = size(signal);
%     if b==length(freq_array)
%         signal = signal';
%     end
    repmat = repelem(freq_array,1,length(freq_array));
    ratio = repmat(:,ii)./freq_array(ii)';
    clear repmat
    rat = log10(ratio);
    clear ratio
    w = (sin(rat*smooth_coeff)./(rat*smooth_coeff)).^4;
    clear rat
    w(isnan(w)) = 0;
    s = sum(w(:,ii))';
    avg = w'./s;
    clear w s
    [~, b] = size(signal);
    jj = (1:b)';
    zz = (1:length(freq_array))';
    y(zz,jj) = avg(zz,ii)*signal(ii,jj);
    %y = (w'./ s)*signal;    
    clear ii jj avg signal
end