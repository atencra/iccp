function params = iccp_strf_parameters(strf, trigger)
%iccp_strf_parameters - STRF parameters from significant strf
%
% Calculates temporal, spectral modulation transfer functions as well 
% as the separability index for all strfs.
%
% params = iccp_strf_parameters(strf, trigger)
%
% strf : struct array of spectrotemporal receptive fields
%
% trigger : ripple stimulus triggers, in sample number.
%
% The output 'params' is a struct array having length = length(strf)
%
% params has the following form:
% 
%    params.exp -> experiment
%    params.site -> recording site
%    params.chan -> channel of neuron
%    params.model -> sorted model number
%    params.depth -> probe depth
%    params.position -> neuron depth
%    params.stim -> dmr1, dmr2, ripple noise, etc
%    params.atten -> attenuation of sound
%    params.sm -> max spectral modulation in stimulus
%    params.tm -> max temporal modulation in stimulus
%    params.mdb -> modulation depth of stimulus
%    params.spl -> sound level of stimulus
%    params.n0 -> number of spikes 
%    params.w0 -> firing rate
%    params.percent_energy -> different power levels for STRF analysis
%    params.rfenergy -> energy in STRF = root mean square (RMS)
%    params.tmf -> temporal modulation frequency axis
%    params.xmf -> spectral modulation frequency axis
%    params.rtf -> ripple transfer function
%    params.singvals -> singular vals from STRF decomposition
%    params.eigvals -> eigenvalues for decomposition
%    params.tci -> temporal correlation index
%    params.sci -> spectral correlation index
%    params.pli -> phase-locking index
% 


params = struct(...
'exp',            [], ...
'site',           [], ...
'chan',           [], ...
'model',          [], ...
'depth',          [], ...
'position',       [], ...
'stim',           [], ...
'atten',          [], ...
'sm',             [], ...
'tm',             [], ...
'mdb',            [], ...
'spl',            [], ...
'n0',             [], ...
'w0',             [], ...
'percent_energy', [], ...
'rfenergy',       [], ...
'tmf',            [], ...
'xmf',            [], ...
'rtf',            [], ...
'singvals',       [], ...
'eigvals',        [], ...
'tci',            [], ...
'sci',            [], ...
'pli',            []);


percent_energy = [80 85 90 95 97.5 100];
energy = zeros(size(percent_energy));


for i = 1:length(strf)
i
strf(i).chan
strf(i).model

    % Define some initial parameters
    p = 0.001;
    n0 = strf(i).n0contra;
    w0 = strf(i).w0contra;
    sm = strf(i).sm;
    tm = strf(i).tm;
    mdb = strf(i).mdb;
    stim = strf(i).stim;
    fs = strf(i).fs;
    t = strf(i).taxis;
    x = log2(strf(i).faxis ./ min(strf(i).faxis));
    dur = (trigger(end)-trigger(1)) / fs; % duration in seconds
    
    % Get the significant strf
    [rfsig] = significant_strf(strf(i).rfcontra, p, n0, mdb, dur);
   
    % take the singular value decomposition to decompose
    % the strf into separable subunits
    [u,s,v] = svd(rfsig);
    singvals = sum(s);
    eigvals = singvals .^ 2;
    
    % Get phase-locking index
    pli = phase_locking_index(rfsig, mdb, w0, stim);
    
    % Get temporal correlation index
    deltat = t(2)-t(1);
    [rt, tshift] = temporal_correlation_index(rfsig, deltat);
    
    % Get spectral correlation index
    deltax = x(2)-x(1);
    [rx, xshift] = spectral_correlation_index(rfsig, deltax);
    
    % get strf energy at different power levels
    cumpower = 100 * cumsum(eigvals) / sum(eigvals+eps);
    
    rfsvd = zeros(size(rfsig,1), size(rfsig,2), length(percent_energy));
    
    for j = 1:length(percent_energy)
        
        % Get rid of noise for plotting appearance
        index_percent_power = find(cumpower <= percent_energy(j));
        num_sing_vals = min([max(index_percent_power) size(u,2) length(singvals)]);
        
        rftemp = zeros(size(rfsig));
        for k = 1:num_sing_vals
            rftemp = rftemp + singvals(k) * u(:,k) * v(:,k)';
        end % (for k)

        % Get STRF energy
        [energy(j)] = strf_energy(rftemp, mdb, stim);
        
        % Save RFs at various SVD power levels
        rfsvd(:,:,j) = rftemp;
        
    end % (for j)
    
    
    dt = diff(t);
    dt = dt(1); % strf temporal resolution
    dx = diff(x);
    dx = dx(1); % strf spectral resolution
    
    maxtmf = ceil( 1 / dt );
    maxsmf = ceil( 1 / dx );
    
    ntbins = 161; % this will give us 12.2 Hz tmtf resolution % generalize
    nfbins = ceil(maxsmf / 0.15); % make resolution 0.15 cycles / octave
    
    dtfreq = maxtmf / ntbins; %size(rtftemp,2); % fft temporal frequency resolution
    dxfreq = maxsmf / nfbins; %size(rtftemp,1); % fft spectral frequency resolution
    
    
    for k = 1:size(rfsvd,3)
        
        % get the ripple transfer function/s
        rfk = fft2(rfsvd(:,:,k), nfbins, ntbins);
        rtftemp = fftshift(abs(rfk));
        
        % Get the tm frequency vector - will be in cycles per second
        if ( mod(size(rtftemp,2),2) )
            tmf = [-(size(rtftemp,2)-1)/2:(size(rtftemp,2)-1)/2]*dtfreq;
        else
            tmf = [-size(rtftemp,2)/2:(size(rtftemp,2)/2-1)]*dtfreq;
        end
        
        itmf0 = find(tmf==0);
        itmfp40 = find(tmf<=tm);
        ditmf = max(itmfp40)-itmf0;
        itmf = [itmf0-ditmf:itmf0+ditmf];
        tmf = tmf(itmf);
        
        % Get the sm frequency vector - will be in cycles per octave
        if ( mod(size(rtftemp,1),2) )
            xmf = [-(size(rtftemp,1)-1)/2:(size(rtftemp,1)-1)/2]*dxfreq;
        else
            xmf = [-size(rtftemp,1)/2:(size(rtftemp,1)/2-1)]*dxfreq;
        end
        
        ixmf0 = find(xmf == 0);
        ixmf4 = find(xmf <= 4);
        ixmf = ixmf0:max(ixmf4);
        xmf = xmf(ixmf);
        
        rtf(:,:,k) = rtftemp(ixmf, itmf);
        
    end % (for k)
    
    
    % assign parameters to output struct array
    params(i).exp = strf(i).exp;
    params(i).site = strf(i).site;
    params(i).chan = strf(i).chan;
    params(i).model = strf(i).model;
    params(i).depth = strf(i).depth;
    params(i).position = strf(i).position;
    params(i).stim = strf(i).stim;
    params(i).atten = strf(i).atten;
    params(i).sm = strf(i).sm;
    params(i).tm = strf(i).tm;
    params(i).mdb = strf(i).mdb;
    params(i).spl = strf(i).spl;
    params(i).n0 = n0;
    params(i).w0 = w0;
    params(i).percent_energy = percent_energy(:);
    params(i).rfenergy = energy(:);
    params(i).tmf = tmf;
    params(i).xmf = xmf;
    params(i).rtf = rtf;
    params(i).singvals = singvals;
    params(i).eigvals = eigvals;
    params(i).tci = [tshift(:) rt(:)];
    params(i).sci = [xshift(:) rx(:)];
    params(i).pli = pli;
    
end % (for i)

 
