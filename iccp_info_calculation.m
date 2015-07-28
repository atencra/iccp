function [info_nmse,sigma_nmse,nmse] = vs_info_calculation(resp,graphs,popgraphs)
% NEEDS TO BE MODULARIZED

% info_calculation: Calculates information values and nmse values from an
% STRF response file. (saved as resp_final)

% [info_nmse, sigma_nmse, nmse] = info_calculation(resp, graphs, popgraphs)

% Inputs:
% -----------------------------
% 1) resp: STRF response file, saved as -*final-resp. 
% 2) graphs: boolean input. 1 if graphs are desired, 0 if not.
% 3) popgraphs: boolean input. 1 if population graphs are desired, 0 if
%       not.

% Outputs:
%   ----------------------------
%   IF GRAPHS is 1 (plots are for each individual neuron):
%
%   1) Generates a plot of Information vs. all the dt values.
%   2) Generates a PSTH, raster plot, and image of the stimulus 
%             spectrogram, using the value of dt where the I vs dt graph 
%             becomes nonlinear.
%   3) Plots all of its information values, divided into intervals of 
%             width [[0.1 0.5 1 2 3 4 6 12], against log(1/those T values).
%   4) Plots sigma against 1/sqrt(those T values).

%   IF POPGRAPHS is 1 (plots are for all the neurons)
%     1) Plots the population histogram of normalized mean square error, 
%             information, and picked dt values.

%  Finally, performs pairwise analysis on information, sigma values, and
%  generates scatterplots of those.

T = 12;

dt = linspace(0.0001,0.020,100);

L = length(resp);
info_pop = [];

info_nmse = cell(1,L);
sigma_nmse = cell(1,L);

nmse_values = zeros(1,L);
info_general = zeros(1,L);
intercepts = zeros(1,L);
dt_pick = zeros(1,L);

for k = 1:L
    
    if graphs
        figure;
    end
    
    I_all = zeros(1,100);
    for n = 1:100 %length(dt)
        
        [r] = vs_plot_psth(resp(k),dt(n),0);
        
        rt = r{1}./150./dt(n);
        rbar = sum(r{1})./150./T;
        
        R = length(rt);
        
        I = 0;
        
        for m = 1:R
            if rt(m) == 0
                I = I;
            else
                I = I + 1/T.*(rt(m)./rbar).*log2(rt(m)./rbar).*dt(n);
            end
        end
        
        I_all(n) = I;
        
        if graphs
            hold on
            plot(dt(n),I,'ko') % I vs. dt
            xlabel('dt (seconds)')
            ylabel('Information')
            title('I vs. dt')
        end
        box off
        tickpref
    end
    info_pop = [info_pop I_all];
    
    slope = (I_all(100)-I_all(51))/(dt(100)-dt(51));
    intercepts(k) = I_all(100)-slope*dt(100);
    
    I_predicted = slope.*dt+intercepts(k);
    error = abs((I_all-I_predicted)./I_all);
    dt_index = min(find(error <= 0.01));
    dt_pick(k) = dt(dt_index);
    
    if graphs
        % fitting a line on the I vs. dt graph
        plot([0 dt(51) dt(100)],[intercepts(k) I_all(51) I_all(100)]);
        
        xlabel('dt(seconds)');
        ylabel('I(dt)');
        pause;
        close
    end
    
    
    % Value of dt used: Value where I vs dt graph becomes nonlinear,
    % individualized per channel
    
    T_values = [0.1 0.5 1 2 3 4 6 12];
    tL = length(T_values);
    
   [R1] = vs_plot_psth(resp(k),dt_pick(k),0);
%     if graphs
%         [R1] = vs_plot_psth(resp(k),dt_pick(k),1);
%     else
%         [R1] = vs_plot_psth(resp(k),dt_pick(k),0);
%     end
    
    RT1 = R1{1}./150./dt_pick(k);
    
    I_values = cell(1,tL);
    I_values2 = cell(1,tL);
    
    if graphs
        figure;
    end
    
    sd = zeros(1,tL);
    sd2 = zeros(1,tL);
    
    for j = 1:tL %length(T_values)
        
        %Divide overall interval into intervals of length T_values(j),
        %extract corresponding R values into arrays representing those
        %intervals
        
        numInt = 12./T_values(j);
        index = floor((length(R1{1})-1)./numInt);
        
        R_T = cell(1,numInt);
        R_trunc = cell(1,numInt);
        RBAR = cell(1,numInt);
        
        I_new = zeros(1,numInt);
        
        for s = 1:numInt
            R_T{s} = RT1((s-1)*index+1:s*index);
            R_trunc{s} = R1{1}((s-1)*index+1:s*index);
            
            RBAR{s} = sum(R_trunc{s})./150./T_values(j);
            
            y = length(R_T{s});
            
            for i = 1:y
                if R_T{s}(i) == 0
                    I_new(s) = I_new(s);
                else
                    I_new(s) = I_new(s) + 1/T_values(j).*(R_T{s}(i)./RBAR{s}).*log2(R_T{s}(i)./RBAR{s}).*dt_pick(k);
                end
            end
            
            
        end
        
        I_values_total = sort(I_new);
        
        I_values{j} = I_values_total;
        
        % Removing outliers
        if j <= 3
            I_values2{j} = I_values_total(2:end-1);
            
        else
            I_values2{j} = I_values_total;
            
        end
        
        sd(j) = std(I_values{j});
        sd2(j) = std(I_values2{j});
        
        T_axis = 1./T_values(j).*ones(1,length(I_values{j}));
        
        T_axis2 = 1./T_values(j).*ones(1,length(I_values2{j}));
        
        if graphs
            subplot(2,1,1)
            
            hold on
            
            plot(log(T_axis),I_values{j},'ko')
            
            xlabel('log(1/T) (s^-1)')
            ylabel('I(dt)')
            
            subplot(2,1,2)
            
            hold on
            
            plot(log(T_axis2),I_values2{j},'bo')
            
        end
        box off
        tickpref
    end
    
    if graphs
        pause;
        close
    end
    
    
    T_axis3 = 1./sqrt(T_values);
    T_axis4 = 1./sqrt(T_values);
    
    % Removing SD Outliers
    index = find(T_axis3 == max(T_axis3));
    T_axis3(index) = [];
    T_axis4(index) = [];
    sd(index) = [];
    sd2(index) = [];
    
    if graphs
        figure;
        subplot(2,1,1)
        plot(T_axis3,sd,'ko')
        
        subplot(2,1,2)
        plot(T_axis4,sd2,'bo')
        
        box off
        tickpref
    end
    
    R = corrcoef(T_axis3,sd);
    R2 = corrcoef(T_axis4,sd2);
    
    r = R(2)^2;
    r2 = R2(2)^2;
    
    P = polyfit(T_axis3,sd,1);
    P2 = polyfit(T_axis4,sd2,1);
    
    sigma_predicted = P(1).*T_axis3 + P(2);
    sigma_predicted2 = P2(1).*T_axis4 + P2(2);
    
    [nmse] = gfit(sd,sigma_predicted,'2');
    [nmse2] = gfit(sd2,sigma_predicted2,'2');
    
    if graphs
        subplot(2,1,1)
        lsline
        xlabel('1/sqrt(T)(s^-.5)')
        ylabel('Sigma')
        title(sprintf('Sigma vs. 1/sqrt(T), nmse = %.3f, r = %.2f, r^2 = %.2f',nmse,R(2),r));
        
        subplot(2,1,2)
        lsline
        xlabel('1/sqrt(T)(s^-.5)')
        ylabel('Sigma')
        title(sprintf('Sigma vs. 1/sqrt(T), nmse = %.3f, r = %.2f, r^2 = %.2f',nmse2,R2(2),r2));
        box off
        tickpref
        
        pause
        close
    end
    
    nmse_values(k) = nmse;
    
    if nmse < 0.15            % we picked this value to find dt
        info_nmse{k} = I_values;
        sigma_nmse{k} = sd;
    end
    info_general(k) = I_values{end};
end

if graphs
    dt_pop = repmat(dt,1,L);
    figure
    subplot(1,[1 2])
    
    scatter(dt_pop,info_pop,'ko')
    
    subplot(2,[1 2])
    scatter(dt_pop,info_pop,'k-')
    
    title('I vs. dt, Population')
    box off
    tickpref
    pause
    close
end
% % Using chosen value of 8 milliseconds to calculate information
% [r] = plot_psth(resp(k),0.008,0);
%
% rt = r{1}./150./0.008;
% rbar = sum(r{1})./150./12;
%
% R = length(rt);
%
% I = 0;
%
% for m = 1:R
%     if rt(m) == 0
%         I = I;
%     else
%         I = I + 1/T.*(rt(m)./rbar).*log2(rt(m)./rbar).*dt(n);
%     end
% end


if popgraphs
    figure
    edges = linspace(0,1,20);
    x = histc(nmse_values, edges);
    x = x./sum(x);
    hb = bar(edges,x,'histc');
    set(hb, 'FaceColor', [0.6 0.6 0.6]);
    xlabel('NMSE')
    ylabel('Proportion')
    title('NMSE Population Data')
    box off
    tickpref
    pause
    close
    
    figure
    edges = linspace(min(info_general),max(info_general),20);
    w = histc(info_general, edges);
    size(info_general)
    w = w./sum(w);
    hb = bar(edges,w,'histc');
    set(hb, 'FaceColor', [0.6 0.6 0.6]);
    xlabel('Information')
    ylabel('Proportion')
    title('Information Population Data')
    box off
    tickpref
    pause
    close
    
    figure
    edges = linspace(min(dt_pick),max(dt_pick),20);
    z = histc(dt_pick, edges);
    z = z./sum(z);
    hb = bar(edges,z,'histc');
    set(hb, 'FaceColor', [0.6 0.6 0.6]);
    xlabel('dt')
    ylabel('Proportion')
    title('Picked dt Population Data')
    box off
    tickpref
    pause
    close
    
end

%Pairwise Analysis
info_nmse
sigma_nmse

info_extracted = zeros(L,172);  % Comes from total number of information values
%                                   (120+24+12+6+4+3+2+1 = 172) that reflects the
%                                   time intervals chosen. See variable
%                                   T_values above.

sigma_extracted = zeros(L,7);    % Comes from the 8 different interval
%                                    lengths picked. There are only sigma
%                                    values because one time interval
%                                    length was the full 12 seconds.
channels = [];
for q = 1:L
    if ~isempty(info_nmse{q})
        info_extracted(q,:) = [info_nmse{q}{1} info_nmse{q}{2} info_nmse{q}{3} info_nmse{q}{4} info_nmse{q}{5} info_nmse{q}{6} info_nmse{q}{7} info_nmse{q}{8}];
        sigma_extracted(q,:) = sigma_nmse{q};
        channels = [channels q];
    end
end

I1 = [];
I2 = [];
S1 = [];
S2 = [];

if length(channels) > 1
    
    cmb = nchoosek(channels, 2); % determine all possible pairwise combinations
    
    [nr, nc] = size(cmb);
    
    for j = 1:nr
        
        index1 = cmb(j,1);
        index2 = cmb(j,2);
        
        I1 = [I1; info_extracted(index1,:)];
        I2 = [I2; info_extracted(index2,:)];
        
        S1 = [S1; sigma_extracted(index1,:)];
        S2 = [S2; sigma_extracted(index2,:)];
    end
end
I1
I2
S1
S2

scatter(I1,I2)
box off
tickpref
scatter(S1,S2)
box off
tickpref
end