function [dtoptim] = vs_info_pick_dt(idt, dt)

dt_error = [];
idt_error = [];

for i = 1:length(dt)-25

   index = [i:length(dt)];
   dtsubset = dt(index);
   idtsubset = idt(index);

   beta = polyfit(dtsubset, idtsubset, 1);
   idtsubset_fit = polyval(beta, dtsubset);
   idtfit = polyval(beta, dt);

   e = abs( (idtsubset - idtsubset_fit) ./ idtsubset );
   idt_error = [idt_error max(e)];
   dt_error = [dt_error min(dtsubset)];

%    clf;
%    hold on;
%    plot(dt, idt, 'ko', 'markersize', 2);
%    plot(dtsubset, idtsubset_fit, 'k-');
%    plot(dt(min(index)), idt(min(index)), 'ko', 'markerfacecolor', 'k');
%    pause;

end % (for i )

index = find(idt_error > 0.2);
dtoptim = dt_error(max(index));

return;



% if ( 0 )
%    % fitting a line on the I vs. dt graph
%    plot([0 dt(51) dt(100)],[intercepts(k) I_all(51) I_all(100)]);
% 
%    xlabel('dt(seconds)');
%    ylabel('I(dt)');
%    pause;
%    close
% end


% slope = (idt(100)-idt(51))/(dt(100)-dt(51));
% intercept = idt(100)-slope*dt(100);
% 
% I_predicted = slope.*dt+intercept;
% error = abs((idt-I_predicted)./idt);
% dt_index = min(find(error <= 0.01));
% dt_pick = dt(dt_index);
% dtoptim = dt_pick;




