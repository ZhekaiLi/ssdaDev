function ssda_Data = ssda_graph(ssda_Data,option)

plot_small = 1;
if plot_small
   blueMarker = '.';   blueMarkersize = 10; blueLineWidth = 2; 
   redMarker  = '^';   redMarkersize  = 4;  redLineWidth  = 1; 
   genMarker  = '*';   genMarkersize  = 4;  genLineWidth  = 0.5; 
else
   blueMarker = 'x';   blueMarkersize = 8;  blueLineWidth = 1.5; 
   redMarker  = '^';   redMarkersize  = 8;  redLineWidth  = 2.5; 
   genMarker  = '*';   genMarkersize  = 8;  genLineWidth  = 1.5; 
end

auto_gca = 0;
if ~auto_gca
   XLim = [-5 5];
   YLim = [-5 5];
   XTick = [-5 0 5];
   YTick = [-5 0 5];
end

% plot_thresholds = 0;
% if plot_thresholds
%     "..."
%     u1gth = linspace(XLim(1),XLim(2),1000);
%     for istep = 1:length(data.y)
%        u2gth(istep,:) = "g(u1gth)" - data.y(istep);
%     end
%     u1g0 = u1gth;
%     u2g0 = "g(u1g0)";
% end

switch option
    
   case 0 % Samples lower than threshold value
       
      figid = 2*ssda_Data.Nb_step+1;
       
      if ssda_Data.Nb_step == 0
         ssda_Data.U1Lim = [ 0 0 ];
         ssda_Data.U2Lim = [ 0 0 ];
      end
      
      U1Lim = [ min(ssda_Data.U(1,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2))) max(ssda_Data.U(1,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2))) ];
      U2Lim = [ min(ssda_Data.U(2,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2))) max(ssda_Data.U(2,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2))) ];
      if U1Lim(1) < ssda_Data.U1Lim(1), ssda_Data.U1Lim(1) = floor(U1Lim(1)*2)/2; end  
      if U1Lim(2) > ssda_Data.U1Lim(2), ssda_Data.U1Lim(2) = ceil(U1Lim(2)*2)/2; end
      if U2Lim(1) < ssda_Data.U2Lim(1), ssda_Data.U2Lim(1) = floor(U2Lim(1)*2)/2; end
      if U2Lim(2) > ssda_Data.U2Lim(2), ssda_Data.U2Lim(2) = ceil(U2Lim(2)*2)/2; end
      for i = 1:figid-1
         if auto_gca
            set(gca(i),'XLim',ssda_Data.U1Lim,'YLim',ssda_Data.U2Lim);
         end
         axis square
      end
      
      figure(figid)
      title([ 'Step #' num2str(ssda_Data.Nb_step) ],'FontSize',14)
      hold on
      h1 = plot(ssda_Data.U(1,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2)),ssda_Data.U(2,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2)),blueMarker);
      set(h1,'Markersize',blueMarkersize,'LineWidth',blueLineWidth,'Color','b','MarkerFaceColor','b');
      imax_Indgerm = find(ssda_Data.Indgerm(end,:)==0);
      if isempty(imax_Indgerm)
         imax_Indgerm = size(ssda_Data.Indgerm(end,:),2);
      else
         imax_Indgerm = imax_Indgerm(1)-1;
      end
      h2 = plot(ssda_Data.U(1,ssda_Data.Indgerm(end,1:imax_Indgerm)),ssda_Data.U(2,ssda_Data.Indgerm(end,1:imax_Indgerm)),redMarker);
      set(h2,'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color','r','MarkerFaceColor','r');
%       if plot_thresholds
%          hgth = plot(u1gth,u2gth(end,:),'-');
%          set(hgth,'LineWidth',2,'Color','r')
%          hg0 = plot(u1g0,u2g0,'-');
%          set(hg0,'LineWidth',2,'Color','k')
%       end
      if auto_gca
         set(gca,'FontSize',14,'XLim',ssda_Data.U1Lim,'YLim',ssda_Data.U2Lim);
      else
         set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
         grid on
         box on
      end
      axis square
      
      xlabel('{\itu}_1','FontSize',14);
      ylabel('{\itu}_2','FontSize',14);
      
      disp('Press a key to continue')
      pause
      
    case 1 % MCMC, all states

      figid = 2*ssda_Data.Nb_step;

      figure(figid)
      clf
      title([ 'Step #' num2str(ssda_Data.Nb_step) ],'FontSize',14)
      hold on
      h1 = plot(ssda_Data.U(1,ssda_Data.Indices(end-1,1):ssda_Data.Indices(end-1,2)),ssda_Data.U(2,ssda_Data.Indices(end-1,1):ssda_Data.Indices(end-1,2)),blueMarker);
      set(h1,'Markersize',blueMarkersize,'LineWidth',blueLineWidth,'Color','b','MarkerFaceColor','b');
      h3 = plot(ssda_Data.U(1,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2)),ssda_Data.U(2,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2)),genMarker);
      set(h3,'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','g');
      imax_Indgerm = find(ssda_Data.Indgerm(end,:)==0);
      if isempty(imax_Indgerm)
         imax_Indgerm = size(ssda_Data.Indgerm(end,:),2);
      else
         imax_Indgerm = imax_Indgerm(1)-1;
      end
      h2 = plot(ssda_Data.U(1,ssda_Data.Indgerm(end,1:imax_Indgerm)),ssda_Data.U(2,ssda_Data.Indgerm(end,1:imax_Indgerm)),redMarker);
      set(h2,'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color','r','MarkerFaceColor','r');
%       if plot_thresholds
%          hgth = plot(u1gth,u2gth(end,:),'-');
%          set(hgth,'LineWidth',2,'Color','r')
%          hg0 = plot(u1g0,u2g0,'-');
%          set(hg0,'LineWidth',2,'Color','k')
%       end
      if auto_gca
         set(gca,'FontSize',14,'XLim',ssda_Data.U1Lim,'YLim',ssda_Data.U2Lim);
      else
         set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
         grid on
         box on
      end
      axis square
      
      xlabel('{\itu}_1','FontSize',14);
      ylabel('{\itu}_2','FontSize',14);
      
      disp('Press a key to continue')
      pause

   case 2 % MCMC, one single state
      
      figid = 2*ssda_Data.Nb_step;

      figure(figid)
      hold on
      
      U1Lim = [ min(ssda_Data.Usubtemp(1,:)) max(ssda_Data.Usubtemp(1,:)) ];
      U2Lim = [ min(ssda_Data.Usubtemp(2,:)) max(ssda_Data.Usubtemp(2,:)) ];
      if U1Lim(1) < ssda_Data.U1Lim(1), ssda_Data.U1Lim(1) = floor(U1Lim(1)*2)/2; end  
      if U1Lim(2) > ssda_Data.U1Lim(2), ssda_Data.U1Lim(2) = ceil(U1Lim(2)*2)/2; end
      if U2Lim(1) < ssda_Data.U2Lim(1), ssda_Data.U2Lim(1) = floor(U2Lim(1)*2)/2; end
      if U2Lim(2) > ssda_Data.U2Lim(2), ssda_Data.U2Lim(2) = ceil(U2Lim(2)*2)/2; end
      for i = 1:figid-1
         if auto_gca
             set(gca(i),'XLim',ssda_Data.U1Lim,'YLim',ssda_Data.U2Lim);
         end
         axis square
         xlabel('{\itu}_1','FontSize',14);
         ylabel('{\itu}_2','FontSize',14);
      end
      
      if ssda_Data.Nb_generation == 1
          
         title([ 'Step #' num2str(ssda_Data.Nb_step)-1 ],'FontSize',14);
         % Crude MC samples
         h1 = plot(ssda_Data.U(1,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2)),...
             ssda_Data.U(2,ssda_Data.Indices(end,1):ssda_Data.Indices(end,2)),blueMarker);
         set(h1,'Markersize',blueMarkersize,'LineWidth',blueLineWidth,'Color','b','MarkerFaceColor','b');
         
         imax_Indgerm = find(ssda_Data.Indgerm(end,:)==0);
         if isempty(imax_Indgerm)
            imax_Indgerm = size(ssda_Data.Indgerm(end,:),2);
         else
            imax_Indgerm = imax_Indgerm(1)-1;
         end
         
         % failure domain samples
         h2 = plot(ssda_Data.U(1,ssda_Data.Indgerm(end,1:imax_Indgerm)),...
             ssda_Data.U(2,ssda_Data.Indgerm(end,1:imax_Indgerm)),redMarker);
         set(h2,'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color','r','MarkerFaceColor','r');
         if auto_gca
            set(gca,'FontSize',14,'XLim',ssda_Data.U1Lim,'YLim',ssda_Data.U2Lim);
         else
            set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
            grid on
            box on
         end
         axis square
         
         pause(1)
         
      end
      
      for i = 1:(ssda_Data.Nb_generation-1)
         ht = ssda_Data.ht;
         hthl = ssda_Data.hthl;
         set(ht{i},'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','g');
         if hthl{i} ~= 0
            set(hthl{i},'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','k');
         end
      end

      title([ 'Step #' num2str(ssda_Data.Nb_step) ' - Generation #' num2str(ssda_Data.Nb_generation) ])
      % plot all the generated samples in current generation
      ht{ssda_Data.Nb_generation} = plot(ssda_Data.Usubtemp(1,:),ssda_Data.Usubtemp(2,:),genMarker);
      set(ht{ssda_Data.Nb_generation},'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','y');
      if ~isempty(ssda_Data.ind)
         % plot those samples the safe domain in current generation
         hthl{ssda_Data.Nb_generation} = plot(ssda_Data.Usubtemp(1,ssda_Data.ind),ssda_Data.Usubtemp(2,ssda_Data.ind),genMarker);
         set(hthl{ssda_Data.Nb_generation},'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','m');
      else
         hthl{ssda_Data.Nb_generation} = 0;
      end
%       if plot_thresholds
%          hgth = plot(u1gth,u2gth(end,:),'-');
%          set(hgth,'LineWidth',2,'Color','r')
%          hg0 = plot(u1g0,u2g0,'-');
%          set(hg0,'LineWidth',2,'Color','k')
%       end
      if auto_gca
         set(gca,'FontSize',14,'XLim',ssda_Data.U1Lim,'YLim',ssda_Data.U2Lim);
      else
         set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
         grid on
         box on
      end
      axis square
      
      ssda_Data.ht = ht;
      ssda_Data.hthl = hthl;
      
      xlabel('{\itu}_1','FontSize',14);
      ylabel('{\itu}_2','FontSize',14);
      
      pause(1)

    case 3 % Final plot (all steps)

      figid = 2*ssda_Data.Nb_step+2;
      
      figure(figid)
      title('Final','FontSize',14)
      hold on

      map = colormap(jet(ssda_Data.Nb_step+1));
      lestr = '{';
      for i = 0:ssda_Data.Nb_step
         imax_Indgerm = find(ssda_Data.Indgerm(i+1,:)==0);
         if isempty(imax_Indgerm)
            imax_Indgerm = size(ssda_Data.Indgerm(i+1,:),2);
         else
            imax_Indgerm = imax_Indgerm(1)-1;
         end
         h(i+1) = plot(ssda_Data.U(1,ssda_Data.Indgerm(i+1,1:imax_Indgerm)),ssda_Data.U(2,ssda_Data.Indgerm(i+1,1:imax_Indgerm)),redMarker);
         set(h(i+1),'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color',map(i+1,:),'MarkerFaceColor',map(i+1,:));
         lestr = [ lestr ' ''Step #' num2str(i) ''' ' ];
      end
      lestr = [ lestr '}' ];
%       if plot_thresholds
%          for istep = 1:length(data.y)
%             hgth(istep) = plot(u1gth,u2gth(istep,:),'-');
%             set(hgth(istep),'LineWidth',2,'Color',map(istep,:))
%          end
%       end
      if auto_gca
         set(gca,'FontSize',14,'XLim',ssda_Data.U1Lim,'YLim',ssda_Data.U2Lim);
      else
         set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
         grid on
         box on
      end
      axis square
      eval([ 'legend(' lestr ');' ]);
      
      xlabel('{\itu}_1','FontSize',14);
      ylabel('{\itu}_2','FontSize',14);

end