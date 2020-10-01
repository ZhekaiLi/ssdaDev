function varargout = demogp_class_gui(varargin)
%DEMOGP_CLASS_GUI Graphical frontend for Onlinde Gaussian Process classification.
%
%	Description
%
%	DEMOGP_CLASS_GUI - provides a graphical user interface for Gaussian
%	Process classification.
%
%	The classification is for two-dimensional inputs, these inputs are
%	provided by the user.  All OGP operations can be reached using the
%	control buttons on the left of the demo window.
%
%	The addition of data to be classified is by using the left/right
%	button of the mouse over the axis. One can then click on the LEARN
%	button to infer the OGP parameters.
%
%	One can set (or reset) the model hyperparameters (such as kernel type
%	or kernel parameter) but in this case the GP inference is restarted.
%
%	This program is designed for a first contact with GP and the idea of
%	sparsity within a full probabilistic inference. Sparsity i controlled
%	by imposing an upper limit to the number of Basis Vectors or setting
%	the threshold value.
%
%	See also
%	OGP, OGPTRAIN, DEMOGP_CLASS, DEMOGP_REG_GUI
%

%	Copyright (c) Lehel Csato (2001-2004)

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

if nargin == 0  % LAUNCH GUI

  fig = figure(1234);   % generate the figure
  [X,Y] = gui_create_data(fig);
  % and the additional structures required
  setappdata(fig,'net',net);
  setappdata(fig,'X_train',X);
  setappdata(fig,'Y_train',Y);
  setappdata(fig,'ep',ep);
  setappdata(fig,'isep',1);
  [h_plot h_hyp] = gui_draw(fig,net);

  % Use system color scheme for figure:
  set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'), ...
	  'Menubar','none');
    
  % Stores the handle for the axes!
  setappdata(fig,'h_plot', h_plot);
  setappdata(fig,'h_hyp' , h_hyp );

    
  if nargout > 0
    varargout{1} = fig;
  end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

  fig     = get(gcbo,'Parent');
  net     = getappdata(fig,'net');
  X_train = getappdata(fig,'X_train');
  Y_train = getappdata(fig,'Y_train');
  ep      = getappdata(fig,'ep');

  try
    [X_train, Y_train] = ...
	feval(varargin{:},X_train,Y_train); % FEVAL switchyard
  catch
    disp(lasterr);
  end;

  % displaying the result
  show_result(getappdata(fig,'h_plot'),X_train,Y_train);
  % storing the changes
  setappdata(fig,'net',net);
  setappdata(fig,'X_train',X_train);
  setappdata(fig,'Y_train',Y_train);
  setappdata(fig,'ep',ep);
end
%%%++++++++++++++++++ END MAIN FUNCTION ++++++++++++++++++++++

%%--!-!-!-!-!-- ADDITIONAL FUNCTIONS REQUIRED FOR GUI --!-!-!-!-!--%%

% --------------------------------------------------------------------
function [X, Y] = gui_create_data(fig);
% Creates the variable NET of class ogp and initialises the necessary 
% variables for the graphical interface 

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

% building a GP
ogp(2,1, ...	      	      % input/output  dimensions
    'SQEXP', ...	      % kernel type
    log([1/10 1/10 5 3]));    % kernel parameters
ogpinit(@c_class_bin,...      % likelihood function
	1);	      	      % likelihood function parameter
% SETTING the PRIOR ::
ogphypcovpar(2e-1);
net.hypmean(end) = net.kpar(end);
net.thresh    = 1e-3;	      % the admission threshold for new BVs
net.maxBV     = 100;	      % maximal number of BV-s
net.bias      = -500;	      % the log! of bias

X  = ones(0,net.nin);
Y  = ones(0,net.nout);
ep = [];
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X,Y]=gui_learn_ogp(h_curr,X,Y);
% function that iterates the training procedure for the online GP. The
% function performs a single iteration.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

% building the GPOPT vector
gpopt              = defoptions;
gpopt.postopt.isep = getappdata(gcf,'isep');
gpopt.postopt.itn   = 2;
gpopt.postopt.fixitn= 2;

gpopt.covopt.opt   = 10;
gpopt.covopt.fnopt = 'conjgrad';

% performing learning
if size(X,1);
  ogpreset;
  ogptrain(X,Y);
  ogpreset;
  ogppost(X,Y);
end;
% update the fields displaying the MODEL parameters
gui_par_update(ogppak);
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_reset_ogp(h_curr, X, Y);
% function for resetting the onlie Gaussian Process NET.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

ogpreset;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function gui_par_update(pars);
% function that updates the graphical TAGs corresponding to the
% hyperparameters.
  
h_hyp = getappdata(gcf,'h_hyp');
pars  = exp(pars);

% in pars: (Bias ... kpar(1), kpar(2), kpar(3), kpar(4))
% in h_hyp:  (kpar(3), kpar(4), kpar(1), kpar(2));
set(h_hyp(1),'String',num2str(pars(4)));
set(h_hyp(2),'String',num2str(pars(5)));
set(h_hyp(3),'String',num2str(pars(2)));
set(h_hyp(4),'String',num2str(pars(3)));
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_add_point(h_curr,X,Y);
% function that adds a point to the training data. If left button is
% pressed, then the new data belongs to CLASS1, otherwise to CLASS2

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

N_train = length(Y);

fprintf(['\nAdding training DATA:\n\nLeft button over the axis adds ' ...
	 'the mouse position to CLASS_1,\n right button to CLASS_2.' ...
	 ' \n\nInput addition is ended if ANY key is pressed\n\n'])

hold on;
h_title = get(gca,'title');
old_str = get(h_title,'string');
set(h_title,...
    'string','Buttons: Left- cl.1; Right- cl.2; any key - FINISH');
set(h_title,'FontSize',14);

X_new = []; Y_new = [];
lt = 0;
while ~waitforbuttonpress;
  lt = lt + 1;
  tp = get(gca,'CurrentPoint');
  if max(abs(tp))<5;
    X_new(lt,:) = tp(1,1:2);  
  end;
  Button = upper(get(gcf,'SelectionType'));
  if Button(1) == 'N';
    Y_new(lt,1) = -1;
    ck        = 'r+';
  else
    Y_new(lt,1) = 1;
    ck        = 'bd';
  end;
  plot(X_new(lt,1),X_new(lt,2),ck,'Linewidth',1);
end;

% updating EP structure -- if there is one
if ~isempty(ep) & isfield(ep,'lamP');
  lOld = length(ep.aP);
  lNew = lt + lOld;
  ep.aP((lOld+1):lNew) = 0;
  ep.lamP((lOld+1):lNew,(lOld+1):lNew) = 0;
  ep.projP((lOld+1):lNew,:) = 0;
  ep.X((lOld+1):lNew,:)     = X_new;
end;
% concatenation
X = [X; X_new];
Y = [Y; Y_new];
% putting the title back.
set(h_title,'string',old_str);

%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X,Y] = gui_clear_points(h_curr,X,Y);
% function that removes all data points from the application.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

fprintf(['\nClearing ! ' num2str(length(Y)) ' ! training data\n\n']);

X  = zeros([0 2]);
Y  = zeros([0 1]);
ep = [];
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X,Y] = gui_newfig(h_curr,X,Y);
% function that removes all data points from the application.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

figure; h_plot = axes; box on;
set(h_plot,'Xlim',[-5 5]);
set(h_plot,'Ylim',[-5 5]);
show_result(h_plot,X,Y);
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X,Y] = gui_surfprob(h_curr,X,Y);
% function that removes all data points from the application.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

figure(234);
% constants to be used
xG = 151; yG = 151;
xLim = [-5 5]; yLim = [-5 5];
xP      = linspace(xLim(1), xLim(2), xG);
yP      = linspace(yLim(1), yLim(2), yG);
[xT,yT] = meshgrid(xP,yP);
[pP,pV] = ogpfwd([xT(:) yT(:)]);
pP      = reshape(pP./pV,[xG yG]);
pP      = (1 + erf(pP./sqrt(2)))./2;
h       = surf(xP,yP,pP);
set(h,'Linewidth',0.1,'EdgeAlpha',0.2);
axis tight;
xlabel('X_1'); ylabel('X_2');
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function show_result(h_plot,X,Y);
% show_res(handle,net,X,Y) - plots the training inputs (two classes)
% together with a decision boundary (drawn from the posterior) into the
% axes object with handle h_plot

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

% constants to be used
xG = 51; yG = 51;
xLim = get(h_plot,'Xlim');
yLim = get(h_plot,'Ylim');

% rendering the new plot
axes(h_plot); cla; hold on;
pos = find(Y>0); neg = find(Y<0);
h_pos = scatter(X(pos,1),X(pos,2),'bd');
h_neg = scatter(X(neg,1),X(neg,2),'r+');
set(h_pos,'Linewidth',3);
set(h_neg,'Linewidth',3);
if size(net.BV,1);	      % contour plot only if needed
  % generating surface for post. class-prob.
  xP      = linspace(xLim(1), xLim(2), xG);
  yP      = linspace(yLim(1), yLim(2), yG);
  [xT,yT] = meshgrid(xP,yP);
  [pP,pV] = ogpfwd([xT(:) yT(:)]);
  pP      = reshape(pP./pV,[xG yG]);
  pP      = (1 + erf(pP./sqrt(2)))./2;
  [cs,ch] = contour(xP,yP,pP,[0.5 0.5]);
  if length(ch);
    clabel(cs,ch,'Fontsize',13);
    set(ch,'Linewidth',4,'EdgeColor',[0.5 0.1 0.1]);
  end;
  [cs,ch] = contour(xP,yP,pP,[0.3 0.7]);
  if length(ch);
    clabel(cs,ch,'Fontsize',10);
    set(ch,'Linewidth',2,'EdgeColor',[0.3 0.8 0.3]);
  end;
  h_BV    = scatter(net.BV(:,1),net.BV(:,2),'ko');
  set(h_BV,'Linewidth',2,'Markersize',12);
end;
title(['Data#: ' num2str(length(Y)) ', #BV: ' ... 
       num2str(size(net.BV,1)) ', Neg.Log.Evid: ' ...
       num2str(ogpevid([]))]);
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [h_plot, h_allH] = gui_draw(fig,net)
% LOOOONG function for defining the outlook of the graphical user
% interface. The function returns the handle for the axis where to plot
% the points and the result of the online GP classification

bg_col = [0.65 0.6 0.5];

set(fig,'Units','normalized','NumberTitle','off', ...
	'Name','Online Gaussian Process Classification');
gui_frame2= uicontrol(fig,'Style','frame',...
		     'Units','normalized', ... 
		     'Position',[0.005 0.005 0.24 0.99], ...
		     'BackgroundColor',bg_col);
gui_frame = uicontrol(fig,'Style','frame',...
		     'Units','normalized', ...
		     'Position',[0.01 0.34 0.23 0.59], ...
		     'BackgroundColor',bg_col);
gui_text  = uicontrol(fig,'Style','text','String','Parameters', ...
		     'BackgroundColor',bg_col);
set(gui_text,'Units','normalized','Position',[0.04 0.93 0.17 0.06], ...
	    'Fontsize',14,'HorizontalAlignment','center');

% TYPE of kernel
gui_ker_text= uicontrol(fig,'Style','text','String','Type:', ...
		       'Units','normalized', ...
		       'Position',[0.015 0.86 0.07 0.04], ...
		       'Fontsize',12,'HorizontalAlignment','left', ...
		       'BackgroundColor',bg_col);
allKers    = strvcat('SQEXP','POLY','RATQUAD','MATERN');
ker_val    = strmatch(upper(net.covarfn),allKers);
gui_kertype = ...
    uicontrol(fig,'Style','popup','String',allKers, ... 
	      'Value',ker_val,'Units','normalized', ... 
	      'Position',[0.08 0.86 0.15 0.05],...
	      'Callback','demogp_class_gui(''gui_set_ktype'',gcbo)', ...
	      'BackgroundColor',bg_col);
gui_par_text = ...
    uicontrol(fig,'Style','text','String','K. parameters', ...
	      'Units','normalized', 'Backgroundcolor','red', ...
	      'Position',[0.02 0.80 0.21 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','center',...
	      'BackgroundColor',bg_col);

% AMPLITUDE of kernel
gui_amp_text = ...
    uicontrol(fig,'Style','text','String','Amp:', ...
	      'Units','normalized', ...
	      'Position',[0.02 0.75 0.10 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left', ...
	      'BackgroundColor',bg_col);
h_allH(1)  = ...
    uicontrol(fig,'Style','edit', ...
	      'String',num2str(exp(net.kpar(3))), ...
	      'Units','normalized', ...
	      'Position',[0.02 0.70 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_class_gui(''gui_set_amp'',gcbo)');

% ORDER of polynomial
gui_ord_text = ...
    uicontrol(fig,'Style','text','String','Ord:', ...
	      'Units','normalized', ...
	      'Position',[0.12 0.75 0.10 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'BackgroundColor',bg_col);
h_allH(2)  = ...
    uicontrol(fig,'Style','edit','String',num2str(exp(net.kpar(4))), ...
	      'Units','normalized', ...
	      'Position',[0.13 0.70 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_class_gui(''gui_set_ord'',gcbo)');
% SCALING of the inputs
gui_sc_text = ...
    uicontrol(fig,'Style','text','String','Input weight:', ...
	      'Units','normalized', ...
	      'Position',[0.02 0.65 0.20 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'BackgroundColor',bg_col);
h_allH(3) = ...
    uicontrol(fig,'Style','edit', ... 
	      'String',num2str(exp(net.kpar(1))), ...
	      'Units','normalized', ...
	      'Position',[0.02 0.60 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_class_gui(''gui_set_sc1'',gcbo)');
h_allH(4) = ...
    uicontrol(fig,'Style','edit', ... 
	      'String',num2str(exp(net.kpar(2))), ...
	      'Units','normalized', ...
	      'Position',[0.13 0.60 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_class_gui(''gui_set_sc2'',gcbo)');

% MAXIMAL # of BVs
gui_mbv_text = ...
    uicontrol(fig,'Style','text','String','BV#:', ...
	      'Units','normalized', ...
	      'Position',[0.02 0.55 0.10 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'BackgroundColor',bg_col);
gui_mbv_val  = ...
    uicontrol(fig,'Style','edit','String',num2str(net.maxBV), ...
	      'Units','normalized', ...
	      'Position',[0.02 0.50 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_class_gui(''gui_set_mbv'',gcbo)');

% THRESHOLD for adding new BVs
gui_thr_text = ...
     uicontrol(fig,'Style','text','String','Thrs.:', ...
	       'Units','normalized', ...
	       'Position',[0.13 0.55 0.10 0.04], ...
	       'Fontsize',12,'HorizontalAlignment','left',...
	       'BackgroundColor',bg_col);
gui_thr_val  = ...
    uicontrol(fig,'Style','edit','String',num2str(net.thresh), ...
	      'Units','normalized', ...
	      'Position',[0.13 0.50 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left', ...
	      'Callback','demogp_class_gui(''gui_set_thr'',gcbo)');
% OUTPUT NOISE 
gui_out_text = ...
    uicontrol(fig,'Style','text','String','Out. noise:', ...
	      'Units','normalized', ...
	      'Position',[0.02 0.45 0.13 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'BackgroundColor',bg_col);
gui_out_val  = ...
    uicontrol(fig,'Style','edit','String',num2str(net.likpar(1)), ...
	      'Units','normalized', ...
	      'Position',[0.145 0.45 0.085 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_class_gui(''gui_set_out'',gcbo)');

% EPCONTROL button
isep       = getappdata(gcf,'isep');
gui_ep_val  = ...
    uicontrol(fig,'Style','checkbox','String','TAP/EP', ...
	      'Units','normalized','Value',isep, ...
	      'Position',[0.02 0.40 0.21 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_class_gui(''gui_set_ep'',gcbo)');

% FIXEDBVCONTROL button
gui_ep_val  = ...
    uicontrol(fig,'Style','checkbox','String','Fixed BV Set', ...
	      'Units','normalized','Value',net.isBVfixed, ...
	      'Position',[0.02 0.35 0.21 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_class_gui(''gui_set_fixed'',gcbo)');

%%%%%%%%%%%% EXECUTION CONTROL %%%%%%%%%%%%

gui_add_point = ...
    uicontrol(fig,'Style','pushbutton','String','Add pts.','Fontsize',12, ...
	      'Units','normalized','Position',[0.01 0.245 0.115 0.06], ...
	      'HorizontalAlignment','center', ...
	      'Callback','demogp_class_gui(''gui_add_point'',gcbo)', ...
	      'TooltipString',...
	      'New training points: left (red) class 1; right (blue) class 2');

gui_del_point = ...
    uicontrol(fig,'Style','pushbutton','String','Clear pts.',...
	      'Units','normalized', ...
	      'Position',[0.125 0.245 0.115 0.06], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','demogp_class_gui(''gui_clear_points'',gcbo)',...
	      'ToolTipString',['Removes ALL TRAINING data' ...
		    ' without resetting the Gaussian Process']);

gui_learn = ...
    uicontrol(fig,'Style','pushbutton','String','Learn!!',...
	      'Units','normalized', ...
	      'Position',[0.01 0.125 0.115 0.12], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','demogp_class_gui(''gui_learn_ogp'',gcbo)', ...
	      'BackgroundColor','cyan', ...
	      'ToolTipString',...
	      sprintf(['Training GP: ' ...
		    'obtaining a posterior\n and' ...
		    ' adjusting hyperparameters']));

gui_reset = ...
    uicontrol(fig,'Style','pushbutton','String','Reset',...
	      'Units','normalized', ...
	      'Position',[0.125 0.125 0.115 0.12], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','demogp_class_gui(''gui_reset_ogp'',gcbo)', ...
	      'ToolTipString','Resets the state of the OGP');

gui_newfig = ...
    uicontrol(fig,'Style','pushbutton','String','Sep. Fig',...
	      'Units','normalized', ...
	      'Position',[0.01 0.065 0.115 0.06], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','demogp_class_gui(''gui_newfig'',gcbo)', ...
	      'ToolTipString','Plots the figure in a separate window');

gui_surfprob = ...
    uicontrol(fig,'Style','pushbutton','String','Surf. Prob',...
	      'Units','normalized','Fontsize',12, ...
	      'Position',[0.125 0.065 0.115 0.06], ...
	      'HorizontalAlignment','center', ...
	      'Callback','demogp_class_gui(''gui_surfprob'',gcbo)', ...
	      'ToolTipString',...
	      'Plots the class-cond. probability surface');

gui_close = ...
    uicontrol(fig,'Style','pushbutton','String','Close',...
	      'Units','normalized', ...
	      'Position',[0.01 0.01 0.23 0.055], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','close',...
	      'BackgroundColor', [0.8 0.4 0.4], ...
	      'ToolTipString','Closes the DEMO');

% DEFINING the axes ...
h_plot = axes('Position',[0.29 0.05 0.67 0.89], ...
	      'Box','on','Units','normalized', ...
	      'Fontsize',12, ...
	      'Xlim',[-5 5],'Ylim',[-5 5]);
title('BV#: 0');
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

%%
%% Section for setting the various OGP parameters.
%%

% --------------------------------------------------------------------
function [X, Y] = gui_set_ktype(h_curr,X,Y);
% setting kernel type

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

oldcovfn = net.covarfn;
% retrieving value of POPUP
switch get(h_curr,'Value')
 case 1;
  net.covarfn = 'SQEXP';
 case 2,
  net.covarfn = 'POLY';
 case 3;
  net.covarfn = 'RATQUAD';
 case 4;
  net.covarfn = 'USER';
  net.kfnaddr   = @cov_matern;
  net.gradkaddr = @covgrad_matern;
end;

%!!! if the new and old covarfns are different, then the OGP
%!!! structure needs resetting
if ~strcmp(upper(net.covarfn),upper(oldcovfn));
  ogpreset;
  if strcmp(upper(oldcovfn),'USER');
    net.kfnaddr   = [];
    net.gradkaddr = [];
  end;
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_amp(h_curr,X,Y);
% setting kernel amplitude

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

valStr   = get(h_curr,'String');
amp      = str2num(valStr);
if amp>0;
  net.kpar(3) = log(amp);
  ogpreset;
else
  tStr = num2str(net.kpar(3));
  set(h_curr,'String',tStr);
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_ord(h_curr,X,Y);
% Setting the order parameter of the kernel (order of the polynomial
% and the rational quadratic kernel

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

valStr   = get(h_curr,'String');
ord      = str2num(valStr);
if ord>=0;
  if strcmp(upper(net.covarfn),'POLY');
    ord = ceil(ord);
  end;
  net.kpar(4) = log(ord);
  ogpreset;
else
  tStr = num2str(net.kpar(2));
  set(h_curr,'String',tStr);
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_sc1(h_curr,X,Y);
% setting horizontal scaling in the input space

valStr   = get(h_curr,'String');
sc      = str2num(valStr);
if sc>0;
  net.kpar(1) = log(sc);
  ogpreset;
else
  set(h_curr,'String',num2str(exp(net.kpar(1))));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_sc2(h_curr,X,Y);
% setting vertical scaling of the input space

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

valStr   = get(h_curr,'String');
sc      = str2num(valStr);
if sc>0;
  net.kpar(2) = log(sc);
  ogpreset;
else
  set(h_curr,'String',num2str(exp(net.kpar(2))));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_mbv(h_curr,X,Y);
% setting the maximal number of BVs for the OGP

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

valStr   = get(h_curr,'String');
sc      = str2num(valStr);
if sc>=1;
  net.maxBV = ceil(sc);
else
  set(h_curr,'String',num2str(net.maxBV));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_thr(h_curr,X,Y);
% setting the threshold for new input addition in the OGP

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

valStr   = get(h_curr,'String');
sc      = str2num(valStr);
if sc>0;
  net.thresh = sc;
else
  set(h_curr,'String',num2str(net.thresh));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_out(h_curr,X,Y);
% Setting output noise level

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

valStr   = get(h_curr,'String');
sc      = str2num(valStr);
if sc>=0;
  net.likpar(1) = sc;
else
  set(h_curr,'String',num2str(net.likpar(1)));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_ep(h_curr,X,Y);
% Setting the EP-FLAG: whether to perform or not the EP iterations

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

val      = get(h_curr,'Value');
isep = val;
% saving the indicator
setappdata(gcf,'isep',isep);
if ~isep;		      % removing the EP structure
  ep = [];
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_fixed(h_curr,X,Y);
% Setting the EP-FLAG: whether too perform or not the EP iterations

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

val      = get(h_curr,'Value');
net.isBVfixed = val;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++
