function varargout = demogp_reg_gui(varargin)
%DEMOGP_REG_GUI Graphical frontend to Online Gaussian Process regression.
%
%	Description
%
%	DEMOGP_REG_GUI - provides a graphical user interface for Online
%	Gaussian Process regression.
%
%	This demo illustrates the flexibility of the Online Gaussian
%	Processes for regression using additive noise with Gaussian, Laplace,
%	or positive-exponential distributions. All OGP operations can be
%	reached using the control buttons on the left of the window.
%
%	One can set (or reset) the model hyperparameters (such as kernel type
%	or kernel parameter) but in this case the GP inference is restarted.
%
%	This program is designed for a first contact with GP and the idea of
%	sparsity within a full probabilistic inference.
%
%	See also
%	OGP, OGPTRAIN, DEMOGP_REG, DEMOGP_CLASS_GUI
%

%	Copyright (c) Lehel Csato (2001-2004)

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

if nargin == 0  % LAUNCH GUI

  fig = figure(4321);   % generate the figure

  % initialise global variables 
  [X,Y] = gui_create_data(fig);
  t_noise      = 'GAUSS';
  t_var        = 1;
  t_num        = 10;
  x_range      = 10;
  y_range      = 10;

  % and the additional structures required
  setappdata(fig,'net',net);
  setappdata(fig,'X_train',X);
  setappdata(fig,'Y_train',Y);
  setappdata(fig,'ep',ep);
  setappdata(fig,'isep',1);
  setappdata(fig,'isem',1);
  setappdata(fig,'t_noise',t_noise);
  setappdata(fig,'t_var',t_var);
  setappdata(fig,'t_num',t_num);
  setappdata(fig,'x_range',x_range);
  setappdata(fig,'y_range',y_range);

  [h_plot, h_hyp] = gui_draw(fig);

  % Use system color scheme for figure:
  set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'), ...
	  'Menubar','none');

  % Stores the handle for the axes!
  setappdata(fig,'h_plot', h_plot);
  setappdata(fig,'h_hyp' , h_hyp );

  show_result(h_plot,X,Y);

  if nargout > 0
    varargout{1} = fig;
  end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

  fig     = get(gcbo,'Parent');
  net     = getappdata(fig,'net');
  X_train = getappdata(fig,'X_train');
  Y_train = getappdata(fig,'Y_train');
  ep      = getappdata(fig,'ep');
    
%  try
    [X_train, Y_train] = ...
	feval(varargin{:},X_train,Y_train); % FEVAL switchyard
%  catch
%    disp(lasterr);
%  end;
    
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
ogp(1,1, ...	    	      % output  dimension
    'SQEXP', ...	      % kernel type
    log([1 10 3]));  	      % kernel parameter
ogpinit(@c_reg_gauss,...      % likelihood function
	1,...	      	      % likelihood function parameter
	@em_gauss);           % lik.opt. parameters
net.thresh    = 1e-3;	      % the admission threshold for new BVs
net.maxBV     = 40;	      % maximal number of BV-s
% SETTING the PRIOR ::
ogphypcovpar(3e-1);

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
gpopt               = defoptions;
gpopt.postopt.isep  = getappdata(gcf,'isep');
gpopt.postopt.itn   = 2;
gpopt.postopt.fixitn= 2;

isem = getappdata(gcf,'isem');
if ~isem;
  saveAddr     = net.likoptfn;
  net.likoptfn = [];
end;

% performing learning
if size(X,1);
  ogpreset;
  ogptrain(X,Y);
  ogpreset;
  ogppost(X,Y);
end;

if ~isem;
  net.likoptfn = saveAddr;
end;
% update the fields displaying the MODEL parameters
gui_par_update([exp(ogppak), net.likpar(1)]);
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [net, X, Y, ep] = gui_reset_ogp(h_curr, net, X, Y, ep);
% function for resetting the onlie Gaussian Process NET.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

ogpreset;

%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function gui_par_update(pars);
% function that updates the graphical TAGs corresponding to the
% hyperparameters.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

h_hyp = getappdata(gcf,'h_hyp');

% in pars: (Bias ... kpar(1) kpar(2), kpar(3))
% in h_hyp:  (kpar(2), kpar(3), kpar(1));

set(h_hyp(1),'String',num2str(pars(3)));
set(h_hyp(2),'String',num2str(pars(4)));
set(h_hyp(3),'String',num2str(pars(2)));
set(h_hyp(4),'String',num2str(pars(5)));
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_add_point(h_curr,X,Y);
% function that adds a point to the training data. If left button is
% pressed, then the new data belongs to CLASS1, otherwise to CLASS2

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

t_num   = getappdata(gcf,'t_num');
x_range = getappdata(gcf,'x_range');
y_range = getappdata(gcf,'y_range');
t_var   = getappdata(gcf,'t_var');
t_noise = getappdata(gcf,'t_noise');

% generating the data
[X_new,Y_new] = sincdata(t_noise,t_num,t_var,[],[],x_range,y_range);

% updating EP structure -- if there is one
if ~isempty(ep) & isfield(ep,'lamP');
  lOld = length(ep.aP);
  lNew = t_num + lOld;
  ep.aP((lOld+1):lNew) = 0;
  ep.lamP((lOld+1):lNew,(lOld+1):lNew) = 1e-20*speye(t_num);
  ep.projP((lOld+1):lNew,:) = 0;
  ep.X((lOld+1):lNew,:)     = X_new;
  ep.logZ((lOld+1):lNew,:)  = 1;
end;
% concatenation
X = [X; X_new];
Y = [Y; Y_new];
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X,Y] = gui_clear_points(h_curr,X,Y);
% function that removes all data points from the application.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

fprintf(['\nClearing ! ' num2str(length(Y)) ' ! training data\n\n']);

X  = zeros([0 1]);
Y  = zeros([0 1]);
ep = [];
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X,Y] = gui_newfig(h_curr,X,Y);
% function that removes all data points from the application.

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

x_range = getappdata(gcf,'x_range');
y_range = getappdata(gcf,'y_range');
figure; h_plot = axes; box on;
set(h_plot,'Xlim',[-x_range x_range]);
set(h_plot,'Ylim',[-y_range*.5 y_range*1.4]);
setappdata(gcf,'x_range',x_range);
setappdata(gcf,'y_range',y_range);
show_result(h_plot,X,Y);
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function show_result(h_plot,X,Y);
% show_res(handle,net,X,Y) - plots the training inputs (two classes)
% together with a decision boundary (drawn from the posterior) into the
% axes object with handle h_plot

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

% resolution of the plot
xG = 201;
x_range     = getappdata(gcf,'x_range');
y_range     = getappdata(gcf,'y_range');
[xP,yTest]  = sincdata('GAUSS',xG,0,[],1,x_range,y_range);

% rendering the new plot
axes(h_plot); cla; hold on;
h_true  = plot(xP,yTest,'k-.','Linewidth',2);
h_train = plot(X,Y,'b+');
set(h_train,'Linewidth',2);

[yP,yV] = ogpfwd(xP);     % computing the predictions
yV      = sqrt(yV);	      % adjusting to get the deviation
h_meanp = plot(xP,yP,'g--','Linewidth',3);
h_stdP  = plot(xP,yP+yV,'r-','Linewidth',2);
          plot(xP,yP-yV,'r-','Linewidth',2);

[bvP, bvV] = ogpfwd(net.BV);
h_BV    = plot(net.BV,bvP,'m^');
set(h_BV,'Linewidth',4,'Markersize',12);
if length(net.BV);
  title(['Data#: ' num2str(size(X,1)) ', BV#: ' ...
	 num2str(size(net.BV,1)) ', Neg.Log.Evid: ' ...
	 num2str(ogpevid([]))]);
else
  title(['Data#: ' num2str(size(X,1)) ', BV#: 0']);
end;

legend([h_true, h_meanp,h_stdP, h_train, h_BV], ...
       'True function','Mean function', ...
       'Bay. Error bars','Noisy points','BV positions',1);
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [h_plot, h_allH] = gui_draw(fig)
% LOOOONG function for defining the outlook of the graphical user
% interface. The function returns the handle for the axis where to plot
% the points and the result of the online GP classification

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

bg_col = [0.65 0.6 0.5];

set(fig,'Units','normalized','NumberTitle','off', ...
	'Name','Online Gaussian Process Regression');
gui_frame2= uicontrol(fig,'Style','frame',...
		     'Units','normalized', ... 
		     'Position',[0.005 0.005 0.24 0.99], ...
		     'BackgroundColor',bg_col);
gui_frame = uicontrol(fig,'Style','frame',...
		     'Units','normalized', ...
		     'Position',[0.01 0.49 0.23 0.44], ...
		     'BackgroundColor',bg_col);
gui_text  = uicontrol(fig,'Style','text','String','Parameters', ...
		     'BackgroundColor',bg_col);
set(gui_text,'Units','normalized','Position',[0.04 0.93 0.17 0.055], ...
	    'Fontsize',14,'HorizontalAlignment','center');

% TYPE of kernel
gui_ker_text= uicontrol(fig,'Style','text','String','Type:', ...
		       'Units','normalized', ...
		       'Position',[0.015 0.86 0.07 0.04], ...
		       'Fontsize',12,'HorizontalAlignment','left', ...
		       'BackgroundColor',bg_col);
allKers    = strvcat('SQEXP','POLY','RATQUAD');
ker_val    = strmatch(net.covarfn,allKers);
gui_kertype = ...
    uicontrol(fig,'Style','popup','String',allKers, ... 
	      'Value',ker_val,'Units','normalized', ... 
	      'Position',[0.08 0.86 0.15 0.05],...
	      'Callback','demogp_reg_gui(''gui_set_ktype'',gcbo)', ...
	      'BackgroundColor',bg_col);

% AMPLITUDE of the kernel
gui_amp_text = ...
    uicontrol(fig,'Style','text','String','Amp:', ...
	      'Units','normalized', ...
	      'Position',[0.02 0.81 0.10 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left', ...
	      'BackgroundColor',bg_col);
h_allH(1)  = ...
    uicontrol(fig,'Style','edit', ...
	      'String',num2str(exp(net.kpar(2))), ...
	      'Units','normalized', ...
	      'Position',[0.02 0.77 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_amp'',gcbo)');

% ORDER of the polynomial
gui_ord_text = ...
    uicontrol(fig,'Style','text','String','Ord:', ...
	      'Units','normalized', ...
	      'Position',[0.12 0.81 0.10 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'BackgroundColor',bg_col);
h_allH(2)  = ...
    uicontrol(fig,'Style','edit','String',num2str(exp(net.kpar(3))), ...
	      'Units','normalized', ...
	      'Position',[0.13 0.77 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_ord'',gcbo)');

% input SCALING 
gui_sc_text = ...
    uicontrol(fig,'Style','text','String','Input sc.:', ...
	      'Units','normalized', ...
	      'Position',[0.02 0.71 0.12 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'BackgroundColor',bg_col);
h_allH(3) = ...
    uicontrol(fig,'Style','edit', ... 
	      'String',num2str(exp(-net.kpar(1))), ...
	      'Units','normalized', ...
	      'Position',[0.13 0.71 0.10 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_sc1'',gcbo)');

% LIKELIHOOD NOISE 

% DEFINITION OF THE NOISE POPUP
all_types = strvcat('GAUSS','LAPLACE','POSEXP');

gui_nout_text = ...
    uicontrol(fig,'Style','text','String','Noise:', ...
	      'Units','normalized', ...
	      'Position',[0.02 0.645 0.21 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'BackgroundColor',bg_col);

% NOISECONTROL button
isem       = getappdata(gcf,'isem');
gui_ep_val  = ...
    uicontrol(fig,'Style','checkbox','String','Est.', ...
	      'Units','normalized','Value',isem, ...
	      'Position',[0.13 0.645 0.1 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_em'',gcbo)');

n_type = func2str(net.likaddr);
l_type = length(n_type); n_type = upper(n_type(7:l_type));
type_val  = strmatch(n_type,all_types);
gui_nout_type = ...
    uicontrol(fig,'Style','popup','String',all_types, ... 
	      'Value',type_val,'Units','normalized', ... 
	      'Position',[0.02 0.595 0.13 0.05],...
	      'Callback','demogp_reg_gui(''gui_set_likout'',gcbo)', ...
	      'BackgroundColor',bg_col);

h_allH(4)  = ...
    uicontrol(fig,'Style','edit','String',num2str(net.likpar(1)), ...
	      'Units','normalized', ...
	      'Position',[0.15 0.595 0.08 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_out'',gcbo)');

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
	      'Callback','demogp_reg_gui(''gui_set_mbv'',gcbo)');

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
	      'Callback','demogp_reg_gui(''gui_set_thr'',gcbo)');

% EPCONTROL button
isep       = getappdata(gcf,'isep');
gui_ep_val  = ...
    uicontrol(fig,'Style','checkbox','String','TAP/EP', ...
	      'Units','normalized','Value',isep, ...
	      'Position',[0.02 0.42 0.21 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_ep'',gcbo)');

% FIXEDBVCONTROL button
gui_ep_val  = ...
    uicontrol(fig,'Style','checkbox','String','Fixed BV Set', ...
	      'Units','normalized','Value',net.isBVfixed, ...
	      'Position',[0.02 0.36 0.21 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_fixed'',gcbo)');

%%%%%%%%%%%% EXECUTION CONTROL %%%%%%%%%%%%

gui_learn = ...
    uicontrol(fig,'Style','pushbutton','String','Learn!!',...
	      'Units','normalized', ...
	      'Position',[0.01 0.250 0.115 0.10], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','demogp_reg_gui(''gui_learn_ogp'',gcbo)', ...
	      'BackgroundColor', [0.4 0.9 0.8], ...
	      'ToolTipString',...
	      sprintf(['Training of the OGP structure:\n' ...
	       'obtaining a posterior and\n' ...
	       ' adjusting the hyperparameters']));

gui_reset = ...
    uicontrol(fig,'Style','pushbutton','String','Reset',...
	      'Units','normalized', ...
	      'Position',[0.125 0.250 0.115 0.10], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','demogp_reg_gui(''gui_reset_ogp'',gcbo)', ...
	      'ToolTipString','Resets the state of the OGP');

%%%%%%%%%%%% DATA GENERATION %%%%%%%%%%%%
bg_data  = [0.7 0.6 0.7];

gui_frame = uicontrol(fig,'Style','frame',...
		     'Units','normalized', ...
		     'Position',[0.01 0.07 0.23 0.17], ...
		     'BackgroundColor',bg_data);
gui_ntrue_text = ...
    uicontrol(fig,'Style','text','String','Data (true) noise:', ...
	      'Units','normalized', ...
	      'Position',[0.02 0.12 0.21 0.04], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'BackgroundColor',bg_data);

n_type    = upper(getappdata(gcf,'t_noise'));
type_val  = strmatch(n_type,all_types);
gui_ntrue_type = ...
    uicontrol(fig,'Style','popup','String',all_types, ... 
	      'Value',type_val,'Units','normalized', ... 
	      'Position',[0.02 0.075 0.13 0.05],...
	      'Callback','demogp_reg_gui(''gui_set_t_noise'',gcbo)', ...
	      'BackgroundColor',bg_data,...
	      'ToolTipString','Noise type');

ttt = getappdata(gcf,'t_var');
gui_ntrue_val  = ...
    uicontrol(fig,'Style','edit','String',num2str(ttt), ...
	      'Units','normalized', ...
	      'Position',[0.15 0.075 0.08 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_t_var'',gcbo)');

gui_add_point = ...
    uicontrol(fig,'Style','pushbutton','String','Add', ...
	      'Fontsize',12, 'Units','normalized',...
	      'Position',[0.02 0.18 0.08 0.05], ...
	      'HorizontalAlignment','center', ...
	      'Callback','demogp_reg_gui(''gui_add_point'',gcbo)', ...
	      'TooltipString',...
	      sprintf('Adds noisy training inputs\nto the training set'), ...
	      'BackgroundColor',bg_data);

gui_del_point = ...
    uicontrol(fig,'Style','pushbutton','String','Clear',...
	      'Units','normalized', ...
	      'Position',[0.16 0.18 0.07 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','demogp_reg_gui(''gui_clear_points'',gcbo)',...
	      'ToolTipString',['Removes ALL TRAINING data' ...
		    ' without resetting the Gaussian Process'], ...
	      'BackgroundColor',bg_data);

ttt = getappdata(gcf,'t_num');
gui_ntrain_val  = ...
    uicontrol(fig,'Style','edit','String',num2str(ttt), ...
	      'Units','normalized', ...
	      'Position',[0.10 0.18 0.06 0.05], ...
	      'Fontsize',12,'HorizontalAlignment','left',...
	      'Callback','demogp_reg_gui(''gui_set_t_num'',gcbo)');

% Button to display the result in a different window
gui_newfig = ...
    uicontrol(fig,'Style','pushbutton','String','Sep. Fig',...
	      'Units','normalized', ...
	      'Position',[0.01 0.01 0.115 0.055], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','demogp_reg_gui(''gui_newfig'',gcbo)', ...
	      'ToolTipString','Plots the figure in a separate window');

% DUMMY BUTTON TO CLOSE THE WINDOW
gui_close = ...
    uicontrol(fig,'Style','pushbutton','String','Close',...
	      'Units','normalized', ...
	      'Position',[0.125 0.01 0.115 0.055], ...
	      'Fontsize',12,'HorizontalAlignment','center', ...
	      'Callback','close',...
	      'BackgroundColor', [0.8 0.4 0.4], ...
	      'ToolTipString','Closes the DEMO');

% DEFINING the axes ...
x_range = getappdata(gcf,'x_range');
y_range = getappdata(gcf,'y_range');
h_plot  = axes('Position',[0.29 0.05 0.67 0.89], ...
	       'Box','on','Units','normalized', ...
	       'Fontsize',12, ...
	       'Xlim',[-x_range x_range],...
	       'Ylim',[-y_range*.5 y_range*1.4]);

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
c_val       = get(h_curr,'Value');
net.covarfn = get(h_curr,'String');
net.covarfn = deblank(net.covarfn(c_val,:));
%!!! if new and old covarfns differ, the OGP needs resetting
if ~strcmp(upper(net.covarfn),upper(oldcovfn));
  ogpreset;
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_amp(h_curr,X,Y);
% setting kernel amplitude

valStr   = get(h_curr,'String');
amp      = str2num(valStr);
if amp>0;
  net.kpar(2) = log(amp);
  ogpreset;
else
  set(h_curr,'String',num2str(exp(net.kpar(2))));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_ord(h_curr,X,Y);
% Setting the order parameter of the kernel (for the polynomial
% and rational quadratic kernels)

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

valStr   = get(h_curr,'String');
ord      = str2num(valStr);
if ord>=0;
  if strcmp(upper(net.covarfn),'POLY');
    ord = ceil(ord);
  end;
  net.kpar(3) = log(ord);
  ogpreset;
else
  set(h_curr,'String',num2str(exp(net.kpar(3))));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_sc1(h_curr,X,Y);
% setting horizontal scaling in the input space

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

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
% Setting the EP-FLAG: whether too perform or not the EP iterations

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
function [X, Y] = gui_set_em(h_curr,X,Y);
% Setting the EP-FLAG: whether too perform or not the EP iterations

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

val      = get(h_curr,'Value');
isem = val;
% saving the indicator
setappdata(gcf,'isem',isem);
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_fixed(h_curr,X,Y);
% Setting the indicator for FIXED BV set

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

val      = get(h_curr,'Value');
net.isBVfixed = val;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_t_noise(h_curr,X,Y);
% Setting the true noise TYPE (Gauss/Lapl/Posexp)

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

val      = get(h_curr,'Value');
all_T    = get(h_curr,'String');
setappdata(gcf,'t_noise',deblank(all_T(val,:)));
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_t_var(h_curr,X,Y);
% Setting the "variance" of the noise 

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

val      = str2num(get(h_curr,'String'));
if val>0
  setappdata(gcf,'t_var',val);
else
  set(h_curr,'String',getappdata(gcf,'t_var'));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_t_num(h_curr,X,Y);
% Setting the number of data inputs to be added

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

val      = floor(str2num(get(h_curr,'String')));
if val>0
  setappdata(gcf,'t_num',val);
else
  set(h_curr,'String',getappdata(gcf,'t_num'));
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++

% --------------------------------------------------------------------
function [X, Y] = gui_set_likout(h_curr,X,Y);
% Setting the type of the assumed noise

% NET and EP are global variables, no need to be function arguments 
global net gpopt ep;

oldfn = func2str(net.likaddr);
val = get(h_curr,'Value');
% SWITCH because the EM-function needs to be set up!
switch val;
 case 1;		      % GAUSS
  net.likaddr  = @c_reg_gauss;
  net.likoptfn = @em_gauss;
 case 2;		      % LAPLACE
  net.likaddr  = @c_reg_lapl;
  net.likoptfn = @em_lapl;
 case 3;		      % POSEXP
  net.likaddr  = @c_reg_exp;
  net.likoptfn = @em_exp;
end;

if ~strcmp(upper(oldfn),upper(func2str(net.likaddr)));
  ogpreset;
end;
%%%+++++++++++++++++++++ END +++++++++++++++++++++++++++
