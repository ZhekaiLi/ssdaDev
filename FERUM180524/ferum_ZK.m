%     Finite Element Reliability Using Matlab Modified by ZK, FERUM_ZK, Version 1.0, 2020 
%
%     Compare to the initial ferum.m 
%     ferum_ZK.m retain
%        case 0: Exit
%        case 1: Help
%        case 23: Subset Simulation Analysis
%        case 24: Subset Simulation Analysis w/ delayed acceptance (inital SSDA)
%     add
%        case 25: ssda + store GP model in each run to find the prediction with the least variance
%        case 26: ssda + use Matlab RGP to replace the initial ogp method
%        case 27: ssda modified  



if exist('analysisopt','var')

   if isfield(analysisopt,'echo_flag')
      echo_flag = analysisopt.echo_flag;
   else
      echo_flag = 1;
   end
   
end


if echo_flag

   clc, format short e

   disp('   ________________________________________________________________________');
   disp('  | Welcome to FERUM_ZK Version 1.0 (Finite Element Reliability Using Matlab Modified by ZK) |');
   disp('  | Note: All the analysis options below assumes that necessary data       |');
   disp('  |       are available in the current Matlab workspace.                   |');
   disp('   ________________________________________________________________________');
   disp(' ');
   disp('   0: Exit');
   disp('   1: Help');
   disp('  23: Subset Simulation Analysis');
   disp('  24: Subset Simulation Analysis w/ delayed acceptance');
   disp('  25: ssda + store GP model in each run to find the prediction with the least variance');
   disp('  26: ssda + use Matlab RGP to replace the initial ogp method');
   disp('  27: ssda modified');
   disp(' ');
   analysistype = input('  CHOOSE OPTION FROM THE LIST ABOVE: ');
   analysisopt.analysistype = analysistype;

end


if analysisopt.analysistype > 2
   
   if ~isfield(analysisopt,'already_updated')
      % This function updates probdata and gfundata before any analysis (must be run only once)
      if echo_flag
         [probdata,gfundata,analysisopt] = update_data(1,probdata,analysisopt,gfundata,femodel);
         dummy=input('Hit a key to continue\n');
      else
         [probdata,gfundata,analysisopt] = update_data(1,probdata,analysisopt,gfundata,femodel);
      end
   end
   
   % This function completely determines and updates parameters, mean and standard deviation 
   % associated with the distribution of each random variable
   probdata.marg  = distribution_parameter(probdata.marg);
   
end


tic


switch analysisopt.analysistype
   

   case 0 % ---- EXIT ---------------------------------------------------------------------------------------------------
      
      disp(' ');
      disp('  Bye, bye.');
      disp(' ');
     
   
   case 1 % ---- HELP ---------------------------------------------------------------------------------------------------

      clc
      disp(' ');
      disp('  FERUM HELP');
      disp(' ');
      disp('  If you are new to FERUM, you are recommended to visit the web page:');
      disp('  http://www.ifma.fr/FERUM');
      disp(' ');
      disp('  To run FERUM, do the following:');
      disp('  1. Specify necessary parameters in your current Matlab workspace. ');
      disp('     (The format of the input can be found in the inputfile_template.m file.');
      disp('     If you want to try one of the provided example inputfiles,');
      disp('     simply read the file into your workspace by writing the file name');
      disp('     without .m extension and press enter. ');
      disp('  2. Start FERUM (the shell program) by issuing the command ferum.');
      disp('  3. Choose the alternative that fits your purpose.');
      disp('  4. Main results as displayed on the screen and more detailed results are usually');
      disp('     stored in a specific data structure (e.g. formresults for a FORM analysis).');
      disp(' ');
   
   
      
   
  
      
      
   case 23 % ---- SUBSET SIMULATIONS ------------------------------------------------------------------------------------

      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('SUBSET SIMULATION is running, please wait... (Ctrl+C breaks)')
         disp([' '])
      end

      % Run simulation analysis
      [ subsetsimulationresults, probdata ] = subset_simulation(1,probdata,analysisopt,gfundata,femodel,randomfield);
   
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING SUBSET SIMULATION' ])
         disp([' '])
         disp('Probability of failure: '), disp(subsetsimulationresults.pf)
         if isfield(subsetsimulationresults,'cov_pf')
            disp('Coefficient of variation of failure probability (lower bound):'), disp(subsetsimulationresults.cov_pf)
         end
         disp('Reliability index beta: '), disp(subsetsimulationresults.beta)
         disp('Number of calls to the limit-state function: '), disp(subsetsimulationresults.nfun)
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   subsetsimulationresults.pf         = Failure probability from subset simulations')
         if isfield(subsetsimulationresults,'cov_pf')
            disp('   subsetsimulationresults.cov_pf     = Coefficient of variation for the failure probability (lower bound)')
         end
         disp('   subsetsimulationresults.beta       = Generalized reliability index beta from this simulation')
         disp('   subsetsimulationresults.SubsetData = Subset data structure')
         disp('   subsetsimulationresults.nfun       = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
      
   case 24 % ---- SUBSET SIMULATIONS W/ DELAYED ACCEPTANCE---------------------------------------------------------------

      if echo_flag
         % Clear screen and display message
         disp(' ');
         disp('SUBSET SIMULATION W/ DELYED ACCEPTANCE is running, please wait... (Ctrl+C breaks)')
         disp(' ');
      end

      % Run simulation analysis
      [ssda_result, probdata ] = ss_da(1,probdata,analysisopt,gfundata,femodel,randomfield);
   
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING SUBSET SIMULATION W/ DELAYED ACCEPTANCE' ])
         disp([' '])
         disp('Probability of failure: '), disp(ssda_result.pf)
         if isfield(ssda_result,'cov_pf')
            disp('Coefficient of variation of failure probability (lower bound):'), disp(ssda_result.cov_pf)
         end
         disp('Reliability index beta: '), disp(ssda_result.beta)
         disp('Number of calls to the limit-state function: '), disp(ssda_result.nfun)
         disp('Accumulated GP accept rate'), disp(sum(ssda_result.ssda_Data.AccRate(2,:)));
         disp('Accumulated SS accept rate'), disp(sum(ssda_result.ssda_Data.AccRate(3,:)));
         disp('GP rate / SS rate'), disp(sum(ssda_result.ssda_Data.AccRate(2,:)) / sum(ssda_result.ssda_Data.AccRate(3,:)));
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   ssda_result.pf         = Failure probability from subset simulations')
         if isfield(ssda_result,'cov_pf')
            disp('   ssda_result.cov_pf     = Coefficient of variation for the failure probability (lower bound, upper bound)')
         end
         disp('   ssda_result.beta       = Generalized reliability index beta from this simulation')
         disp('   ssda_result.SubsetData = Subset data structure')
         disp('   ssda_result.nfun       = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
      
   case 25 % ---- ssda + store GP model in each run to find the prediction with the least variance ----------------------

      global NETs;
      NETs = {};
      if echo_flag
         % Clear screen and display message
         disp(' ');
         disp('SUBSET SIMULATION W/ DELYED ACCEPTANCE is running, please wait... (Ctrl+C breaks)')
         disp(' ');
      end

      % Run simulation analysis
      [ssda_result, probdata ] = ss_da_compareGPMs(1,probdata,analysisopt,gfundata,femodel,randomfield);
   
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING SUBSET SIMULATION W/ DELAYED ACCEPTANCE' ])
         disp([' '])
         disp('Probability of failure: '), disp(ssda_result.pf)
         if isfield(ssda_result,'cov_pf')
            disp('Coefficient of variation of failure probability (lower bound):'), disp(ssda_result.cov_pf)
         end
         disp('Reliability index beta: '), disp(ssda_result.beta)
         disp('Number of calls to the limit-state function: '), disp(ssda_result.nfun)
         disp('Reliability index beta: '), disp(ssda_result.beta)
         disp('Number of calls to the limit-state function: '), disp(ssda_result.nfun)
         disp('Accumulated GP accept rate'), disp(sum(ssda_result.ssda_Data.AccRate(2,:)));
         disp('Accumulated SS accept rate'), disp(sum(ssda_result.ssda_Data.AccRate(3,:)));
         disp('GP rate / SS rate'), disp(sum(ssda_result.ssda_Data.AccRate(2,:)) / sum(ssda_result.ssda_Data.AccRate(3,:)));
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   ssda_result.pf         = Failure probability from subset simulations')
         if isfield(ssda_result,'cov_pf')
            disp('   ssda_result.cov_pf     = Coefficient of variation for the failure probability (lower bound, upper bound)')
         end
         disp('   ssda_result.beta       = Generalized reliability index beta from this simulation')
         disp('   ssda_result.SubsetData = Subset data structure')
         disp('   ssda_result.nfun       = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
   
   case 26 % ---- ssda + use Matlab RGP to replace the initial ogp method -----------------------------------------------
      global data;
      global gpm;

      if echo_flag
         % Clear screen and display message
         disp(' ');
         disp('SUBSET SIMULATION W/ DELYED ACCEPTANCE is running, please wait... (Ctrl+C breaks)')
         disp(' ');
      end

      % Run simulation analysis
      [ssda_result, probdata ] = ss_da_GPinMatlab(1,probdata,analysisopt,gfundata,femodel,randomfield);
   
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING SUBSET SIMULATION W/ DELAYED ACCEPTANCE' ])
         disp([' '])
         disp('Probability of failure: '), disp(ssda_result.pf)
         if isfield(ssda_result,'cov_pf')
            disp('Coefficient of variation of failure probability (lower bound):'), disp(ssda_result.cov_pf)
         end
         disp('Reliability index beta: '), disp(ssda_result.beta)
         disp('Number of calls to the limit-state function: '), disp(ssda_result.nfun)
         disp('Reliability index beta: '), disp(ssda_result.beta)
         disp('Number of calls to the limit-state function: '), disp(ssda_result.nfun)
         disp('Accumulated GP accept rate'), disp(sum(ssda_result.ssda_Data.AccRate(2,:)));
         disp('Accumulated SS accept rate'), disp(sum(ssda_result.ssda_Data.AccRate(3,:)));
         disp('GP rate / SS rate'), disp(sum(ssda_result.ssda_Data.AccRate(2,:)) / sum(ssda_result.ssda_Data.AccRate(3,:)));
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   ssda_result.pf         = Failure probability from subset simulations')
         if isfield(ssda_result,'cov_pf')
            disp('   ssda_result.cov_pf     = Coefficient of variation for the failure probability (lower bound, upper bound)')
         end
         disp('   ssda_result.beta       = Generalized reliability index beta from this simulation')
         disp('   ssda_result.SubsetData = Subset data structure')
         disp('   ssda_result.nfun       = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
  
   case 27 % ---- ssda modified -----------------------------------------------------------------------------------------

      global NETs;
      NETs = {};
      if echo_flag
         % Clear screen and display message
         disp(' ');
         disp('SUBSET SIMULATION W/ DELYED ACCEPTANCE is running, please wait... (Ctrl+C breaks)')
         disp(' ');
      end

      % Run simulation analysis
      [ssda_result, probdata ] = ss_da_modified(1,probdata,analysisopt,gfundata,femodel,randomfield);
   
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING SUBSET SIMULATION W/ DELAYED ACCEPTANCE' ])
         disp([' '])
         disp('Probability of failure: '), disp(ssda_result.pf)
         if isfield(ssda_result,'cov_pf')
            disp('Coefficient of variation of failure probability (lower bound):'), disp(ssda_result.cov_pf)
         end
         disp('Reliability index beta: '), disp(ssda_result.beta)
         disp('Number of calls to the limit-state function: '), disp(ssda_result.nfun)
         disp('Reliability index beta: '), disp(ssda_result.beta)
         disp('Number of calls to the limit-state function: '), disp(ssda_result.nfun)
         disp('Accumulated GP accept rate'), disp(sum(ssda_result.ssda_Data.AccRate(2,:)));
         disp('Accumulated SS accept rate'), disp(sum(ssda_result.ssda_Data.AccRate(3,:)));
         disp('GP rate / SS rate'), disp(sum(ssda_result.ssda_Data.AccRate(2,:)) / sum(ssda_result.ssda_Data.AccRate(3,:)));
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   ssda_result.pf         = Failure probability from subset simulations')
         if isfield(ssda_result,'cov_pf')
            disp('   ssda_result.cov_pf     = Coefficient of variation for the failure probability (lower bound, upper bound)')
         end
         disp('   ssda_result.beta       = Generalized reliability index beta from this simulation')
         disp('   ssda_result.SubsetData = Subset data structure')
         disp('   ssda_result.nfun       = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
   
   
   otherwise % --------------------------------------------------------------------------------------
      
      disp(' ');
      disp('  An invalid choice was entered.');
      disp(' ');
  
end % --------------------------------------------------------------------------------------


if echo_flag
   toc
end
