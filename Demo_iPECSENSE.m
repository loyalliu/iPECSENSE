%% Setup environment
clear all
close all
clc

%% Initialize parameters
PEC_GRAPPA_set_path;
%% Set parameters
PEC_GRAPPA_initial_para;

ncalib =128;
% ----- Directories -------------------------------------------------------
para.data_name = 'DTI_MB4';
para.data_dir = ['.\data\']; % Directory for the data
para.result_dir =['.' filesep fullfile('result', para.data_name) filesep]; % Directory for interim and final results
tmp_dir = fullfile(para.result_dir, 'tmp'); mkdir(tmp_dir);

% ----- VC-SAKE recon -----------------------------------------------------
para.VCSAKE.threshold_list = 0.01;
para.VCSAKE.nIter = 100;
para.VCSAKE.ksize = [1 1]*3;
para.VCSAKE.ncalib = ncalib;

% ----- CSM and PEM estimation --------------------------------------------
para.bart.cmd_for_PEM_calculation = {0.15, 6, ncalib};
para.bart.cmd_for_CSM_calculation = {0.03 , 6, ncalib};
para.recon.PEM_option.model='2dFreeESPIRiT';  %'2dFreeVcSake';; %'2dFreeVcSake';  % '2dFreeESPIRiT', '2dFreeBART', '1dLpc';
para.recon.CSM_option.source='segEPI SAKE';
%para.recon.CSM_option.source = 'gre';
%para.recon.CSM_option.source='segEPI LPC'

para.recon.slice2show = [1:48];
%para.recon.CSM_option.mask_type = 'CSM';
para.recon.CSM_option.mask_type = 'IMG';
para.recon.CSM_option.mask = 0.9;
para.recon.CSM_option.mask = 0.3;
para.recon.CSM_option.mask_dilate = 7;

para.CSM_calc.ncalib = ncalib;

% ----- iPEC-SENSE recon ------------------------------------
para.recon.filter = 4;
para.recon.reg_factor = 0.0001; % regularization factor for SENSE

para.recon.calib_size = [2 1]*ncalib;
para.recon.ksize = [4, 2];
para.recon.ksize = [6, 4];

% ----- Other information -------------------------------------------------
para.MB_factor = 4;
para.mon_size = [4 12];
para.skip.VCSAKE = 1; % skip VC-SAKE if already done
para.skip.recon = 0; % skip PEC-SENSE/PEC-GRAPPA/PEC-SP-SG reconstruction if already done

para.verbose.gmap = 1; % output g-factor map for evaluation
para.disp.gmap_range = [0.9 3];

para.disp.std_max = 60;

ndir=35; %number of directions
para.nIter = 200;
para.wavWeight = 0.000001;
para.prefix ='DTI_';
reps = [0 ndir];
para.rep_list = 1:ndir;

para_tmp = para;

%% Batch processing
% copyfile(para.data_dir, para.result_dir);
% PEC_GRAPPA_step2_lpcOnRefData(para);
% PEC_GRAPPA_step3_VCSAKE(para);
PEC_GRAPPA_step4_calc_CSM(para);

run_iPEC_SENSE = 1;

for mbf = 4
    para_tmp.MB_factor = mbf;
    para_tmp.rep_list = 1:ndir;

    % ----- iPEC-SENSE recon -----------------------------------------------  ----
    if run_iPEC_SENSE==1
        para = para_tmp; para.skip.PEM_calc = 1;
        para.recon.Recon_option = 'PEC-SENSE'; para.recon.CSM_option.filter='ESPIRiT';

        para.recon.flag_calc_tSNR = 0; SMS_EPI_Recon(para);  % Reconstruction for 1st frame
        para.recon.flag_calc_tSNR = 1; SMS_EPI_Recon(para);
    end

end
% end
%}