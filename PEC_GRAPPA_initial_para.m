para.data_name = 'phantom_1';
para.data_dir = ['.\sample_data\' para.data_name];
para.result_dir = ['.\results\' para.data_name];
para.nshot = 2;
para.mon_size = [5 16];
para.MB_factor = 2;
para.rep_list = 1:2;

para.crop_ky_range = 1:128; % used in LPC correction for ref data

% For VC-SAKE
para.VCSAKE.p = 0.1;
para.VCSAKE.ncalib = 48;
para.VCSAKE.ksize = [3, 3];
para.VCSAKE.nIter = 100;
para.VCSAKE.verbose_ind_coil = [];
% para.VCSAKE.verbose_ind_coil = [11 21];
para.VCSAKE.threshold_list=[4*ones(1,10),4.5*ones(1,15),5*ones(1,30),4.5*ones(1,15),4*ones(1,10)];

para.CSM_calc.ncalib = 48;

% For final reconstruction
para.recon.slice2show = [1:80];
para.recon.reg_factor = 0.001; % regularization factor for both SENSE and GRAPPA
para.recon.ksize = [4, 4]; % For GRAPPA reconstruction
para.recon.calib_size = [120, 60]; % For GRAPPA reconstruction
para.recon.CSM_option.source='segEPI SAKE';
para.recon.CSM_option.filter='ESPIRiT';
para.recon.CSM_option.mask_type = 'IMG';
para.recon.CSM_option.mask = 0.7; % Threshold for masking, 0 refers to no masking
para.recon.CSM_option.msize = [128 130];

para.recon.PEM_option.model='2dFreeBART';
para.recon.PEM_option.ref='ext';
para.recon.PEM_option.csm='sake';
para.recon.PEM_size = [128 128];
para.recon.SMS_EPI_Get_DataAndParmeters = @SMS_EPI_Get_DataAndParmeters;
para.recon.Recon_option = 'PEC-SENSE';
para.recon.nCoilReduced = 64;
para.recon.flag_calc_tSNR = 0;
para.recon.save_unmixing_map = 0;

% Skip some steps?
para.skip.VCSAKE = 0;
para.skip.recon = 0;
para.skip.PEM_calc = 0;
para.skip.regen_pic = 0;
% For verbose output ...
para.verbose.LPC = 1;
para.verbose.VCSAKE = 1;
para.verbose.recon = 1;
para.verbose.cmrr = 0; % 1 for outputing default CMRR reconstruction

% For display
para.disp.SNR_max = 50;
para.disp.PEM_range = [-1 1]*pi;
para.disp.ROT_flag = 1;
para.disp.img_surfix = [];
para.disp.mf = '-m2'; % magnify ...

% For BART
para.bart.cmd_for_PEM_calculation = 'ecalib -c0 -r 48 -k 9 -t0.005'; % for oblique phantom ...
para.bart.cmd_for_CSM_calculation = 'ecalib -c0 -m 1 -r 48';

para.slice_info.ns = 48;
para.slice_info.s2x2 = [1:20 41:60];
para.slice_info.s2x4 = [1:10 41:50];
para.slice_info.s2x5 = [1:8 41:48];