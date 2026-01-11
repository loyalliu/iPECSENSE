%%%%%  This function calculates CSM from different sources. 

function PEC_GRAPPA_step4_calc_CSM(para)
import LmyUtility.*

%% 
load_path = para.result_dir;
ncalib = para.CSM_calc.ncalib;
calib_size = [1 1]*ncalib;
msize = para.recon.CSM_option.msize;

% try % For Windows system
    switch para.recon.CSM_option.source
        case 'segEPI SAKE'
            load(fullfile(load_path, ls(fullfile(load_path,'SbData_sake_*'))),'epi_kxkyzc_sakeCor')
        case 'segEPI LPC'
            load(fullfile(load_path, ls(fullfile(load_path,'SbData_lpc.*'))),'epi_kxkyzc_lpcCor')
            epi_kxkyzc_sakeCor = epi_kxkyzc_lpcCor;
        case 'gre'
            load(fullfile(load_path, ls(fullfile(load_path,'GreData_128kx.*'))),'gre_kxkyzc')
            epi_kxkyzc_sakeCor = gre_kxkyzc;
 
    end

% catch % For Linux system
%     switch para.recon.CSM_option.source
%         case 'segEPI SAKE'
%             filename = fullfile(ls(fullfile(load_path,'SbData_sake_*')));
%             load(filename(1:end-1),'epi_kxkyzc_sakeCor');
%         case 'segEPI LPC'
%             filename = fullfile(ls(fullfile(load_path,'SbData_lpc_*')));
%             load(filename(1:end-1),'epi_kxkyzc_lpcCor')
%             epi_kxkyzc_sakeCor = epi_kxkyzc_lpcCor;
%         case 'gre'
%             filename = fullfile(ls(fullfile(load_path,'GreData_128kx.*')));
%             load(filename(1:end-1),'gre_kxkyzc')
%             epi_kxkyzc_sakeCor = gre_kxkyzc;
%     end
% end

% load(fullfile(load_path, ls(fullfile(load_path,'SbData_48kx*'))),'epi_kxkyzc_2shot_sakeCor')
% epi_kxkyzc_sakeCor = epi_kxkyzc_2shot_sakeCor;

epi_xyzc_sakeCor=ifft2c_MN(epi_kxkyzc_sakeCor, msize(1), msize(2));
[nFE,nPE,nSlice,nCoil]=size(epi_xyzc_sakeCor);
CSM_xyzc_epi_sakeCor=zeros(nFE,nPE,nSlice,nCoil);
CSMweight_xyz_epi_sakeCor=zeros(nFE,nPE,nSlice);

img_ref = sos(ifft2c(epi_kxkyzc_sakeCor));

%% using BART
bart_setting = para.bart.cmd_for_CSM_calculation;
for iSlice=1:nSlice
    if ~iscell(bart_setting)
        [a, b] = bart(bart_setting,fft2c(epi_xyzc_sakeCor(:,:,iSlice,:)));
        CSM_xyzc_epi_sakeCor(:,:,iSlice,:) = a(:,:,1,:,1);
        CSMweight_xyz_epi_sakeCor(:,:,iSlice) = b(:,:,1,1,1);
    else
        [a, b] = gp_Calculate_CSM(squeeze(epi_xyzc_sakeCor(:,:,iSlice,:)), 'ESPIRiT', bart_setting);
        CSM_xyzc_epi_sakeCor(:,:,iSlice,:) = a;
        CSMweight_xyz_epi_sakeCor(:,:,iSlice) = b(:,:,end);
%         imshow(b(:, :, end), [0 1]); pause;
        disp(['Calculate CSM for Slice #' num2str(iSlice)]);
    end
    
end

% Save computed CSM
switch para.recon.CSM_option.source
    case 'segEPI SAKE'
        CSM_xyzc_epi_sakeCor=single(CSM_xyzc_epi_sakeCor);
        CSMweight_xyz_epi_sakeCor=single(CSMweight_xyz_epi_sakeCor);
        save(fullfile(load_path,['CSM_xyzc_epi_sakeCor_bart']),'CSM_xyzc_epi_sakeCor','CSMweight_xyz_epi_sakeCor', 'img_ref', 'calib_size','-v7.3')%
    case 'segEPI LPC'
        CSM_xyzc_epi_lpcCor=single(CSM_xyzc_epi_sakeCor);
        CSMweight_xyz_epi_lpcCor=single(CSMweight_xyz_epi_sakeCor);
        save(fullfile(load_path,['CSM_xyzc_epi_lpcCor_bart']),'CSM_xyzc_epi_lpcCor','CSMweight_xyz_epi_lpcCor', 'img_ref', 'calib_size','-v7.3')%
    case 'gre'
        CSM_xyzc_gre=single(CSM_xyzc_epi_sakeCor);
        CSMweight_xyz_gre=single(CSMweight_xyz_epi_sakeCor);
        save(fullfile(load_path,['CSM_xyzc_gre']),'CSM_xyzc_gre','CSMweight_xyz_gre', 'img_ref', 'calib_size','-v7.3')%
end

end % end of the function