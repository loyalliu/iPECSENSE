% 1D LPC correction for reference EPI data
function PEC_GRAPPA_step2_lpcOnRefData(para)

import LmyGhostCorrection.*
import LmyUtility.*
load_folder = para.result_dir;
if isfield(para, 'refscan') && para.refscan.flag == 1
    load(fullfile(load_folder,'SbData'),'epi_sampling_mask_kys','epi_fovY','epi_rangeKy');
    load(fullfile(load_folder, para.refscan.filename),'MB_refscan');
    epi_kxkyzc = MB_refscan;
else
    load(fullfile(load_folder,'SbData'),'epi_kxkyzc','epi_sampling_mask_kys','epi_fovY','epi_rangeKy')
end

nShot = para.nshot;
ky_range = para.crop_ky_range;

if max(ky_range) < size(epi_kxkyzc, 2)
    epi_kxkyzc=single(epi_kxkyzc(:,ky_range,:,:));
    epi_rangeKy = length(ky_range);
    disp('cropped');
else
    if max(ky_range)-10 > size(epi_kxkyzc, 2)
        epi_rangeKy = size(epi_kxkyzc, 2);
        disp('unchanged');
    else
        [nx, ny, ns, nc] = size(epi_kxkyzc);
        tmp = zeros([nx, max(ky_range), ns, nc]);
        tmp(:, 1:ny, :, :) = epi_kxkyzc;
        epi_kxkyzc = tmp;
        epi_rangeKy = length(ky_range);
        disp('zeropadded');
    end
end

LPC_method = 1;
if ~exist(fullfile(load_folder,'SbData_lpc.mat'), 'file')
% if 1 % ~exist(fullfile(load_folder,'SbData_lpc.mat'), 'file')
    switch LPC_method
        case {'Entropy_min', 1}
            profile = sos(epi_kxkyzc(:, :, :));
            profile = sos(profile');
            ind_kys = find(profile>0);
            if mod(length(ind_kys), 2) == 1
                ind_kys = [ind_kys; ind_kys(end)+nShot];
            end
            % Detect if the single band data is undersampled
            if sum(profile==0) > 5
                tmp = epi_kxkyzc(:, ind_kys, :, :);
                % No significant displacement found in Siemens data
%                 [tmp,even_from_odd] = oneDimlinearCorr_bulk(tmp,1); 
                [tmp, phasepara_evenOddMajar_zp] = oneDimLinearCorr_entropy(tmp, 1);
                epi_kxkyzc_LPC = epi_kxkyzc;
                epi_kxkyzc_LPC(:, ind_kys, :, :) = tmp;
            else
                [epi_kxkyzc_LPC, phasepara_evenOddMajar_zp] = oneDimLinearCorr_entropy(epi_kxkyzc, nShot);          % LPC for even/odd echo phase error
            end
            if isfield(para, 'refscan') && para.refscan.flag == 1 && isfield(para.refscan, 'phasepara_mapping')
                phasepara_mapping = para.refscan.phasepara_mapping;
                phasepara_evenOddMajar_zp(:, 1) = phasepara_evenOddMajar_zp(:, 1)*phasepara_mapping(1, 1)+phasepara_mapping(1, 2);
                phasepara_evenOddMajar_zp(:, 2) = phasepara_evenOddMajar_zp(:, 2)*phasepara_mapping(2, 1)+phasepara_mapping(2, 2);
            end
            if isfield(para, 'refscan') && para.refscan.flag == 1 && isfield(para.refscan, 'phasepara')
                phasepara_evenOddMajar_zp = para.refscan.phasepara;
            end
            
            %
            if sum(profile==0) < 5
                if nShot > 1 % For multishot cases, inter-shot phase error needs to be corrected
                    [epi_kxkyzc_CPC, phasepara_interShot_zp] = oneDimConstCorr_entropy(epi_kxkyzc_LPC, nShot); % LPC for inter-shot phase error
                    [epi_kxkyzc_LPC1, evenOddResidual_zp] = oneDimLinearCorr_entropy(epi_kxkyzc_CPC, nShot);    % LPC for even/odd echo phase error
                    epi_kxkyzc_lpcCor=single(epi_kxkyzc_LPC1);
                    %
                    save(fullfile(load_folder,'SbData_lpc'),'epi_kxkyzc_lpcCor', 'epi_kxkyzc_LPC', 'epi_kxkyzc_LPC1', ...
                        'epi_kxkyzc_CPC', 'phasepara_*','epi_fovY','epi_rangeKy','-v7.3');
                else
                    epi_kxkyzc_lpcCor=single(epi_kxkyzc_LPC);
                    save(fullfile(load_folder,'SbData_lpc'),'epi_kxkyzc_lpcCor','phasepara_*','epi_fovY','epi_rangeKy','-v7.3'); % Save raw data
                end
            else
                epi_kxkyzc_CPC = epi_kxkyzc_LPC;
                epi_kxkyzc_LPC1 = epi_kxkyzc_LPC;
                epi_kxkyzc_lpcCor = epi_kxkyzc_LPC;
                save(fullfile(load_folder,'SbData_lpc'),'epi_kxkyzc_lpcCor','phasepara_*','epi_fovY','epi_rangeKy','-v7.3'); % Save raw data
            end
            %}
        case 'SLR_min'
            
        case {'LR', 3}
%             kSize = [1 1]*6; wnRank = 0.5;
            kSize = [1 1]*3; wnRank = 0.7;
%             kSize = [1 1]*3; wnRank = 1.5;
            [nx, ny, ns, nc] = size(epi_kxkyzc);
            epi_kxkyzc_LPC = epi_kxkyzc;
            x_trend = {};
            for ind_s = 1:ns
                tmp0 = squeeze(epi_kxkyzc_LPC(:, :, ind_s, :));
                itmp0 = ifftc(tmp0, 1);
                tmp = tmp0;
                tic
                switch ind_s
                    case {4, 8, 12, 16}
                        niter = 100;
                    otherwise
                        niter = 100;
                end
                for ind_iter = 1:niter
                     % reorder data to get Hankel structure. 
                    tmp = im2row(tmp, kSize); 
                    [tsx,tsy,tsz] = size(tmp);
                    A = reshape(tmp,tsx,tsy*tsz);
                    % SVD thresholding
                    R = floor(wnRank*prod(kSize));
%                     [U,S,V] = gp_givefastSVD(A, R, 0);
                    [U,S,V] = svd(A,'econ');
                    keep = 1:R;
                    A = U(:,keep)*S(keep,keep)*V(:,keep)';
                    % Enforce Hankel structure
                    A = reshape(A,tsx,tsy,tsz);
                    tmp = row2im(A,[nx,ny,nc],kSize);
                    
                    itmp = ifftc(tmp, 1);
                    itmp = itmp.*conj(itmp0);
                    for ind_seg = 1:2*nShot
                        itmp1 = itmp(:, ind_seg:2*nShot:end, :);
                        itmp1 = itmp1(:, :);
                        itmp1 = sum(itmp1, 2);
%                         itmp1 = itmp1.';
                        tmpa = angle(itmp1);
                        tmpm = abs(itmp1);
                        Aa = [ones(nx, 1), (1:nx)'];
%                         Aa = [ones(nx, 1), (1:nx)' (1:nx)'.^2];
                        x = lscov(Aa, tmpa, tmpm);
%                         disp(x'*1000);
                        tmpa_fitted = Aa*x;
                        
%                         plot(tmpa, 'LineWidth', 3);  hold on;
%                         plot(tmpa_fitted, 'r-', 'LineWidth', 3); hold off;
%                         legend('Measured', 'Linearly fitted');
%                         xlabel('Pixel');
%                         ylabel('Phase Difference (Rad)');
%                         xlim([1 length(tmpa)])
%                         set(gca,'FontSize',20)
%                         
%                         title(ind_seg);
%                         pause(0.1);
%                         pause;
%                         
                        itmp(:, ind_seg:2*nShot:end, :) = ...
                            itmp0(:, ind_seg:2*nShot:end, :).*...
                            repmat(exp(1j*tmpa_fitted), [1 ny/2/nShot nc]);
                        
                        if ind_iter == 1
                            x_trend{ind_s}(ind_iter, :, ind_seg) = x';
                        else
                            x_trend{ind_s}(ind_iter, :, ind_seg) = x_trend{ind_s}(ind_iter-1, :, ind_seg) + x';
                        end
                    end 
                    
%                     disp(' ')
                    itmp0 = itmp;
                    tmp = fftc(itmp, 1);
                    
                    if max(x) < 1e-4
                        disp(['Iteration #' num2str(ind_iter)]);
                        
%                         S = diag(S); % S = S(1:2*prod(kSize));
%                         S = S/max(S);
%                         S = log(S);
%                         St = S; St(10*prod(kSize):end) = min(S)-10; 
%                         plot((1:length(S))/prod(kSize), St, 'r', 'LineWidth', 3); hold on;
%                         plot((1:length(S))/prod(kSize), S, 'b--', 'LineWidth', 3);
%                         xlabel('Window normalized singular vector number');
%                         ylabel('Singular values in log scale');
%                         xlim([0 length(S)/prod(kSize)])
%                         ylim([min(S) max(S)]);
%                         set(gca,'FontSize',20)
%                         pause;
                        
                       break; 
                    end
                    
                end
                toc
%                 disp(num2str(ind_s))
%                 pause;
                epi_kxkyzc_LPC(:, :, ind_s, :) = reshape(tmp, [nx, ny, 1, nc]);
            end
            
            save LR_x_trend.mat x_trend
            epi_kxkyzc_CPC = epi_kxkyzc_LPC;
            epi_kxkyzc_LPC1 = epi_kxkyzc_LPC;
    end
    
else
    % Load existing data
    load(fullfile(load_folder,'SbData_lpc'));

end
%% Exporting LPC corrected images
if para.verbose.LPC
    export_dir = [fullfile(load_folder, 'LPC') filesep]; mkdir(export_dir);
    para.mon_size = [4 12];
    
    ns = size(epi_kxkyzc, 3);
    slice_inds = 1:ns;
    
    epi_kxkyzc = rot90(epi_kxkyzc(:, :, slice_inds, :), 1);
    epi_kxkyzc_CPC = rot90(epi_kxkyzc_CPC(:, :, slice_inds, :), 1);
    epi_kxkyzc_LPC = rot90(epi_kxkyzc_LPC(:, :, slice_inds, :), 1);
    epi_kxkyzc_LPC1 = rot90(epi_kxkyzc_LPC1(:, :, slice_inds, :), 1);
    para.mon_size = [4 12];

% ss = 10;
% epi_kxkyzc = rot90(epi_kxkyzc(:, :, ss, :), 1);
% epi_kxkyzc_CPC = rot90(epi_kxkyzc_CPC(:, :, ss, :), 1);
% epi_kxkyzc_LPC = rot90(epi_kxkyzc_LPC(:, :, ss, :), 1);
% epi_kxkyzc_LPC1 = rot90(epi_kxkyzc_LPC1(:, :, ss, :), 1);
% para.mon_size = [1 1];
% 
% para.mon_size = [3 6];
    
    % ----- Display with original brightness ------------------------------
    mon_size = para.mon_size;
    temp=sos(ifft2c(epi_kxkyzc));
    disp_max999=prctile(abs(temp(:)),99.9);
    
    im2show = sos(ifft2c(epi_kxkyzc));
    figure(1); MY_montage(im2show,'size',mon_size,'displayrange',[0 disp_max999]); % title('before correction');
    export_fig([export_dir 'LPC_step0_No_Correction'],'-m2','-tif');
    im2show = sos(ifft2c(epi_kxkyzc_LPC));
    figure(1); MY_montage(im2show,'size',mon_size,'displayrange',[0 disp_max999]); % title('Step 1: correct as 2 shots');
    export_fig([export_dir 'LPC_step1'],'-m2','-tif');
    if nShot > 1
        im2show = sos(ifft2c(epi_kxkyzc_CPC));
        figure(1); MY_montage(im2show,'size',mon_size,'displayrange',[0 disp_max999]); % title('Step 2: correct as 1 shot');
        export_fig([export_dir 'LPC_step2'],'-m2','-tif');
        im2show = sos(ifft2c(epi_kxkyzc_LPC1));
        figure(1); MY_montage(im2show,'size',mon_size,'displayrange',[0 disp_max999]); % title('Step 3: correct as 2 shots again');
        export_fig([export_dir 'LPC_step3'],'-m2','-tif');
    end
    
    % ----- Display with enhanced brightness ------------------------------
    disp_max999=prctile(abs(temp(:)),99.9)/5;
    
    im2show = sos(ifft2c(epi_kxkyzc));
    figure(1); MY_montage(im2show,'size',mon_size,'displayrange',[0 disp_max999]); % title('before correction');
    export_fig([export_dir 'LPC_x10_step0_No_Correction'],'-m2','-tif');
    im2show = sos(ifft2c(epi_kxkyzc_LPC));
    figure(1); MY_montage(im2show,'size',mon_size,'displayrange',[0 disp_max999]); % title('Step 1: correct as 2 shots');
    export_fig([export_dir 'LPC_x10_step1'],'-m2','-tif');
    if nShot > 1
        im2show = sos(ifft2c(epi_kxkyzc_CPC));
        figure(1); MY_montage(im2show,'size',mon_size,'displayrange',[0 disp_max999]); % title('Step 2: correct as 1 shot');
        export_fig([export_dir 'LPC_x10_step2'],'-m2','-tif');
        im2show = sos(ifft2c(epi_kxkyzc_LPC1));
        figure(1); MY_montage(im2show,'size',mon_size,'displayrange',[0 disp_max999]); % title('Step 3: correct as 2 shots again');
        export_fig([export_dir 'LPC_x10_step3'],'-m2','-tif');
    end
    
end % Ending for exporting LPC corrected images

end