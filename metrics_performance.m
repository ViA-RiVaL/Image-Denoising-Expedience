clc
close all
clear
%%
addpath(genpath('FR IQA')); % Full-reference metrics
addpath(genpath('Dependencies'));
%%
imNameRef = [1 3 5 7 8 17 19 23 25:32];   % Reference images 
sigma = [3 5 10 15 20 25 30]; % noise std
filters = {'DCTF','BM3D'}; 
metrics_gray = {'PSNR','PSNR-HVS','PSNR-HVS-M','PSNR-HMA','PSNR-HA','MAD index', ...
    'VIF','IFC','UQI','SSIM','MS-SSIM','FSIM','IW-PSNR','IW-SSIM','DSI','GMSD',...
    'HaarPSI','RFSIM','SR-SIM','GSM','NQM','WSNR','VIFp','ADD-SSIM','ADD-GSIM',...
    'IGM','CW-SSIM','ADM','DSS','M-SVD','QILV','IISI','CVSSI','WASH','ESSIM',...
    'MCSD'};

metrics_color = {'VSI','SFF','PSIM','MDSI','IFS','UNIQUE','MS-UNIQUE'};
%%
folder_ref = 'images\ref';
folder_n = 'images\noisy';
folder_f = 'images\filtered';
%%
score = zeros(numel(imNameRef), numel(sigma));
%% Noisy images processing
for z = 1:1:numel(metrics_gray)
    for i = 1:1:numel(imNameRef)
        % read a reference image
        img_ref = double(imread([pwd '\' folder_ref '\' ...
            num2str(imNameRef(i)) '.bmp']));
        % read a noisy image
        for j = 1:1:numel(sigma)
            
            img_n = double(imread([pwd '\' folder_n '\' ...
                'NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
            
            score(i,j) = IQA_score(img_ref, img_n, metrics_gray{z});       

        end
    end
    
    clc
    disp(['Metrics: ' num2str(100*z/numel(metrics_gray)) '%...'])
    
    save(strcat(metrics_gray{z},'_n'), 'score');
end
%% color
for z = 1:1:numel(metrics_color)
    for i = 1:1:numel(imNameRef)
        % read a reference image
        img_ref = double(imread([pwd '\' folder_ref '\' ...
            num2str(imNameRef(i)) '.bmp']));
        % read a noisy image
        img_ref = repmat(img_ref, [1 1 3]);
        for j = 1:1:numel(sigma)
            
            img_n = double(imread([pwd '\' folder_n '\' ...
                'NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
            
            img_n = repmat(img_n, [1 1 3]);
            
            score(i,j) = IQA_score(img_ref, img_n, metrics_color{z});       

        end
    end
    
    clc
    disp(['Metrics: ' num2str(100*z/numel(metrics_color)) '%...'])
    
    save(strcat(metrics_color{z},'_n'), 'score');
end
%% Filtered images processing DCTF
for z = 1:1:numel(metrics_gray)
    for i = 1:1:numel(imNameRef)
        % read a reference image
        img_ref = double(imread([pwd '\' folder_ref '\' ...
            num2str(imNameRef(i)) '.bmp']));
        % read a filtered image
        for j = 1:1:numel(sigma)
            

            img_f = double(imread([pwd '\' folder_f '\' ...
                'DCTF_NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
            
            score(i,j) = IQA_score(img_ref, img_f, metrics_gray{z});       

        end
    end
    
    clc
    disp(['Metrics: ' num2str(100*z/numel(metrics_gray)) '%...'])
    
    save(strcat(metrics_gray{z},'_DCTF_f'), 'score');
end
%% color DCTF
for z = 1:1:numel(metrics_color)
    for i = 1:1:numel(imNameRef)
        % read a reference image
        img_ref = double(imread([pwd '\' folder_ref '\' ...
            num2str(imNameRef(i)) '.bmp']));
        img_ref = repmat(img_ref, [1 1 3]);
        % read a filtered image
        
        for j = 1:1:numel(sigma)
            

            img_f = double(imread([pwd '\' folder_f '\' ...
                'DCTF_NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
            
            img_f = repmat(img_f, [1 1 3]);
            
            score(i,j) = IQA_score(img_ref, img_f, metrics_color{z});       

        end
    end
    
    clc
    disp(['Metrics: ' num2str(100*z/numel(metrics_color)) '%...'])
    
    save(strcat(metrics_color{z},'_DCTF_f'), 'score');
end
%%
%% color BM3D
for z = 1:1:numel(metrics_color)
    for i = 1:1:numel(imNameRef)
        % read a reference image
        img_ref = double(imread([pwd '\' folder_ref '\' ...
            num2str(imNameRef(i)) '.bmp']));
        img_ref = repmat(img_ref, [1 1 3]);
        % read a filtered image
        
        for j = 1:1:numel(sigma)
            

            img_f = double(imread([pwd '\' folder_f '\' ...
                'BM3D_NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
            
            img_f = repmat(img_f, [1 1 3]);
            
            score(i,j) = IQA_score(img_ref, img_f, metrics_color{z});       

        end
    end
    
    clc
    disp(['Metrics: ' num2str(100*z/numel(metrics_color)) '%...'])
    
    save(strcat(metrics_color{z},'_BM3D_f'), 'score');
end
%%  Filtered images processing BM3D
for z = 1:1:numel(metrics_gray)
    for i = 1:1:numel(imNameRef)
        % read a reference image
        img_ref = double(imread([pwd '\' folder_ref '\' ...
            num2str(imNameRef(i)) '.bmp']));
        % read a filtered image
        for j = 1:1:numel(sigma)
            

            img_f = double(imread([pwd '\' folder_f '\' ...
                'BM3D_NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
            
            score(i,j) = IQA_score(img_ref, img_f, metrics_gray{z});       

        end
    end
    
    clc
    disp(['Metrics: ' num2str(100*z/numel(metrics_gray)) '%...'])
    
    save(strcat(metrics_gray{z},'_BM3D_f'), 'score');
end
%%
%%
% iqa_n = 'IQA Results\n';
% iqa_dctf = 'IQA Results\f_dctf';
% iqa_bm3d = 'IQA Results\f_bm3d';
%%
% metrics = [metrics_gray, metrics_color];
%%
%%
%%
matrics_add = {'SSIM4','CSSIM','CSSIM4'};
%%
for z = 1:1:numel(matrics_add)
    for i = 1:1:numel(imNameRef)
        % read a reference image
        img_ref = double(imread([pwd '\' folder_ref '\' ...
            num2str(imNameRef(i)) '.bmp']));
        % read a noisy image
        for j = 1:1:numel(sigma)
            
            img_n = double(imread([pwd '\' folder_n '\' ...
                'NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
            
            score(i,j) = IQA_score(img_ref, img_n, matrics_add{z});       

        end
    end
    
    clc
    disp(['Metrics: ' num2str(100*z/numel(matrics_add)) '%...'])
    
    save(strcat(matrics_add{z},'_n'), 'score');
end
%%
%% Filtered images processing 
for z = 1:1:numel(matrics_add)
    for i = 1:1:numel(imNameRef)
        % read a reference image
        img_ref = double(imread([pwd '\' folder_ref '\' ...
            num2str(imNameRef(i)) '.bmp']));
        % read a filtered image
        for j = 1:1:numel(sigma)
            
            img_f = double(imread([pwd '\' folder_f '\' ...
                'BM3D_NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
%             img_f = double(imread([pwd '\' folder_f '\' ...
%                 'DCTF_NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
            
            score(i,j) = IQA_score(img_ref, img_f, matrics_add{z});       

        end
    end
    
    clc
    disp(['Metrics: ' num2str(100*z/numel(matrics_add)) '%...'])
    
%     save(strcat(matrics_add{z},'_DCTF_f'), 'score');
    save(strcat(matrics_add{z},'_BM3D_f'), 'score');
end
%%
%%
%%
iqa_n = 'IQA Results\n';
iqa_dctf = 'IQA Results\f_dctf';
iqa_bm3d = 'IQA Results\f_bm3d';

metrics = [metrics_gray, metrics_color, matrics_add];
%%
addpath(genpath('IQA Results'));
load('Vote_DCTF.mat');
load('Vote_BM3D.mat');
%%
dctf_srocc = zeros(1, numel(metrics));
bm3d_srocc = zeros(1, numel(metrics));
%%
for i = 1:1:numel(metrics)
    score_n = load([metrics{i},'_n.mat']);
    score_dctf = load([metrics{i},'_DCTF_f.mat']);
    score_bm3d = load([metrics{i},'_BM3D_f.mat']);
    
    score_n = score_n.score';
    score_dctf = score_dctf.score';
    score_bm3d = score_bm3d.score';
    
    impr_dctf = score_dctf - score_n;
    impr_bm3d = score_bm3d - score_n;
    
    dctf_srocc(i) = corr(impr_dctf(:), Vote_DCTF(:), 'type', 'Spearman');
    bm3d_srocc(i) = corr(impr_bm3d(:), Vote_BM3D(:), 'type', 'Spearman');
end
%%
%%
clc
close all
clear
%%
imNameRef = [1 3 5 7 8 17 19 23 25:32];   % Reference images 
sigma = [3 5 10 15 20 25 30]; % noise std
folder_n = 'images\noisy';
%%
Trans2D = dctmtx(8);
ind = 1;
fb1_features = zeros(112,4);
fb2_features = zeros(112,4);
fb3_features = zeros(112,4);
fb4_features = zeros(112,4);
%%
for i = 1:1:numel(imNameRef)
    for j = 1:1:numel(sigma)
        
        img_n = double(imread([pwd '\' folder_n '\' ...
            'NSD_' num2str(sigma(j)) '_' num2str(imNameRef(i)) '.bmp']));
        
        Est = dct2Dblind_ratio( img_n, Trans2D, 'stripped_latitude' );
        fb1 = Est(:,:,1);
        fb2 = Est(:,:,2);
        fb3 = Est(:,:,3);
        fb4 = Est(:,:,4);
        
        fb1_features(ind, :) = [mean(fb1(:)), var(fb1(:)), kurtosis(fb1(:)), skewness(fb1(:))];
        fb2_features(ind, :) = [mean(fb2(:)), var(fb2(:)), kurtosis(fb2(:)), skewness(fb2(:))];
        fb3_features(ind, :) = [mean(fb3(:)), var(fb3(:)), kurtosis(fb3(:)), skewness(fb3(:))];
        fb4_features(ind, :) = [mean(fb4(:)), var(fb4(:)), kurtosis(fb4(:)), skewness(fb4(:))];
        
        ind = ind + 1;
        
    end
end


%%



