%% Coarse to Fine registration written by Kim(2019.11.21)
%% Add path
addpath('E:\Codes\envi');
addpath('E:\Codes\functions');
addpath('E:\Data\Sangju\UAV');
addpath('E:\Taeheon\UAV_KSRs');
addpath('E:\Taeheon\UAV_sanju\orthophoto');
addpath('E:\Taeheon\UAV_sanju\orthophoto\gcp_assessment\geotiff');
addpath('E:\Taeheon\UAV_sanju\orthophoto\gcp_assessment');
clc;clear;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data with header and image
cd('E:\Taeheon\UAV_sanju\orthophoto');
[f1,path] = uigetfile({'*.tif; *.jpg; *.bmp; *.gif'}, 'Please Select The Reference image');
if path == 0
    return;
end

cd('E:\Taeheon\UAV_sanju\orthophoto');
[f2,path] = uigetfile({'*.tif; *.jpg; *.bmp; *.gif'}, 'Please Select The Subject image');
if path == 0
    return;
end

tic
image1=imread(f1);
image2=imread(f2);

[row1,col1,band1] = size(image1);

[row2,col2,band2] = size(image2);
info1 = imfinfo(f1); % Reference image 헤더파일
info2 = imfinfo(f2); % Subject image 헤더파일
psize1 = info1(1).ModelPixelScaleTag(1:2);
psize2 = info2(1).ModelPixelScaleTag(1:2);
sratio = mean(psize1./psize2) % scale 계산
format long
GCP_coordinate1 = info1(1).ModelTiepointTag(4:5)    % upper-left
GCP_coordinate2 = info2(1).ModelTiepointTag(4:5)    % upper-left
%% transform the scale and translation of Sensed image
[Row,Col,Band]=size(image1);
img2 = imresize(image2, 1/sratio); % img2's scale -> img1's scale (higher res. should be reduced.)
psize2 = psize2*sratio; 

dxy = (GCP_coordinate2 - GCP_coordinate1)./psize2;
T=[1 0 0; 0 1 0; dxy(1) -dxy(2) sratio]/sratio;
tform = affine2d(T);
RA = imref2d(size(image1));

img2 = imwarp(img2,tform,'OutputView',RA);
mosaic_image_generation_3color(uint8(image1(:,:,1:3)),uint8(img2(:,:,[3 2 1])),min(Row,Col)/10); title('checkboard by original data(original)');
% figure,imshowpair(image1(:,:,1),image2(:,:,1))
%  figure,imshowpair(image1(:,:,1),img2(:,:,1))

% addpath('D:\Taeheon\UAV_sanju\orthophoto\ENVI');   % data path (K3A image)
% image1_hdr='20190801_reference';
% [image112, info1_1]=enviread(image1_hdr);
% cd 'D:\Taeheon\UAV_sanju\orthophoto\non_registraed';
% filename = '20191004_orthophoto_non';
% enviwrite(img2,info1_1,filename);

% img2_dem=img2(:,:,5);
% %% Finding overlapping blocks
% image1(image1<0)=0;
% img2_dem(-10<img2_dem)=1;
% img1 = image1;
% img_overlap =double(img1(:,:,5)).*double(img2_dem);

%% Finding overlapping blocks
image1(image1<0)=0;
img2(img2<0)=0;
img1 = image1;
img_overlap =double(img1(:,:,5)).*double(img2(:,:,5));
bw = img_overlap>0;
bw = imfill(bw,'holes');
ov  = regionprops(double(bw),'ConvexHull');
ov(1).ConvexHull = fix(ov(1).ConvexHull)
area=round(ov(1).ConvexHull);
area_sort=sort(area,2);
imshowpair(img1(:,:,1),img2(:,:,1)); hold on; plot(area(:,1),area(:,2),'k','LineWidth',3);

ov  = regionprops(bw,'BoundingBox');
ov(1).BoundingBox = round(ov(1).BoundingBox);
img1r = img1(ov(1).BoundingBox(2):ov(1).BoundingBox(2)+ov(1).BoundingBox(4),ov(1).BoundingBox(1):ov(1).BoundingBox(1)+ov(1).BoundingBox(3),:);
img2r = img2(ov(1).BoundingBox(2):ov(1).BoundingBox(2)+ov(1).BoundingBox(4),ov(1).BoundingBox(1):ov(1).BoundingBox(1)+ov(1).BoundingBox(3),:);
bw = bw(ov(1).BoundingBox(2):ov(1).BoundingBox(2)+ov(1).BoundingBox(4),ov(1).BoundingBox(1):ov(1).BoundingBox(1)+ov(1).BoundingBox(3));

%% Working on the downsampling images(계산 효율을 위해)
plevel = 1;
img1_down = imresize(img1r,1/2^plevel);
img2_down = imresize(img2r,1/2^plevel);
%% Outliers removal by iterative affine Transformation estimation
% Std-Mean Linear Stretch 
% Std-Mean Linear Stretch는 최대한 많은 특징점을 추출하기 위해 활용
Ref_use=double(img1_down(:,:,1));
Sub_use=double(img2_down(:,:,1));
Ref_mean=mean(mean(Ref_use));
Ref_std=std(Ref_use(:));
Ref_max=Ref_mean+Ref_std;
Ref_min=Ref_mean-Ref_std;
Ref_k=((Ref_use-Ref_min)./(Ref_max-Ref_min))*255;
Ref_use=uint8(Ref_k);
Sub_mean=mean(mean(Sub_use)); 
Sub_std=std(Sub_use(:));
Sub_max=Sub_mean+Sub_std;
Sub_min=Sub_mean-Sub_std;
Sub_k=((Sub_use-Sub_min)./(Sub_max-Sub_min))*255;
Sub_use=uint8(Sub_k);
%% Feature extraction using SURF method and searching space %%
% find the SURF features
points1 = detectSURFFeatures(Ref_use);
[f1, vpts1] = extractFeatures(Ref_use, points1,'Upright',true);
points2 = detectSURFFeatures(Sub_use);
[f2, vpts2] = extractFeatures(Sub_use, points2,'Upright',true);
fpoint_ref=vpts1.Location;ref_loc_x=fpoint_ref(:,1);ref_loc_y=fpoint_ref(:,2);
fpoint_sub=vpts2.Location;sub_loc_x=fpoint_sub(:,1);sub_loc_y=fpoint_sub(:,2);

% find the matchedfeatures using searching space
bl_size=200;
bl_move=fix(bl_size/2)
idx=1;k=0;

for i=1:length(fpoint_ref)
    k=k+1;
    CPs_r=[ref_loc_x(i),ref_loc_y(i)];
    ver_size=[CPs_r(:,1)-bl_move,CPs_r(:,1)+bl_move];
    hor_size=[CPs_r(:,2)-bl_move,CPs_r(:,2)+bl_move];
    CPs_sub_loc=find(min(ver_size)<=sub_loc_x & sub_loc_x<=max(ver_size) & min(hor_size)<=sub_loc_y...
        & sub_loc_y<=max(hor_size));
    [indexPairs, matchmetric] = matchFeatures(f1(i,:), f2(CPs_sub_loc,:), 'Method','Approximate','Metric', 'SSD', 'MaxRatio',0.6);  
        if(length(indexPairs)~=0)
            sub_data=CPs_sub_loc(indexPairs(:,2));
            ref_CPs=vpts1(i);
            sub_CPs=vpts2(sub_data);
            final_rCPs(idx,:)=ref_CPs.Location;
            final_sCPs(idx,:)=sub_CPs.Location;
            idx=idx+1;
        end

    clearvars indexPairs
    clearvars matchmetric 
end


img1r_vi=img1r; img1r_vi(img1r_vi==255)=nan;
img2r_vi=img2r; img2r_vi(img2r_vi==255)=nan;
img1_down(img1_down==255)=nan;
img2_down(img2_down==255)=nan;
% figure; showMatchedFeatures(uint8(img1_down(:,:,1:3)),uint8(img2_down(:,:,1:3)),final_rCPs,final_sCPs,'montage');
% legend('Reference tie-points','Sensed tie-points');

% Outliers elimination using RMSE and affine 
TH_R =10;
matchedRefPoints1=double(final_rCPs);
matchedSubPoints2=double(final_sCPs);
[matchedSubPoints_final, matchedRefPoints_final] = outlier_removal_iter_affine(matchedSubPoints2,matchedRefPoints1, TH_R);
[matchedSubPoints_final, matchedRefPoints_final] = PreProcessCp2tform(double(matchedSubPoints_final), double(matchedRefPoints_final));
% 
% [tform,inlierIdx] = estimateGeometricTransform2D(matchedSubPoints2,matchedRefPoints1,'affine','MaxNumTrials',2000,'MaxDistance',5);
% matchedRefPoints_final=matchedRefPoints1(inlierIdx,:);
% matchedSubPoints_final=matchedSubPoints2(inlierIdx,:);

matchedRefPoints_final=matchedRefPoints_final*2^plevel;
matchedSubPoints_final=matchedSubPoints_final*2^plevel;

figure; showMatchedFeatures(uint8(img1r_vi(:,:,1:3)),uint8(img2r_vi(:,:,1:3)),matchedRefPoints_final,matchedSubPoints_final,'montage');
% legend('Reference tie-points','Sensed tie-points'); 
mytform = cp2tform(matchedSubPoints_final, matchedRefPoints_final, 'affine');
registered_coarse = imtransform(img2r, mytform, 'nearest', 'FillValues', 255, 'XData', [1 size(img1r,2)],'YData', [1 size(img1r,1)]);
%   mosaic_image_generation_3color(uint8(img1r(:,:,1:3)),uint8(img2r(:,:,[3 2 1])),1500); title('original image');
%   mosaic_image_generation_3color(uint8(img1r(:,:,1:3)),uint8(registered_coarse(:,:,[3 2 1])),1500); title('co-registered image');

%% Coarse Warping for full scene
matchedRefPoints_final=double(matchedRefPoints_final);
matchedSubPoints_final=double(matchedSubPoints_final);

matchedRefPoints_final_full = 1+(matchedRefPoints_final+ov(1).BoundingBox(1:2)-2);
matchedSubPoints_final_full = 1+(matchedSubPoints_final+ov(1).BoundingBox(1:2)-2);

% figure; showMatchedFeatures(uint8(image1(:,:,1:3)),uint8(img2(:,:,1:3)),matchedSubPoints_final_full,matchedRefPoints_final_full,'montage');
% legend('Reference tie-points','Subject tie-points'); 

mytform_full = cp2tform(matchedSubPoints_final_full, matchedRefPoints_final_full, 'affine');
coarse_tform=affine2d(mytform_full.tdata.T);
% Translation을 sensed 영상 스케일로 변환(이부분은 안들어가는게 맞는거 같은데 나중에 다시한번 확인!!)
% mytform_full.tdata.Tinv(3,1:2) = mytform_full.tdata.Tinv(3,1:2)*sratio; %
registered_coarse_full = imtransform(img2, mytform_full, 'nearest', 'FillValues', 0, 'XData', [1 size(image1,2)],'YData', [1 size(image1,1)]);
mosaic_image_generation_3color(uint8(image1(:,:,1:3)),uint8(img2(:,:,[3 2 1])),min(Row,Col)/10); title('checkboard by original data(original)');
mosaic_image_generation_3color(uint8(image1(:,:,1:3)),uint8(registered_coarse_full(:,:,[3 2 1])),min(Row,Col)/10); title('checkboard by coarse registered data(registered)');

 CCover = find(bw~=0);
 [CC_original] = corrcoef(double(img1r(CCover)),double(img2r(CCover)));
 [CC_ERN] = corrcoef(double(registered_coarse(CCover)),double(img1r(CCover)));
 disp(['CC_original: ', num2str(CC_original(1,2))])
 disp(['CC_Coarse: ',num2str(CC_ERN(1,2))])
%% Accuracy of coarse registration result using CC(removing regions where DN=0)
registered_coarse_full(registered_coarse_full<0)=0;
img_overlap =double(img1(:,:,5)).*double(registered_coarse_full(:,:,5));
bw2 = img_overlap>0;
bw2 = imfill(bw2,'holes');
ov2  = regionprops(bw2,'BoundingBox');
ov2(1).BoundingBox = round(ov2(1).BoundingBox);
img1r_fine = img1(ov2(1).BoundingBox(2):ov2(1).BoundingBox(2)+ov2(1).BoundingBox(4),ov2(1).BoundingBox(1):ov2(1).BoundingBox(1)+ov2(1).BoundingBox(3),:);
coarse = registered_coarse_full(ov2(1).BoundingBox(2):ov2(1).BoundingBox(2)+ov2(1).BoundingBox(4),ov2(1).BoundingBox(1):ov2(1).BoundingBox(1)+ov2(1).BoundingBox(3),:);
bw2 = bw2(ov2(1).BoundingBox(2):ov2(1).BoundingBox(2)+ov2(1).BoundingBox(4),ov2(1).BoundingBox(1):ov2(1).BoundingBox(1)+ov2(1).BoundingBox(3));
% figure,imshow(bw2,[])

%% Fine registration using local MI(Code by Kim)
[Row, Col, Band] = size(img1r);
tsize = 450;
ssize = 200;
num = fix(Row/ssize)*fix(Col/ssize);
CP_s = zeros(num,2);
CP_r = zeros(num,2);
ind = zeros(num,1);
k=0;
 [optimizer,metric] = imregconfig('multimodal');
 optimizer.InitialRadius = optimizer.InitialRadius/3.5;
%  optimizer.MaximumIterations = 200;

    for i = 1 : ssize : Row-ssize
        for j = 1 : ssize : Col-ssize
            k=k+1;
            if(min(min(bw(i:i+ssize-1,j:j+ssize-1))) && i+tsize<Row && j+tsize<Col)% && i+ssize>tsize &&j+ssize>tsize)
                img1t = img1r(i:i+tsize-1,j:j+tsize-1,1);
                img2t = registered_coarse(i:i+tsize-1,j:j+tsize-1,1);
                %                         tform = imregcorr(img1t,img2t,'translation');
                tform = imregtform(img1t,img2t,'translation',optimizer,metric);
                xy = [tform.T(3,1), tform.T(3,2)];
                if(xy~=0)
                    ind(k) = 1;
                    CP_r(k,:) = [(2*j+tsize-1)/2,(2*i+tsize-1)/2];
                    CP_s(k,:) = [(2*j+tsize-1)/2+xy(1),(2*i+tsize-1)/2+xy(2)];
                end
            end
        end
    end


CP_rf = CP_r(ind==1,:); 
CP_sf = CP_s(ind==1,:);
%% Outliers removal by iterative affine Transformation estimation
  TH_R = 10;
  [CP_s_final, CP_r_final] = outlier_removal_iter_affine(CP_sf,CP_rf, TH_R);
  [CP_s_final, CP_r_final] = PreProcessCp2tform(CP_s_final, CP_r_final);
% 
% [tform,inlierIdx] = estimateGeometricTransform2D(CP_sf,CP_rf,'affine','MaxNumTrials',2000,'MaxDistance',5);
% CP_r_final=CP_rf(inlierIdx,:);
% CP_s_final=CP_sf(inlierIdx,:);
%% Visual assessment
img1r_vi2=img1r;img1r_vi2(img1r_vi2==255)=nan;
coarse_vi=registered_coarse; coarse_vi(coarse_vi==255)=nan;
figure; showMatchedFeatures(uint8(img1r_vi2(:,:,1:3)),uint8(coarse_vi(:,:,1:3)),double(CP_rf),double(CP_sf),'montage');
legend('Reference tie-points','Subject tie-points'); title('Matching-points')
% 
 figure; showMatchedFeatures(uint8(img1r_vi2(:,:,1:3)),uint8(coarse_vi(:,:,1:3)),CP_r_final,CP_s_final,'montage');
% legend('Reference tie-points','Subject tie-points'); title('Matching-points')

%% Sensed image Warping (affine transformation)
mytform = cp2tform(CP_s_final, CP_r_final, 'affine');
registered_lpc = imtransform(registered_coarse, mytform, 'bilinear', 'FillValues', 0, 'XData', [1 size(registered_coarse,2)],'YData', [1 size(registered_coarse,1)]);

%  mosaic_image_generation_3color(uint8(img1r(:,:,1:3)),uint8(registered_coarse(:,:,[3 2 1])),min(Row,Col)/10); title('original image');
%  mosaic_image_generation_3color(uint8(img1r(:,:,1:3)),uint8(registered_lpc(:,:,[3 2 1])),min(Row,Col)/10); title('co-registered image');
%% Fine Warping for full scene

CP_r_final_full = 1+(CP_r_final+ov(1).BoundingBox(1:2)-2);
CP_s_final_full = 1+(CP_s_final+ov(1).BoundingBox(1:2)-2);

% figure; showMatchedFeatures(uint8(image1(:,:,1:3)),uint8(registered_coarse_full(:,:,1:3)),CP_r_final_full,CP_s_final_full,'montage');
% legend('Reference tie-points','Subject tie-points'); title('Matching-points')

 fine_mytform_full = cp2tform(CP_s_final_full,CP_r_final_full, 'affine');
 fine_tform=affine2d(fine_mytform_full.tdata.T);
 registered_full_fine = imtransform(registered_coarse_full, fine_mytform_full, 'nearest', 'FillValues', 0, 'XData', [1 size(registered_coarse_full,2)],'YData', [1 size(registered_coarse_full,1)]);

mosaic_image_generation_3color(uint8(image1(:,:,1:3)),uint8(img2(:,:,[3 2 1])),min(Row,Col)/10); title('checkboard by coarse registered data');
mosaic_image_generation_3color(uint8(image1(:,:,1:3)),uint8(registered_full_fine(:,:,[3 2 1])),min(Row,Col)/10); title('checkboard by fine registered data');


% Only SURF
clearvars points1; clearvars points2; clearvars f1; clearvars f2; clearvars vpts1; clearvars vpts2;
points1 = detectSURFFeatures(Ref_use);
[f1, vpts1] = extractFeatures(Ref_use, points1,'Upright',true);
points2 = detectSURFFeatures(Sub_use);
[f2, vpts2] = extractFeatures(Sub_use, points2,'Upright',true);
[indexPairs, matchmetric] = matchFeatures(f1, f2, 'Method','Approximate','Metric', 'SSD', 'MaxRatio',0.6); 

matchedRefPoints1 = vpts1(indexPairs(:, 1));
matchedSubPoints2 = vpts2(indexPairs(:, 2));

figure; showMatchedFeatures(uint8(img1_down(:,:,1:3)),uint8(img2_down(:,:,1:3)),matchedRefPoints1,matchedSubPoints2,'montage');
legend('Reference tie-points','Sensed tie-points');

[matchedSubPoints_final, matchedRefPoints_final] = outlier_removal_iter_affine(matchedSubPoints2.Location,matchedRefPoints1.Location, TH_R);
[matchedSubPoints_final, matchedRefPoints_final] = PreProcessCp2tform(double(matchedSubPoints_final), double(matchedRefPoints_final));
matchedRefPoints_final=matchedRefPoints_final*2^plevel;
matchedSubPoints_final=matchedSubPoints_final*2^plevel;
% img1r_vi=img1r; img1r_vi(img1r_vi==255)=0;
% img2r_vi=img2r; img2r_vi(img2r_vi==255)=0;
figure; showMatchedFeatures(uint8(img1r_vi(:,:,1:3)),uint8(img2r_vi(:,:,1:3)),matchedRefPoints_final,matchedSubPoints_final,'montage');
legend('Reference tie-points','Sensed tie-points'); 

mytform = cp2tform(matchedSubPoints_final, matchedRefPoints_final, 'affine');
registered_surf = imtransform(img2r, mytform, 'nearest', 'FillValues', 255, 'XData', [1 size(img1r,2)],'YData', [1 size(img1r,1)]);
 mosaic_image_generation_3color(uint8(img1r(:,:,1:3)),uint8(img2r(:,:,[3 2 1])),1500); title('original image');
 mosaic_image_generation_3color(uint8(img1r(:,:,1:3)),uint8(registered_surf(:,:,[3 2 1])),1500); title('co-registered image');
CCover = find(bw~=0);
[CC_ERN] = corrcoef(double(registered_surf(CCover)),double(img1r(CCover)));

disp(['CC_SURF: ',num2str(CC_ERN(1,2))])

 %% Accuracy of fine registration result using CC (removing regions where DN=0)
CCover = find(bw~=0);
[CC_original] = corrcoef(double(img1r(CCover)),double(img2r(CCover)));
[CC_ERN] = corrcoef(double(registered_coarse(CCover)),double(img1r(CCover)));
disp(['CC_original: ', num2str(CC_original(1,2))])
disp(['CC_Coarse: ',num2str(CC_ERN(1,2))])

CCover = find(bw~=0);
[CC_ERN] = corrcoef(double(registered_lpc(CCover)),double(img1r(CCover)));
disp(['CC_Fine: ',num2str(CC_ERN(1,2))])



cd('E:\Taeheon\UAV_sanju\orthophoto\Accuracy assessment\checkpoint');
[f4,path] = uigetfile({'*.csv;*.xlsx;'}, 'Please Select file');
if path == 0
    return;
end
[~,sheet_name]=xlsfinfo(f4);
parfor k=1:numel(sheet_name)
  data{k}=xlsread(f4,sheet_name{k});
end
reference_check=cell2mat(data(1)); reference_check(isnan(reference_check))=0; reference_check=reference_check(3:end,2:end)
sensed_check=cell2mat(data(2)); sensed_check(isnan(sensed_check))=0; sensed_check=sensed_check(3:end,2:end)
dif_check=(reference_check-sensed_check).^2
sum_dif_check=dif_check(:,1)+dif_check(:,2)
RMSE_raw=sqrt(sum(sum_dif_check)/length(dif_check))
coarse_CP_temp = transformPointsForward(coarse_tform,sensed_check);
fine_CP_temp = transformPointsForward(fine_tform,coarse_CP_temp);

dif_surftemp=(reference_check-coarse_CP_temp).^2;
sum_difsurf_temp=dif_surftemp(:,1)+dif_surftemp(:,2);
RMSE_SURF=sqrt(sum(sum_difsurf_temp)/length(dif_check));

dif_temp=(reference_check-fine_CP_temp ).^2;
sum_dif_temp=dif_temp(:,1)+dif_temp(:,2);
RMSE_reg=sqrt(sum(sum_dif_temp)/length(dif_check));
fprintf('==========RMSE========');
RMSE_raw
RMSE_SURF
RMSE_reg
toc
% 
% %% imwrite_non
% registered_full_fine(registered_full_fine<=0)=255;
% img2(img2<=0)=255;
% figure,imshow(uint8(registered_full_fine(:,:,1:3)))
%  addpath('D:\Taeheon\UAV_sanju\orthophoto\ENVI');   % data path (K3A image)
% image1_hdr='20190801_reference';
% [image112, info1_1]=enviread(image1_hdr);
% cd 'D:\Taeheon\UAV_sanju\orthophoto\registered';
% %  info11 = info1(1);
% %  info11.GeoKeyDirectoryTag = [];
% %  imwrite2tif(registered_full_fine,info11,'20190621_orthophoto_reg.tif','double');
% filename = '20191004_orthophoto_reg';
% enviwrite(registered_full_fine,info1_1,filename);
% %  
% % cd 'D:\Taeheon\UAV_sanju\orthophoto\non_registraed';
% % %  info11 = info1(1);
% % %  info11.GeoKeyDirectoryTag = [];
% % %  imwrite2tif(img2,info11,'20190621_orthophoto_non.tif','double');
% % filename = '20190621_orthophoto_non';
% % enviwrite(img2,info1_1,filename);
% % 
% cd 'D:\Taeheon\UAV_sanju\orthophoto\SURF';
% %  info11 = info1(1);
% %  info11.GeoKeyDirectoryTag = [];
% %  imwrite2tif(img2,info11,'20190621_orthophoto_non.tif','double');
% filename = '20190820_orthophoto_surf';
% enviwrite(registered_coarse_full,info1_1,filename);
% % 
% cd 'D:\Taeheon\UAV_sanju\orthophoto\MI';
% %  info11 = info1(1);
% %  info11.GeoKeyDirectoryTag = [];
% %  imwrite2tif(img2,info11,'20190621_orthophoto_non.tif','double');
% filename = '20190820_orthophoto_MI';
% enviwrite(registered_MI,info1_1,filename);
% %% imwrite_1
% 
%  cd 'D:\Taeheon\UAV_sanju\orthophoto\non_registraed';
%  info11 = info1(1);
%  info11.GeoKeyDirectoryTag = [];
%  imwrite2tif(registered_coarse,info11,'20190920_orthophoto_surf.tif','double');
% 
% %{
%  %% imwrite
% 
%  cd 'D:\Taeheon\UAV_sanju\orthophoto\registered';
%  info11 = info1(1);
%  info11.GeoKeyDirectoryTag = [];
%  imwrite2tif(registered_lpc,info11,'20190906_orthophoto.tif','double');
%  %}
%  
%  
%  addpath('D:\Taeheon\UAV_sanju\orthophoto\ENVI');   % data path (K3A image)
% image1_hdr='20190801_reference';
% [image112, info1_1]=enviread(image1_hdr);
