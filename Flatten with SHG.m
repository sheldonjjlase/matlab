clear
close all
tic
currentFolder=pwd;

[FileName,PathName] = uigetfile('*.hdf5;*.h5','Select the file');
cd(currentFolder)
%
trial=0;
optimize=0;


info=h5info(strcat(PathName,FileName));
dataset=info.Datasets.Name;

rawdata=h5read(strcat(PathName,FileName),['/' dataset]);

%% select shg
% timepoint=size(rawdata,5);
% allheight=cell(timepoint,1);
% rawdata=rawdata(1,:,:,:,:);
% rawdata=permute(squeeze(rawdata),[2 1 3 4]);

% timepoint=size(rawdata,5);
timepoint=1;
allheight=cell(timepoint,1);
rawdata=squeeze(rawdata(1,:,:,:));
rawdata=permute(rawdata,[2 1 3]);
% rawdata=permute(squeeze(rawdata),[2 1 3 4]);

%%
close all
trial=0;
for t=1:timepoint
    while(1)
        
        prompt = {'Enter XY resolution (um/pixel):','Enter Z resolution (um/pixel):',...
            'Enter approx. frame # where BM starts','XY blurring (in um, typically 10-25)',...
            'Z blurring (in um, typically 1-4)','Downsample resolution (um/pixel)',...
            '# of frames in each zstack','# of timepoints per z stack'};
        dlg_title = 'Input';
        num_lines = 1;
        
        if trial ==0 && t==1
    
                defaultans = {'1','3',...
                    '15','40',...
                    '3','2',...
                    '45',num2str(timepoint)};                  

            
        else
            defaultans = {num2str(res_xy), num2str(res_z),...
                num2str(BM_start), num2str(blur),...
                num2str(blurz), num2str(dsrate),...
                num2str(stacksize),num2str(timepoint)};
        end
        
        
        if optimize == 0
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            res_xy=str2num(answer{1}); %XY RESOLUTION
            res_z=str2num(answer{2}); %Z RESOLUTION
            BM_start=str2num(answer{3}); %ROUGHLY WHERE THE BASEMENT MEMBRANE STARTS
            blur=str2num(answer{4}); %AMOUNT TO BLUR
            blurz=str2num(answer{5}); %AMOUNT TO BLUR
            dsrate=str2num(answer{6}); %Down sampled resolution
            stacksize=str2num(answer{7}); %number of slices in the zstack
            timepoint=str2num(answer{8}); %SHG channel number
            
        else if optimize == 1
            end
        end
        
                blue=squeeze(rawdata(:,:,:,t));              

        stacksize=size(blue,3);
        movie=zeros(size(blue,1),size(blue,2),3,stacksize);
        
            movie(:,:,1,:)=blue;
            movie(:,:,2,:)=blue;
            movie(:,:,3,:)=blue;

        
        originalsize_x=size(movie,1);
        originalsize_y=size(movie,2);
        dsrate_adjust=(res_xy/dsrate);
        if originalsize_x*dsrate_adjust > floor(originalsize_x*dsrate_adjust)
        dssize_x=floor(originalsize_x*dsrate_adjust)+1;
        else
         dssize_x=floor(originalsize_x*dsrate_adjust);
        end
        if originalsize_y*dsrate_adjust > floor(originalsize_y*dsrate_adjust)
        dssize_y=floor(originalsize_y*dsrate_adjust)+1;
        else
         dssize_y=floor(originalsize_y*dsrate_adjust);  
        end
        
        zslices=size(movie,4);
        filtSHG=zeros(dssize_x,dssize_y,zslices);
        
        for z=1:zslices
            filtSHG(:,:,z)=imgaussfilt(imresize(movie(:,:,3,z),dsrate_adjust),blur);
            
%             j=figure(1); clf;
%             set(gcf,'Position',[225 420 1015 415]);        
%             subplot(1,2,1);
%             titlelabel=sprintf('Timepoint %d / %d; plane %d / %d',t,timepoint,z,zslices);
%             imshow(filtSHG(:,:,z),[0 255],...
%                 'InitialMagnification','fit');
%             title(titlelabel);
%             subplot(1,2,2);
%             imshow(blue(:,:,z),[0 255],'InitialMagnification','fit');
%             title(titlelabel);
            
%             pause(0.01);
        end
        
        % min(min(filtSHG(:,:,z))) max(max(filtSHG(:,:,z)))
        
%         close(j);
        
        threedblur=blurz/res_z;
        filtSHG=imgaussfilt3(filtSHG,threedblur);
        height=zeros(dssize_x,dssize_y);
        
        
        parfor x=1:dssize_x
            for y=1:dssize_y
                locs=[];
                test=squeeze(filtSHG(x,y,BM_start:end));
                [peak,locs]=findpeaks(test,'MinPeakProminence',3)
%                 [peak,locs]=findpeaks(test);
                [mtest,BM_end]=max(test);
                %  [~,locs]=findpeaks(test);
                if isempty(locs) == 1
                    if (BM_end + BM_start) > size(filtSHG,3)-10
                        height(x,y)=(BM_end + BM_start);
                    else;
                        height(x,y)=BM_start;
                        continue
                    end
                else
                    if mtest > 1.5*peak(end) && BM_end ~=1
                        height(x,y)=(BM_end + BM_start);
                    else
%                         [~,mp]=max(peak);
                        height(x,y)=locs(end)+BM_start;
                    end
                end
            end
        end
        
        
        
        boxsize=round(dssize_x/5);
        height2=height;
        %remove outliers
        for x=1:dssize_x
            for y=1:dssize_y
                if height(x,y)==BM_start
                    %check limits
                    if x<(boxsize+1)
                        minx=1;
                    else
                        minx=x-boxsize;
                    end
                    
                    if x>(dssize_x-boxsize)
                        maxx=dssize_x;
                    else
                        maxx=x+boxsize;
                    end
                    
                    if y<(boxsize+1)
                        miny=1;
                    else
                        miny=y-boxsize;
                    end
                    
                    if y>(dssize_y-boxsize)
                        maxy=dssize_y;
                    else
                        maxy=y+boxsize;
                    end
                    
                    height2(x,y)=median(median(height(minx:maxx,miny:maxy)));
                    
                end
            end
        end
        
        
        h=figure(2);
        height3=medfilt2(height2,[blur blur],'symmetric');
        height4=imresize(height3,1/dsrate_adjust);
        height5=imgaussfilt(height4,blur/dsrate_adjust);
        height6=round(height5);
        imagesc(height6);
        c=colorbar;
        ylabel(c,'z-slice');
        axis off
        titlelabel2=sprintf('Timepoint %d / %d',t,timepoint);
        title(titlelabel2);
        
        allheight{t}=height6;
        
        if optimize ==0
            promptMessage = sprintf('Do you want to continue with this heat map ,\n or try with different parameters?');
            button = questdlg(promptMessage, 'Continue?', 'Continue', 'Restart', 'Continue');
        end
        
        if strcmpi(button, 'Continue')
            optimize=1;
            break
        end
        trial=trial+1;
    end

end


%    %% normalize
%     
%     maxheight=max(max(height6));
%     minheight=min(min(height6));
%     heightdiff=maxheight-minheight;
%     newzslices=zslices+heightdiff;
%     newstack=zeros(originalsize_x,originalsize_y,3,newzslices);
%     
%     %pad zeros
%     movie2=cat(4,zeros(originalsize_x,originalsize_y,3,heightdiff),movie);
%     movie3=cat(4,movie2,zeros(originalsize_x,originalsize_y,3,heightdiff));
%     
%     f=waitbar(1,'Normalizing z-stack...');
%     for x=1:originalsize_x
%         waitbar(x/originalsize_x);
%         for y=1:originalsize_y
%             newstack(x,y,:,:)=movie3(x,y,:,height6(x,y)-minheight+[1:newzslices]);
%         end
%     end
%     
%     close(f);
%     
% 
%     
%     for z=1:newzslices
%         zprofile(z)=sum(sum(sum(newstack(:,:,:,z))));
%     end
%     [locmax,loc]=max(zprofile);
%     
%     locmin=20; %how many um below the zero plane (i.e. basement membrane)
%     locmax=loc+find(zprofile(loc:end)<0.25*locmax,1);
%     allmovies{t}=newstack(:,:,:,loc-round(locmin/res_z):locmax);
%     [a,b,c,d]=size(allmovies{t});
%     
%     if t >1
%         while size(allmovies{t},4) > size(allmovies{1},4);
%             allmovies{t}(:,:,:,end)=[];
%         end
%         while size(allmovies{t},4) < size(allmovies{1},4);
%             allmovies{t}(:,:,:,end+1)=zeros(a,b,c);
%         end
%     end
%     
%     if t==1
%         allmoviestiff=zeros([size(allmovies{1}) timepoint],'uint8');
%         originalstiff=zeros([size(movie) timepoint],'uint8');
%         
%           allmoviestiff(:,:,:,:,t)=uint8(allmovies{t}*255/4095);
%           originalstiff(:,:,:,:,t)=uint8(movie*255/4985);
%     end
%     
%     allmoviestiff(:,:,:,:,t)=uint8(allmovies{t}*255/4095);
%     originalstiff(:,:,:,:,t)=uint8(movie*255/4095);
%     
% %     
% % end
%%
% 
% addpath(currentFolder)
% cd(PathName)
% 
% % originalfilename=strcat(FileName(1:end-4),'.tiff');
% % [FileName2,PathName2] = uiputfile(originalfilename,'Save Original Data as');
% % cd(PathName2)
% % bfsave(originalstiff, FileName2,'dimensionOrder','XYCZT');
% 
% defaultfilename=strcat(FileName(1:end-4),'-NORM-BM-',num2str(BM_start),'-XY-',...
%     num2str(blur),'-Z-',num2str(blurz),'-DS-',num2str(dsrate),'.tiff');
% [FileName3,PathName3] = uiputfile(defaultfilename,'Save Normalized Data as');
% cd(PathName3)
% bfsave(allmoviestiff, FileName3,'dimensionOrder','XYCZT');
% 
% cd(currentFolder)
% 
% 
% 
% %
% %  %% 3d registration
% %
% % finalsize=size(allmovies{1},4);
% %
% % fixed=[];
% % moving=[];
% %
% % fixed=squeeze(allmovies{1}(:,:,2,16));
% % moving=squeeze(allmovies{2}(:,:,2,16));
% % moving2=squeeze(allmovies{3}(:,:,2,16));
% %
% % % figure;
% % % imshowpair(fixed, moving,'Scaling','joint');
% % tformEstimate = imregcorr(moving,fixed);
% % Rfixed = imref2d(size(fixed));
% % movingReg = imwarp(moving,tformEstimate,'OutputView',Rfixed);
% % figure;
% % imshowpair(fixed, movingReg,'Scaling','joint')
% %
% % tformEstimate2 = imregcorr(moving2,movingReg);
% % Rfixed2 = imref2d(size(movingReg));
% % movingReg2 = imwarp(moving2,tformEstimate2,'OutputView',Rfixed2);
% % figure;
% % imshowpair(movingReg, movingReg2,'Scaling','joint')
% %
% %
% % regallmovies=cell(length(files),1);
% % regallmovies{1}=allmovies{1};
% % for i=2:length(files)
% %     regallmovies{i}=zeros(size(allmovies{i}));
% % end
% %
% %
% % for z=1:finalsize
% %         regallmovies{2}(:,:,2,z)=imwarp(allmovies{2}(:,:,2,z),tformEstimate,'OutputView',Rfixed);
% %         regallmovies{2}(:,:,3,z)=imwarp(allmovies{2}(:,:,3,z),tformEstimate,'OutputVIew',Rfixed);
% %
% %         regallmovies{3}(:,:,2,z)=imwarp(allmovies{3}(:,:,2,z),tformEstimate2,'OutputView',Rfixed2);
% %         regallmovies{3}(:,:,3,z)=imwarp(allmovies{3}(:,:,3,z),tformEstimate2,'OutputVIew',Rfixed2);
% % end
% %
% % for i=1:length(files)
% % cd(PathName)
% %     [FileName2,PathName2] = uiputfile(defaultfilename,'Save Normalized z-stack as');
% %     v=VideoWriter(strcat(PathName2,FileName2),'Uncompressed AVI');
% %     open(v)
% %     writeVideo(v,regallmovies{i}/4095);
% %     close(v)
% %     cd(currentFolder)
% % end


%%

figure(3); clf;

for i=3:3
    titlelabel=sprintf('Timepoint %d / %d',i,timepoint);
figure(1); clf;
min_height=0.75*min(min(allheight{i}));
max_height=0.75*max(max(allheight{i}));
h=surf(0.75*imresize(allheight{i},3),'FaceLighting','gouraud','FaceColor','interp');
set(h,'LineStyle','none'); 
view([-35,24]);
zlim([20 40]);
% title(titlelabel);
% figure(2);
% plot(i,max_height-min_height,'o','MarkerSize',5); hold on
% xlabel('Time');
% ylabel('Max height - min height');
grid on
xlabel('x');
ylabel('y');
zlabel('z');

end



%%
%NORMALIZE AND SAVE DATA 



load allheight.mat
%
% data=h5read('register2.hdf5','/data');
% data=permute(data,[1 3 2 4 5]);

data=h5read('stitch.hdf5','/data');
data=permute(data,[1 3 2 4]);
%%
interpz=90;
heightstack=zeros(size(data,2),size(data,3),size(data,5));
for i=1:length(allheight)
    heightstack(:,:,i)=imresize(allheight{i},[size(data,2) size(data,3)]);
end

% normalize
xsize=size(data,2);
ysize=size(data,3);

qq=waitbar(0,'Analyzing Data..');
% for t=1:timepoint
    waitbar(t/timepoint,qq);
    movie=squeeze(data(:,:,:,:,t));
    movie=permute(movie,[2 3 1 4]);
    
    for jj=1:size(movie,3)
    temp{jj}=imresize3(squeeze(movie(:,:,jj,:)),[size(movie,1) size(movie,2) interpz]);
    end
    movie2=cat(4,temp{1},temp{2},temp{3});
    movie=[];
    movie=permute(movie2,[1 2 4 3]);
       
    height6=round(interpz*imadjust(heightstack(:,:,t)/interpz,[],[]));
    
    maxheight=max(max(height6));
    minheight=min(min(height6));
    heightdiff=maxheight-minheight;
    newzslices=interpz+heightdiff;
    newstack=zeros(xsize,ysize,3,newzslices);
    
    %pad zeros
    movie2=cat(4,zeros(xsize,ysize,3,heightdiff),movie);
    movie3=cat(4,movie2,zeros(xsize,ysize,3,heightdiff));
    
    f=waitbar(1,'Normalizing z-stack...');
    for x=1:xsize
        waitbar(x/xsize,f);
        for y=1:ysize
            newstack(x,y,:,:)=movie3(x,y,:,height6(x,y)-minheight+[1:newzslices]);
        end
    end
    
    close(f);

        
    for jj=1:size(newstack,3)
%     temp{jj}=imresize3(squeeze(newstack(:,:,jj,:)),[size(newstack,1) size(newstack,2) stacksize]);
        temp{jj}=squeeze(newstack(:,:,jj,:));

    end
    
    newstack2=cat(4,temp{1},temp{2},temp{3});
    newstack=[];
    newstack=permute(newstack2,[1 2 4 3]);
    
    

    allmovies{t}=newstack;
    [a,b,c,d]=size(allmovies{t});
    
    if t >1
        while size(allmovies{t},4) > size(allmovies{1},4);
            allmovies{t}(:,:,:,end)=[];
        end
        while size(allmovies{t},4) < size(allmovies{1},4);
            allmovies{t}(:,:,:,end+1)=zeros(a,b,c);
        end
    end
    
    if t==1
        allmoviestiff=zeros([size(allmovies{1}) timepoint],'uint8');
        allmoviestiff(:,:,:,:,t)=uint8(allmovies{t});
    end
    
%      allmoviestiff(:,:,:,:,t)=imadjustn(uint8(allmovies{t}));
     allmoviestiff(:,:,:,:,t)=uint8(allmovies{t});
% end
close(qq);

%%
savename='norm_interp2.h5';
if exist(savename, 'file')==2
  delete(savename);
end
h5create(savename,'/data',[size(allmoviestiff)],'Datatype','uint8');
h5write(savename,'/data',uint8(allmoviestiff));
