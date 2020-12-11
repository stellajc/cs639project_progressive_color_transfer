styleImage = imread('tar0.png');
styleImage = styleImage(1:224,1:224,:);
imshow(styleImage);
contentImage = imread('in0.png');
contentImage = contentImage(1:224,224:448,:);
imshow(contentImage);
% load pre-trainded vgg
net = vgg19;
% analyzeNetwork(net);
layername = ["relu1_1","relu2_1","relu3_1","relu4_1","relu5_1"];

[sw,sh,~] = size(styleImage);
[cw,ch,~] = size(contentImage);

% all possible combination for x y position in the patch
pairs = combvec(1:9,1:9);
S = contentImage;

for L = 5:-1:1
    
    % get the features for source image
    Fs = activations(net,S,layername(L));
    % get the features for ref image
    Fr = activations(net,styleImage,layername(L));
    % normalize channels
    Fsbar = normalize(Fs,3);
    Frbar = normalize(Fr,3);
    [xw, xh, xd] = size(Fsbar);
    [yw, yh, yd] = size(Frbar);
    Fswin = zeros(9,(xw-2)*(xh-2),xd);
    Frwin = zeros(9,(yw-2)*(yh-2),yd);
    % generate patches for FS
    for i = 1:xd
        Fswin(:,:,i) = im2col(Fsbar(:,:,i), [3 3], 'sliding');
    end
    % generate patches for FR
    for i = 1:yd
        Frwin(:,:,i) = im2col(Frbar(:,:,i), [3 3], 'sliding');
    end
    
    Fsnum = size(Fswin, 2);
    Frnum = size(Frwin, 2);
    cpFswin = Fswin;

    % iteration for BDS convergence only 
    for ite = 1:1000
        % S to R correspondence
        % go over all FS patches
        S2R = zeros(Fsnum,1);
        for i = 1:Fsnum
            % go over all FR patches
            idx = -1;
            minv = 99999;
            for j = 1:Frnum
                sumt = 0;
                for k = 1: 81
                    p = pairs(1,k);
                    q = pairs(2,k);
                    sumt = sumt + norm(reshape(cpFswin(p,i,:),1,xd)-reshape(Frwin(q,j,:),1,yd))^2;
                end
                if sumt < minv
                    idx = j;
                    minv = sumt;
                end
            end
            S2R(i,1) = idx;
        end

        % R to S correspondence
        % go over all FS patches
        R2S = zeros(Frnum,1);
        for i = 1:Frnum
            % go over all FR patches
            idx = -1;
            minv = 99999;
            for j = 1:Fsnum
                sumt = 0;
                for k = 1: 81
                    p = pairs(1,k);
                    q = pairs(2,k);
                    sumt = sumt + norm(reshape(Frwin(p,i,:),1,yd)-reshape(cpFswin(q,j,:),1,xd))^2;
                end
                if sumt < minv
                    idx = j;
                    minv = sumt;
                end
            end
            R2S(i,1) = idx;
        end
        
        % BDS voting
        newcpFswin = BDS(cpFswin,Frwin,S2R,R2S,xw, xh,yw, yh);
        newcpFswin = normalize(newcpFswin,3);
        
        % converge then break (may not converge, than manually select)
%         norm(sum(abs(newcpFswin-cpFswin),3))
        if norm(sum(abs(newcpFswin-cpFswin),3))<1000
            break
        end
        
        
        % get intermediate output G and FG
        FG = zeros(size(Fs,1)-2,size(Fs,2)-2,xd);
        G = zeros(size(S));
        for i = 1:Fsnum
            [srow,scol] = ind2sub([xw-2,xh-2],i);
            [rrow,rcol] = ind2sub([yw-2,yh-2],S2R(i));
            FG(srow,scol,:) = mean(Fr(rrow:rrow+2,rcol:rcol+2,:),[1,2]);
    %             refx = (rrow-1)*ceil(sw/(yw-2))+1:min(rrow*ceil(sw/(yw-2)),sw);
    %             refy = (rcol-1)*ceil(sh/(yh-2))+1:min(rcol*ceil(sh/(yh-2)),sh);
    %             refpatch = styleImage(refx,refy,:);
    %             meanref = mean(refpatch,[1,2]);
            meanref =  styleImage(min(round((2*rrow-1)*ceil(sw/(yw-2))/2),sw),...
                min(round((2*rcol-1)*ceil(sh/(yh-2))/2),sh),:);
            sourx = (srow-1)*ceil(cw/(xw-2))+1:min(srow*ceil(cw/(xw-2)),cw);
            soury = (scol-1)*ceil(ch/(xh-2))+1:min(scol*ceil(ch/(xh-2)),ch);
            for xtmp = sourx
                for ytmp = soury
                    G(xtmp,ytmp,:) = meanref(1,1,:);
                end
            end
        end
        
        % output G
        imwrite(uint8(G),append('G',num2str(L),'.png'));
        FG = activations(net,G,layername(L));
        err = sum((FG-Fs).^2,3);
        
        % local color transfer, update S
        S = local_color_transfer(G,L,contentImage,err);
        
        % output S
        imwrite(uint8(S),append('S',num2str(L),'.png'));

        update Fs patches
        cpFswin = newcpFswin;
    end

end
    