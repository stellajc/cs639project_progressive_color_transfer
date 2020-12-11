function [modified_src] =  local_color_transfer(g_img, level, src_img, matching_err_l)

% downscale S to S(l) to match G; convert rgb to lab
%g_img = im2double(imread('brook1.jpg'));
%src_img = im2double(imread('src3.jpg'));

[H, W, d] = size(g_img);
[h,w, d1] = size(src_img);
sl = rgb2lab(imresize(src_img,[H, W]));
gl = rgb2lab(g_img);

% initialization of Tl
Tl = zeros(H, W, 3, 2); % 3 = # channels; 2 = # coef at a pixel
%patch_states = zeros(H,W,3,4); % 4=#stats to compute coef a, b
e = 0.0002; %color range [0,1]
e1 = 1.3;
for c = 1:3
    for i = 1:H
        for j = 1:W
            [m_s, sd_s, m_g, sd_g] = compute_patch_stats(i,j,sl(:,:,c),gl(:,:,c),H,W);
            al = e1*sd_g/(sd_s + e);
            bl =  m_g - al*m_s;
            Tl(i,j,c,1) = al;
            Tl(i,j,c,2) = bl;
        end
    end
end
m_src = zeros(H,W,3);
for c = 1:3
    for i = 1:H
        for j = 1:W
            m_src(i,j,c) = Tl(i,j,c,1)*sl(i,j,c)+Tl(i,j,c,2);
        end
    end
end

% compute_ed(Tl,sl,gl,level,matching_err_l);
% compute_local(Tl, sl);
modified_src = imresize(lab2rgb(m_src),[h, w]);% this is to convert back
% to it.
%imshow(m_src);
%figure;
%imwrite(m_src,"output.jpg","jpg");

end
function w = wls(e,a,p,q,a_p,a_q,b_p,b_q)
    w = 1/(pow(norm(p-q),a) + e)*(pow(norm((a_p - a_q)),2)+ pow(norm((b_p - b_q)),2));
end
% guidance term
function ed = compute_ed(Tl,sl,gl,level,matching_err_l)
    wl = pow(4, level-1);
    el = matching_err_l;
    E_d = 0;
    for c = 1:3
        for i = 1:H
            for j = 1:W
                x = (Tl(i,j,c,1)*sl(i,j,c)+Tl(i,j,c,2) - gl(i,j,c)).^2;
                E_d = E_d + wl*(1-el(i,j,c))*(pow(x,2));
            end
        end
    end
    ed = E_d;
end
function l = compute_local(Tl,sl)
% local constraint term 
lambda_l = 0.125;
e = 0.0001;
a = 1.2;
E_l = 0;
    for i = 1:H
        for j = 1:W
            el = 0;
            p = sl(i,j);
            if i > 1
                q = sl(i-1, j);
                el = el + wls(e,a,p,q,Tl(i,j,:,1),Tl(i-1,j,:,1),Tl(i,j,:,2),Tl(i-1,j,:,2));
            end
            if j > 1
                q =sl(i, j-1); 
                el = el + wls(e,a,p,q,Tl(i,j,:,1),Tl(i,j-1,:,1),Tl(i,j,:,2),Tl(i,j-1,:,2));
            end
            if j+1 < W
                q = sl(i,j+1);
                el = el + wls(e,a,p,q,Tl(i,j,:,1),Tl(i,j+1,:,1),Tl(i,j,:,2),Tl(i,j+1,:,2));
            end
            if i+1 < H
                q = sl(i+1,j);
                el = el + wls(e,a,p,q,Tl(i,j,:,1),Tl(i+1,j,:,1),Tl(i,j,:,2),Tl(i+1,j,:,2));
            end
            E_l = E_l + el;
        end
    end
    l = lambda_l * E_l;
end 

function [m_s, sd_s, m_g, sd_g] = compute_patch_stats(i,j,s,g,h,w)
% compute the mean and sd for the local patch centered at (i,j)
    patch_s = [];
    patch_g = [];
    %ul
    if i == 1 && j == 1
        patch_s = [s(i,j),s(i+1,j),s(i+1,i+1),s(i,j+1)];
        patch_g = [g(i+1,j),g(i+1,i+1),g(i,j+1)];
    %ur
    elseif i == 1 && j == w
        patch_s = [s(i,j-1),s(i,j),s(i+1,j),s(i+1,j-1)];
        patch_g = [g(i,j-1),g(i,j),g(i+1,j),g(i+1,j-1)];
    %ll
    elseif i == h && j == 1
        patch_s = [s(i,j),s(i-1,j),s(i-1,j+1),s(i,j+1)];
        patch_g = [g(i,j),g(i-1,j),g(i-1,j+1),g(i,j+1)];
    %lr
    elseif i == h && j == w
        patch_s = [s(i,j),s(i,j-1),s(i-1,j),s(i-1,j-1)];
        patch_g = [g(i,j),g(i,j-1),g(i-1,j),g(i-1,j-1)];
    %upper
    elseif i == 1 && j < w
        patch_s = [s(i,j),s(i,j-1),s(i,j+1),s(i+1,j-1),s(i+1,j+1),s(i+1,j)];
        patch_g = [g(i,j-1),g(i,j),g(i+1,j),g(i+1,j-1),g(i+1,j+1),g(i+1,j)];
    %lower
    elseif i == h && j < w
        patch_s = [s(i,j),s(i,j-1),s(i-1,j),s(i-1,j-1),s(i,j+1), s(i-1,j+1)];
        patch_g = [g(i,j),g(i,j-1),g(i-1,j),g(i-1,j-1),g(i,j+1), g(i-1,j+1)];
    %left
    elseif i < h && j == 1
        patch_s = [s(i,j),s(i-1,j),s(i-1,j+1),s(i,j+1),s(i+1,j),s(i+1,j+1)];
        patch_g = [g(i,j),g(i-1,j),g(i-1,j+1),g(i,j+1),g(i+1,j),g(i+1,j+1)];
    %right
    elseif i < h && j == w
        patch_s = [s(i,j),s(i,j-1),s(i-1,j),s(i-1,j-1),s(i+1,j),s(i+1,j-1)];
        patch_g = [g(i,j),g(i,j-1),g(i-1,j),g(i-1,j-1),g(i+1,j),g(i+1,j-1)];
    %normal 9 point patch
    else
        patch_s = [s(i,j),s(i,j-1),s(i-1,j),s(i-1,j-1),s(i+1,j),s(i+1,j-1),s(i,j+1),s(i-1,j+1),s(i+1,j+1)];
        patch_g = [g(i,j),g(i,j-1),g(i-1,j),g(i-1,j-1),g(i+1,j),g(i+1,j-1),g(i,j+1),g(i-1,j+1),g(i+1,j+1)];
    end
        m_s = mean(patch_s,"all");
        m_g = mean(patch_g,"all");
        sd_s = std(patch_s);
        sd_g = std(patch_g);
     
end
