function newcpFswin = BDS(cpFswin,Frwin,S2R,R2S,xw, xh,yw, yh)

    Fsnum = size(S2R,1);
    Frnum = size(R2S,1);
    newcpFswin = zeros(size(cpFswin));
    
    for i = 1:Fsnum
        [srow,scol] = ind2sub([xw-2,xh-2],i);

        m = 0;
        n = 0;
        % go over surounding src feature positions
        for x = max(srow-2,1):min(xw,srow+2)
            for y = max(scol-2,1):min(xh,scol+2)
                if x<= xw-2 && y <= xh-2
                    % surounding src feature idx
                    idx = sub2ind([xw-2,xh-2],x,y);
                    ridx = S2R(idx,1);
                    % find corresponding position in ref feature
                    [tmpx,tmpy]=ind2sub([yw-2,yh-2],ridx);
                    % x,y corresponding position in ref patch
                    s2rx = tmpx-x+srow;
                    s2ry = tmpy-y+scol;
                    if tmpx-x+srow>=1 &&tmpx-x+srow<=yw-2&& ...
                            tmpy-y+scol>=1&&tmpy-y+scol<=yh-2
                        ridx = sub2ind([yw-2,yh-2],s2rx,s2ry);
                        newcpFswin(:,i,:) = newcpFswin(:,i,:)+ Frwin(:,ridx,:)/Fsnum;
                        m = m +1;
                    end
                    % points in ref pixel mapping to surrouding
                    % positions in the src feature
                    r2sidx = find(R2S==idx);
                    for k = 1:size(r2sidx)
                        [r2sx,r2sy] = ind2sub([yw-2,yh-2],r2sidx(k));
                        tx = r2sx+srow-x;
                        ty = r2sy+scol-y;
                        if r2sx+srow-x>=1 && r2sx+srow-x<=yw-2&&...
                            r2sy+scol-y>=1 && r2sy+scol-y<=yh-2
                            ridx2 = sub2ind([yw-2,yh-2],tx,ty);
                            newcpFswin(:,i,:)=newcpFswin(:,i,:)+Frwin(:,ridx2,:)/Frnum;
                            n = n+1;
                        end
                    end
                end
            end
        end
        newcpFswin(:,i,:) = newcpFswin(:,i,:)/(n/Frnum+m/Fsnum);
    end
end
%         for i = 1:Frnum
%             [tmpx,tmpy]=ind2sub([yw-2,yh-2],i);
%             [srow,scol] = ind2sub([xw-2,xh-2],R2S(i,1));
%             for x = max(srow-2,1):min(xw,srow+2)
%                 for y = max(scol-2,1):min(xh,scol+2)
%                     if x<= xw-2 &&y <= xh-2
%                         idx = sub2ind([xw-2,xh-2],x+tmpx-srow,y+tmpy-scol);
%                         newcpFswin(:,idx,:) = newcpFswin(:,idx,:)+ Frwin(:,i,:)/Frnum;
% %                         n = n + 1;
%                     end
%                 end
%             end
%         end