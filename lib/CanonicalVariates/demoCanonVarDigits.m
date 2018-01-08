% ==================================================== 
%  Copyright (C) 2016 Kyle Johnston
%  
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ====================================================

function demoCanonVarDigits
%DEMOCANONVARDIGITS Demo of Canonical Variates on Digits data
figure
p = 800; % number of examples

% To visualize the images in digits:
% imshow(reshape(x(10,:),28,28),[]);
% x(1,:) is a vectorized image.
load digit3; D{1} = x(1:p,:)';  % Puts it the for DxN format for CanonVar, where D is dims and N the number of data points.
load digit5; D{2} = x(1:p,:)';
load digit7; D{3} = x(1:p,:)';

marker={'+','o','d'};

numDims = 3;

% W=CanonVar(D,numDims);
W=CanonVarExtended(D,numDims);

if(numDims == 2)
    V1=W(:,1); V2=W(:,2);
    nclasses=length(D);
    figure; col=get(gca,'colororder');
    hold on
    for c=1:nclasses
        PD{c}(:,1) = (D{c})'*V1;
        PD{c}(:,2) = (D{c})'*V2;
        h = plot(PD{c}(:,1),PD{c}(:,2),marker{c},'markersize',4,'color',col(c,:));
    end
    set(gca,'box','on');axis equal; title('Canonical Variates')

    % PCA:
    X=[];
    for c=1:nclasses
        X=horzcat(X,D{c});
    end
    [u s v]=svds([D{1} D{2}],2);
    V1=u(:,1); V2=u(:,2);

    figure;hold on
    for c=1:nclasses
        PD{c}(:,1) = (D{c})'*V1;
        PD{c}(:,2) = (D{c})'*V2;
        h = plot(PD{c}(:,1),PD{c}(:,2),marker{c},'markersize',4,'color',col(c,:));	
    end
    set(gca,'box','on');axis equal; title('PCA')
elseif(numDims == 3)
    V1=W(:,1); V2=W(:,2); V3=W(:,3);
    nclasses=length(D);
    figure; col=get(gca,'colororder');
    hold on
    for c=1:nclasses
        PD{c}(:,1) = (D{c})'*V1;
        PD{c}(:,2) = (D{c})'*V2;
        PD{c}(:,3) = (D{c})'*V3;
        h = plot3(PD{c}(:,1),PD{c}(:,2),PD{c}(:,3),marker{c},'markersize',4,'color',col(c,:));
    end
    set(gca,'box','on');axis equal; title('Canonical Variates')

    % PCA:
    X=[];
    for c=1:nclasses
        X=horzcat(X,D{c});
    end
    [u s v]=svds([D{1} D{2}],numDims);
    V1=u(:,1); V2=u(:,2); V3=u(:,3);

    figure;hold on
    for c=1:nclasses
        PD{c}(:,1) = (D{c})'*V1;
        PD{c}(:,2) = (D{c})'*V2;
        PD{c}(:,3) = (D{c})'*V3;
        h = plot3(PD{c}(:,1),PD{c}(:,2),PD{c}(:,3),marker{c},'markersize',4,'color',col(c,:));	
    end
    set(gca,'box','on');axis equal; title('PCA')
end