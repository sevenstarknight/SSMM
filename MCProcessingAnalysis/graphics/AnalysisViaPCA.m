 % PCA:
 
nclasses=length(D);
X=[];
for c=1:nclasses
    X=horzcat(X,D{c});
end
[u s v]=svds([D{1} D{2}],numDims);
V1=u(:,1); V2=u(:,2); V3=u(:,3);
col=get(gca,'colororder');

figure;hold on
for c=1:nclasses
    PD{c}(:,1) = (D{c})'*V1;
    PD{c}(:,2) = (D{c})'*V2;
    PD{c}(:,3) = (D{c})'*V3;
    h = plot3(PD{c}(:,1),PD{c}(:,2),PD{c}(:,3),marker{c},'markersize',4,'color',col(c,:));	
end
set(gca,'box','on');axis equal; title('PCA')