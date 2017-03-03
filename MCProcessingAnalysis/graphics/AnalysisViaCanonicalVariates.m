D{1} = dataOut(dataOut(:,69) == 1,1:68)';
D{2} = dataOut(dataOut(:,70) == 1,1:68)';
D{3} = dataOut(dataOut(:,71) == 1,1:68)';
D{4} = dataOut(dataOut(:,72) == 1,1:68)';

figure
marker={'+','o','d', '^'};

numDims = 3;
W=CanonVar(D,3);


V1=W(:,1); V2=W(:,2); V3=W(:,3);
nclasses=length(D);

for c=1:nclasses
    PD{c}(:,1) = (D{c})'*V1;
    PD{c}(:,2) = (D{c})'*V2;
    PD{c}(:,3) = (D{c})'*V3;
end

scatter3(PD{1}(:,1), PD{1}(:,2), PD{1}(:,3), 4, 'r')
hold on
scatter3(PD{2}(:,1), PD{2}(:,2), PD{2}(:,3), 4, 'b')
scatter3(PD{3}(:,1), PD{3}(:,2), PD{3}(:,3), 4, 'g')
% scatter3(PD{4}(:,1), PD{4}(:,2), PD{4}(:,3), 4, 'k')
xlabel('Canonical Variates 1')
ylabel('Canonical Variates 2')
zlabel('Canonical Variates 3')
hold off
legend('SXPh', 'RRab', 'EB')
daspect([1,1,.3]);axis tight; 
OptionZ.FrameRate=15;OptionZ.Duration=7.5;OptionZ.Periodic=true; 
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'CanonicalVariates',OptionZ)