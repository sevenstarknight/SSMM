

S = zeros(1,length(V1));
C = zeros(1,length(V1));

S(V4 == 1) = 15;
S(V5 == 1) = 15;
S(V6 == 1) = 15;
S(V7 == 1) = 15;

C(V4 == 1) = 'r'; 
C(V5 == 1) = 'g';
C(V6 == 1) = 'b';
C(V7 == 1) = 'k';

scatter3(V1(V4 == 1), V2(V4 == 1), V3(V4 == 1), S(V4 == 1), 'r')
hold on
scatter3(V1(V5 == 1), V2(V5 == 1), V3(V5 == 1), S(V5 == 1), 'b')
scatter3(V1(V6 == 1), V2(V6 == 1), V3(V6 == 1), S(V6 == 1), 'g')
scatter3(V1(V7 == 1), V2(V7 == 1), V3(V7 == 1), S(V7 == 1), 'k')
xlabel('SPCA 1')
ylabel('SPCA 2')
zlabel('SPCA 3')
hold off
daspect([1,1,.3]);axis tight; 
OptionZ.FrameRate=15;OptionZ.Duration=7.5;OptionZ.Periodic=true; 
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'AohSPCA',OptionZ)