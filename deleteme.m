rad = 5;
tightness = 1.1; % lattice constant = rad*tightness*2
myScale = 1;
xmax = 200*myScale;
ymax = 200*myScale;
R_grain = 50*myScale;
buffer = rad*tightness; % increasing buffer zone seems to help w border madness
left_lattice = make_hcp(rad*tightness,xmax/2-2,ymax/2+2,ceil(xmax/(2*rad*tightness)),15*pi/180); 
scatter(left_lattice(:,1),left_lattice(:,2),300,'.','MarkerEdgeColor','#000000');
