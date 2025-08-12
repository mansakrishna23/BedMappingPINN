% region or glacier name
glacier = 'Upernavik';
%required Data
accpath   = ['/totten_1/ModelData/Greenland/RACMO2Accumulation/SMBGreenland/smb_1961-1990.mat'];
velpath   = ['/totten_1/ModelData/Greenland/VelMouginot/RignotGreenland2012Vel.mat'];

fprintf(['\n\033[42m  ========  ' upper(glacier) '  ========  \033[0m\n\n']);
% Extract data for this model domain{{{
%Domain Boundary
disp('Getting domain boundaries');
domain=['./' glacier '/' glacier '.exp'];
offset=20*10^3;
A=expread(domain);
domainxlim = [min(A(1).x)-offset max(A(1).x)+offset];
domainylim = [min(A(1).y)-offset max(A(1).y)+offset];

disp(' => Loading ice thickness data');
if strcmp(glacier,'Jakobshavn') | strcmp(glacier,'HelheimStream')
   glacier2=strrep(glacier,'Stream','');
   path = ['CReSIS-TEMP/' glacier2 '_2006_2014_Composite/flightlines/' glacier2 '_2006_2014_Composite_Flightlines.nc'];
   data = ExtractFromGPRIdataset('xlim',domainxlim,'ylim',domainylim,'hemisphere',+1,'exclude',{'CReSIS/GreenlandCReSIS_1993_2019.nc'},'includeextra',{path,1});
elseif strcmp(glacier,'Helheim') | strcmp(glacier,'Kangerdlugssuaq') | strcmp(glacier,'KangerdlugssuaqStream'),
   path = ['CReSIS-TEMP/' strrep(glacier,'Stream','') '_2006_2014_Composite/flightlines/' strrep(glacier,'Stream','') '_2006_2014_Composite_Flightlines.nc'];
   data  = ExtractFromGPRIdataset('xlim',domainxlim,'ylim',domainylim,'hemisphere',+1,'includeonly',{path,1});
   data2 = ExtractFromGPRIdataset('xlim',domainxlim,'ylim',domainylim,'hemisphere',+1);
   distance = Kriging(data.x,data.y,ones(size(data.x)),data2.x,data2.y,'output','distance','treetype',2,'searchradius',1000000);
   pos = find(distance>100);
   data.x=[data.x;data2.x(pos)]; data.y=[data.y;data2.y(pos)]; data.thickness=[data.thickness;data2.thickness(pos)]; data.source = [data.source;data2.source(pos)];
   clear data2;
else
   %Standard...
   data = ExtractFromGPRIdataset('xlim',domainxlim,'ylim',domainylim,'hemisphere',+1);
end
x = data.x; y=data.y; thickness = data.thickness; source = data.source;

disp(' ');
disp(' => Remove positive thickness for regions where there is no ice');
mask = interpBedmachineGreenland(x,y,'mask');
pos = find(mask>1); %0: ocean, 1: ice-free land, 2: grounded ice, 3: floating ice, 4: non-Greenland ice
x=x(pos); y=y(pos); thickness=thickness(pos); source=source(pos);

disp(' => Adding points with zero thickness for BC');
[X,Y] = meshgrid(domainxlim(1):100:domainxlim(2),domainylim(1):100:domainylim(2));
disp('    loading Morlighem''s mask');
   pos=plotboxpos();
   pos=plotboxpos();

mask = interpBedmachineGreenland(x,y,'mask');
pos=find(mask==1);

x=[x;X(pos)];
y=[y;Y(pos)];
thickness=[thickness;zeros(length(pos),1)];
source   =[source;   zeros(length(pos),1)];

disp('Saving');
x=double(x);
y=double(y);
thickness=double(thickness);
source=double(source);

% remove points outside domain for training
numpoints=length(x);
tracks=double([x y [1:1:numpoints]']);

disp('Removing points outside of the domain');
A=expread(domain);
flags=ContourTest(A(1).x,A(1).y,tracks(:,1),tracks(:,2));
for i=2:length(A),
   flags=(flags & ~ContourTest(A(i).x,A(i).y,tracks(:,1),tracks(:,2)));
end
tracks=tracks(find(flags),:);

% dealing with x and y
[xvals, xidx] = intersect(x, tracks(:, 1), 'stable');
[yvals, yidx] = intersect(y, tracks(:, 2), 'stable');
[idx, idxx] = intersect(xidx, yidx, 'stable');

x = x(idx);
y = y(idx);
thickness = thickness(idx);

save(['./' glacier '/ProcessedTracks'],'x','y','thickness','source');
%}}}

%Create Mesh that follows data points{{{
domain=['./' glacier '/' glacier '.exp'];
hmin=400;
hmax=500;

disp('Reading tracks');
load([glacier '/ProcessedTracks.mat']);
numpoints=length(x);
tracks=double([x y [1:1:numpoints]']);

disp('Removing points outside of the domain');
A=expread(domain);
flags=ContourTest(A(1).x,A(1).y,tracks(:,1),tracks(:,2));
for i=2:length(A),
   flags=(flags & ~ContourTest(A(i).x,A(i).y,tracks(:,1),tracks(:,2)));
end
tracks=tracks(find(flags),:);

disp('Coarsening observations for mesh');
tic
pos=GroupObs(tracks(:,1),tracks(:,2),tracks(:,3),hmin);
tracks=tracks(pos,:);
toc
clear x y thickness;

disp('Generating first constrained mesh');
md = bamg(model(), 'domain',domain,'hmax',hmax,'RequiredVertices',tracks);

disp('Interpolating velocities for metric');
[vx vy] = interpRignot2012(md.mesh.x,md.mesh.y);
vel=sqrt(vx.^2+vy.^2);
clear vx vy

disp('Remeshing with metric');
hVertices=NaN*ones(md.mesh.numberofvertices,1);
hVertices(find(md.private.bamg.mesh.Vertices(:,3)>3))=hmin;
md=bamg(md,'hmax',hmax,'hmin',hmin,'hVertices',hVertices,'field',vel,'err',0.8,'grad',1.1,'maxnbv',3.*10^6);
[md.mesh.lat,md.mesh.long]=xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
md.mesh.epsg=3413;
md.miscellaneous.name=glacier;
%}}}
%Parameterize model so that it has the required features for the PINN{{{
load([glacier '/ProcessedTracks.mat']);
bedmachine=interpBedmachineGreenland(md.mesh.x,md.mesh.y,'thickness');
% bedmachine v5 topography for testing
md.geometry.bed=interpBedmachineGreenland(md.mesh.x, md.mesh.y, 'bed');
md.geometry.surface = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'surface');
rheology_B = cuffey(-8+273.15); % -8 degrees C
[md.initialization.vx, md.initialization.vy] = interpJoughinCompositeGreenland(md.mesh.x, md.mesh.y); % m/year
% get SMB from RACMO
md.smb.mass_balance =interpRACMO1km(md.mesh.x,md.mesh.y);
%Get dhdt from ICESat-2
md.balancethickness.thickening_rate =  interpSmith2020(md.mesh.x,md.mesh.y,'gris_filt');
pos1 = isnan(md.balancethickness.thickening_rate);
md.balancethickness.thickening_rate(pos1) = 0;
%}}}

% Invert for basal friction{{{
% parameterize the model
%md.transient.ismovingfront=1;

md=setflowequation(md,'SSA','all');
md=setmask(md,'','');

disp('   Interpolating mask');
mask = int8(interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask'));
md.mask.ice_levelset= -1*ones(md.mesh.numberofvertices,1);
pos = find(mask<1);
md.mask.ice_levelset(pos)=1;

disp('      reading MC bed (assumes no floating ice)');
md.geometry.bed  = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'bed');
md.geometry.base = md.geometry.bed;

disp('      reading Howat surface');
md.geometry.surface=interpBedmachineGreenland(md.mesh.x,md.mesh.y,'surface');
pos = find(md.mask.ice_levelset>0);
md.geometry.surface(pos) = md.geometry.base(pos)+10; %Minimum thickness

md.geometry.thickness = md.geometry.surface - md.geometry.bed;
pos=find(md.geometry.thickness<=10);
md.geometry.surface(pos) = md.geometry.base(pos)+10; %Minimum thickness
md.geometry.thickness = md.geometry.surface - md.geometry.bed;

md.masstransport.min_thickness = 10;

disp('   Adjusting ice mask');
%Tricky part here: we want to offset the mask by one element so that we don't end up with a cliff at the transition
pos = find(max(md.mask.ice_levelset(md.mesh.elements),[],2)>0);
md.mask.ice_levelset(md.mesh.elements(pos,:)) = 1;
% For the region where surface is NaN, set thickness to small value (consistency requires >0)
pos=find((md.mask.ice_levelset<0).*(md.geometry.surface<0));
md.mask.ice_levelset(pos)=1;
pos=find((md.mask.ice_levelset<0).*(isnan(md.geometry.surface)));
md.mask.ice_levelset(pos)=1;

disp('      -- reconstruct thickness');
md.geometry.thickness=md.geometry.surface-md.geometry.base;

disp('      reading velocities ');
[md.inversion.vx_obs md.inversion.vy_obs]=interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
pos=find(isnan(md.inversion.vx_obs) | isnan(md.inversion.vy_obs));
md.inversion.vx_obs(pos)=0;
md.inversion.vy_obs(pos)=0;
md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);
md.initialization.vx  = md.inversion.vx_obs;
md.initialization.vy  = md.inversion.vy_obs;
md.initialization.vz  = zeros(md.mesh.numberofvertices,1);
md.initialization.vel = md.inversion.vel_obs;

disp('   Initialize basal friction using driving stress');
disp('   -- Smooth the ice surface with 20 L2 projections and then compute the surface slopes');
asurf    = averaging(md,md.geometry.surface,20); % maybe executing 20 L2 projection is ok
[sx, sy, s] = slope(md, asurf);
%[sx,sy,s]= slope(md,asurf); % slope 's' comes on elements
sslope   = averaging(md,s,1); % average the slope once on the vertices, because 's' comes on elements, we need this data on vertices
disp('   -- Set the lower bound of velocity, pressure and friction coefficient');
min_velocity = 0;
min_pressure = 0;
min_friction_coef = 0.00;

disp('   -- Process surface velocity data');
vel      = md.inversion.vel_obs;
flags    = (vel==0).*(md.mask.ice_levelset<0); % interpolate on the ice parts
pos1     = find(flags);
pos2     = find(~flags);
vel(pos1)= griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1)); % interpolating the velocities where vel==0
vel      = max(vel,min_velocity); % setting minimum velocity value

disp('   -- Calculate effective pressure and the initial pressure');
Neff                       = (md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base)*md.constants.g;
Neff(find(Neff<=0))        = min_pressure; % setting minimum positve pressure
md.initialization.pressure = md.materials.rho_ice*md.geometry.thickness*md.constants.g; % setting the initial pressure
p = 3;
q = 3;

disp(['   -- Deduce friction coefficient from driving stress, with p=',num2str(p), ', q=', num2str(q)]);
md.friction.p           = p*ones(md.mesh.numberofelements,1);
md.friction.q           = q*ones(md.mesh.numberofelements,1);
driving_stress          = md.materials.rho_ice*md.constants.g*md.geometry.thickness.*(sslope);
md.friction.coefficient = sqrt(driving_stress./((Neff.^(q/p)).*(vel/md.constants.yts).^(1/p)));
md.friction.coefficient = min(md.friction.coefficient,300^(1/p));
md.friction.coupling = 2; % will be a default setting later, coupling=0 will give non-physical water pressure when above sea level.

disp('   -- Extrapolate on ice free regions (using griddata)');
flags = (md.mask.ice_levelset>0); % no ice
pos1  = find(flags);
pos2  = find(~flags);
md.friction.coefficient(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),md.friction.coefficient(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
pos = find(isnan(md.friction.coefficient) | md.friction.coefficient <=0);
md.friction.coefficient(pos)  = min_friction_coef;
% set the no ice area and negative effective pressure area to have minimal friction coef
md.friction.coefficient(pos1)  = min_friction_coef;
md.friction.coefficient(pos1)  = min_friction_coef;
md.friction.coefficient(Neff<=0)  = min_friction_coef;
%md.friction.coefficient = 50*ones(md.mesh.numberofvertices,1); %q = 1.

% flow law
disp('   Creating flow law parameters (assume ice is at -5°C for now)');
md.materials.rheology_n = 3*ones(md.mesh.numberofelements,1);
md.materials.rheology_B = cuffey(273.15 - 8)*ones(md.mesh.numberofvertices,1);
disp('   Loading accumulation rates from RACMO');
load(accpath);
md.smb.mass_balance = InterpFromGridToMesh(x_m,y_m,accumulation,md.mesh.x,md.mesh.y,-99999);
clear x_m y_m accumulation
md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);

disp('   Geothermal flux from Shapiro et al.');
md.basalforcings.geothermalflux=interpSeaRISE(md.mesh.x,md.mesh.y,'bheatflx');

disp('   Setting up thermal model');
md.initialization.temperature=min(0,interpSeaRISE(md.mesh.x,md.mesh.y,'surftemp'))+273.15;
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
md.thermal.spctemperature=md.initialization.temperature;
md.thermal.isenthalpy=1;
md.thermal.isdynamicbasalspc=1;

%Deal with boundary conditions:
disp('   Set Boundary conditions');
md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);
pos=find((md.mask.ice_levelset<0).*(md.mesh.vertexonboundary));
md.stressbalance.spcvx(pos)=md.initialization.vx(pos);
md.stressbalance.spcvy(pos)=md.initialization.vy(pos);
md.stressbalance.spcvz(pos)=0;

% Set the friction law to Weertman
md.friction=frictionweertman();
md.friction.m = 3.0*ones(md.mesh.numberofelements,1);
md.friction.C = 2000*ones(md.mesh.numberofvertices,1);

%No friction on PURELY ocean element
pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
flags=ones(md.mesh.numberofvertices,1);
flags(md.mesh.elements(pos_e,:))=0;
md.friction.C(find(flags))=0.0;

%Control general
md.inversion=m1qn3inversion(md.inversion);
md.inversion.iscontrol=1;
md.verbose=verbose('solution',false,'control',true);
md.transient.amr_frequency = 0;

%Cost functions
costcoeffs = [100, 1, 1e-10];
md.inversion.cost_functions=[101 103 501];
md.inversion.cost_functions_coefficients=zeros(md.mesh.numberofvertices,numel(md.inversion.cost_functions));
md.inversion.cost_functions_coefficients(:,1)=costcoeffs(1);
md.inversion.cost_functions_coefficients(:,2)=costcoeffs(2);
md.inversion.cost_functions_coefficients(:,3)=costcoeffs(3);
pos=find(md.mask.ice_levelset>0);
md.inversion.cost_functions_coefficients(pos,1:2)=0;

%Controls
md.inversion.control_parameters={'FrictionC'};
md.inversion.maxsteps=400;
md.inversion.maxiter =400;
md.inversion.min_parameters=1e-6*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters=5e4*ones(md.mesh.numberofvertices,1);
md.inversion.control_scaling_factors=1;
md.inversion.dxmin = 1e-6;
%Additional parameters
%md.settings.solver_residue_threshold = 1e-3; % added extra line
md.stressbalance.restol=1e-6; %1e-6 earlier
md.stressbalance.reltol=1e-5;
md.stressbalance.abstol=NaN;

md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
%Go solve
%cluster=generic('name','totten','np',20);
cluster=generic('name',oshostname(),'np', 40);
md.cluster=cluster;
md=solve(md,'sb');

% grab results
X = md.mesh.x;
Y = md.mesh.y;
C = md.results.StressbalanceSolution.FrictionC;
%}}}

% save the model
md.friction.C = C;
disp(size(md.mesh.x));

% save PINN data as ISSM model
saveasstruct(md, [glacier '/' lower(glacier) '_md.mat']);
