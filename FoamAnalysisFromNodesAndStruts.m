%%%%% Foam Topology Analysis
%%%%% From Tomography Stack Images.
%%%%% Copyright ICS - 2022

%% structure des données à traiter (b)

% 1) vertices_COM 
%%% vecteur de type "cell" 
%%% vertices_COM{j} contient un vecteur [x y z] contenant les coordonnees (en pixels) du j-ieme vertex
%%% ATTENTION : une grande part de ces vertex sont fantômes suite au processus d'extraction, ce qui est traité dans
%%% la section "pretraitement des vertex".

% 2) struts
%%% vecteur de type "cell"
%%% struts{j} contient un vecteur [v1 v2] où v1 et v2 sont deux index de vertex dans le tableau vertices_COM
%%% Contrairement aux vertex, tous les struts presents sont signifiants

%% structure des données à traiter (e)

%%%%%%%%%%%
%%%% Loading (b)
%%%%%%%%%%%
%{
f                  : Name of the file where the data is.
p                  : True path to the file f.
true name          : The name of the file f without its extension.
struts             : "cell" type where it connects vertex index A to vertex index B.
vertices_COM       : "cell" type where it has [x y z] component of the 3D position of the vertex.
strutsmat          : struts but in a double-column array instead of cell.
nstruts            : total number of struts.
strutsmat_unfolded : strutsmat but with the two columns following each-others.
vertices_alive     : alive vertices.
n_vertices         : total number of vertices.
verticesraw        : vertices_COM but in a three-column array instead of cell.
xg_v               : x coordinate column of verticesraw.
yg_v               : y coordinate column of verticesraw.
zg_v               : z coordinate column of verticesraw.
xmin_v             : minimum x value coordinate of the data set.
xmax_v             : maximum x value coordinate of the data set.
ymin_v             : minimum y value coordinate of the data set.
ymax_v             : maximum y value coordinate of the data set.
zmin_v             : minimum z value coordinate of the data set.
zmax_v             : maximum z value coordinate of the data set.
%}
more off
clear

[f, p]     = uigetfile('*.mat');
eval(sprintf('cd %s',p));
true_name  = regexprep(f,'.mat','','ignorecase');
mkdir(true_name);
filename   = [p f];
fileID     = fopen(filename,'r');
directory  = '';
load(f)

%% For SlicesRight
% This section has to be un-commented if the nodes and struts extraction
% has been made from a "SlicesRight" tomography stack.
% It is so because we need to arrange-back the x,y,z directions following
% the physical ones.

for j = 1:numel(vertices_COM)
    cache = [vertices_COM{j,1}(1),vertices_COM{j,1}(2),vertices_COM{j,1}(3)];
    vertices_COM{j,1}(1) = cache(3);
    vertices_COM{j,1}(2) = cache(1);
    vertices_COM{j,1}(3) = cache(2);
end

%% Struts pre-treatment
strutsmat          = cell2mat(struts);
n_struts           = numel(struts);
strutsmat_unfolded = [strutsmat(:,1);strutsmat(:,2)];

%% Vertex pre-treatment
vertices_alive     = unique(strutsmat_unfolded);
n_vertices         = numel(vertices_alive);
disp(['nombre de vertex=',num2str(n_vertices)])

%%
verticesraw = cell2mat(vertices_COM);
xg_v        = verticesraw(vertices_alive,1);
yg_v        = verticesraw(vertices_alive,2);
zg_v        = verticesraw(vertices_alive,3);
xmin_v      = min(xg_v);
xmax_v      = max(xg_v);
ymin_v      = min(yg_v);
ymax_v      = max(yg_v);
zmin_v      = min(zg_v);
zmax_v      = max(zg_v);

disp('The dimension (in px) of the sample is');
disp([max(verticesraw(:,1)), max(verticesraw(:,2)), max(verticesraw(:,3))])

%% Strut lengths histogram
%{
strutlengths  : length of the struts.
av_st_length  : average value of the struts length.
std_st_length : standard deviation of the struts length values.
%}
strutlengths = zeros(1,n_struts);
for j=1:n_struts
    j1=struts{j}(1);j2=struts{j}(2);
    strutlengths(j)=sqrt((vertices_COM{j1}(1)-vertices_COM{j2}(1))^2+...
(vertices_COM{j1}(2)-vertices_COM{j2}(2))^2+...
(vertices_COM{j1}(3)-vertices_COM{j2}(3))^2);
end

av_st_length  = mean(strutlengths);
std_st_length = std(strutlengths);
fig_slh       = figure('visible','off');
histogram(strutlengths/av_st_length,'Normalization','pdf');
xlabel('$\ell/\langle\ell\rangle$=normalized strut length','Interpreter','latex');
ylabel('proba($\ell/\langle\ell\rangle$)','Interpreter','latex');
an=annotation('textbox','fitboxtotext','on');
an.String={['Average=',num2str(av_st_length)],['Std dev=',num2str(std_st_length)]};
an.Position=[0.4929    0.7595    0.2957    0.1296];
an.FitBoxToText='on';
drawnow
exportgraphics(fig_slh,strcat(true_name,'/strut-length-distribution.png'));

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Dihedral angle fluctuations (109.5 degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
les109degres                    : numerical value : 109.47 degrees.
binning                         : number of bins to discretize axis onto.
xmv                             : discretization of the x axis with linear spacing onto the "binning" value defined before bins. 
ymv                             : discretization of the y axis with linear spacing onto the "binning" value defined before bins.
zmv                             : discretization of the z axis with linear spacing onto the "binning" value defined before bins. 
density109_alongz               : dihedral angle distribution of vertex of order equal to 4.
density109_anomalousv_alongz    : dihedral angle distribution of vertex of order different than 4.
densityv_alongz                 : density of vertices in a binned z of order equal to 4.
density_anomalousv_alongz       : density of vertices in a binned z of order different than 4.
normalised_quantity_density     : number of vertices of order different than 4 divided by the total number of vertices (4 and different than 4) in a given binned z.
zmv.mat                         : save in .mat file extension of zmv.
normalised_quantity_density.mat : save in .mat file extension of normalised_quantity_density.
%}

les109degres                 = acos(-1/3)/pi*180;
binning                      = 30;
xmv                          = linspace(0.9*xmin_v,1.1*xmax_v,binning);
ymv                          = linspace(0.9*ymin_v,1.1*ymax_v,binning);
zmv                          = linspace(0.9*zmin_v,1.1*zmax_v,binning);
density109                   = zeros(binning,binning+1,binning);
density109_alongz            = zeros(binning,2);
density109_anomalousv_alongz = zeros(binning,2);
densityv_alongz              = zeros(binning,1);
density_anomalousv_alongz    = zeros(binning,1);

for j=1:numel(xg_v)
   % if(mod(j,floor(numel(xg_v)/10))==0) disp(num2str(j/n_struts));end
    ux                  = find(xmv>=xg_v(j),1);
    uy                  = find(ymv>=yg_v(j),1); 
    uz                  = find(zmv>=zg_v(j),1);
    vertexabsoluteindex = vertices_alive(j);
	tmplines            = find(strutsmat_unfolded==vertexabsoluteindex);
    vertex_order        = nnz(tmplines);
	if(vertex_order~=4)
        angle                       = 0;
        angle2                      = 0;
		tmplines(tmplines>n_struts) = tmplines(tmplines>n_struts)-n_struts;
		for k=1:vertex_order-1
			j1 = struts{tmplines(k)}(1);
            j2 = struts{tmplines(k)}(2);
			if (j2==vertexabsoluteindex) j2=j1;end
			%on est obligé d'utiliser vertices_COM et non pas xg_v...
			s1x = vertices_COM{j2}(1)-vertices_COM{vertexabsoluteindex}(1);
			s1y = vertices_COM{j2}(2)-vertices_COM{vertexabsoluteindex}(2);
			s1z = vertices_COM{j2}(3)-vertices_COM{vertexabsoluteindex}(3);
			for kk=k+1:vertex_order
				j3=struts{tmplines(kk)}(1);
                j4=struts{tmplines(kk)}(2);
				if (j4==vertexabsoluteindex) j4=j3;end
				s2x    = vertices_COM{j4}(1)-vertices_COM{vertexabsoluteindex}(1);
				s2y    = vertices_COM{j4}(2)-vertices_COM{vertexabsoluteindex}(2);
				s2z    = vertices_COM{j4}(3)-vertices_COM{vertexabsoluteindex}(3);
				lecos  = (s1x*s2x+s1y*s2y+s1z*s2z)/sqrt(s1x^2+s1y^2+s1z^2)/sqrt(s2x^2+s2y^2+s2z^2);
				if(abs(lecos)>1) lecos=sign(lecos)*min(abs(lecos),1);end
				angle  = angle+acos(lecos)/pi*180;
				angle2 = angle2+(acos(lecos)/pi*180)^2;
			end
		end
		density109(ux,uy,uz)               = density109(ux,uy,uz)+sqrt(angle2/6-(angle/6)^2)/(angle/6);
		density109_anomalousv_alongz(uz,1) = density109_anomalousv_alongz(uz,1)+angle/6;
		density109_anomalousv_alongz(uz,2) = density109_anomalousv_alongz(uz,2)+sqrt(angle2/6-(angle/6)^2)/(angle/6);	
        density_anomalousv_alongz(uz)      = density_anomalousv_alongz(uz)+1;
	else
		angle                       = 0;
        angle2                      = 0;
		tmplines(tmplines>n_struts) = tmplines(tmplines>n_struts)-n_struts;
		for k=1:3
			j1  = struts{tmplines(k)}(1);
            j2  = struts{tmplines(k)}(2);
			if (j2==vertexabsoluteindex) j2=j1;end
			%on est obligé d'utiliser vertices_COM et non pas xg_v...
			s1x = vertices_COM{j2}(1)-vertices_COM{vertexabsoluteindex}(1);
			s1y = vertices_COM{j2}(2)-vertices_COM{vertexabsoluteindex}(2);
			s1z = vertices_COM{j2}(3)-vertices_COM{vertexabsoluteindex}(3);
			for kk=k+1:4
				j3     = struts{tmplines(kk)}(1);j4=struts{tmplines(kk)}(2);
				if (j4==vertexabsoluteindex) j4=j3;end
				s2x    = vertices_COM{j4}(1)-vertices_COM{vertexabsoluteindex}(1);
				s2y    = vertices_COM{j4}(2)-vertices_COM{vertexabsoluteindex}(2);
				s2z    = vertices_COM{j4}(3)-vertices_COM{vertexabsoluteindex}(3);
				lecos  = (s1x*s2x+s1y*s2y+s1z*s2z)/sqrt(s1x^2+s1y^2+s1z^2)/sqrt(s2x^2+s2y^2+s2z^2);
				if(abs(lecos)>1) lecos=sign(lecos)*min(abs(lecos),1);end
				angle  = angle+acos(lecos)/pi*180;
				angle2 = angle2+(acos(lecos)/pi*180)^2;
			end
		end
		density109(ux,uy,uz)    = density109(ux,uy,uz)+sqrt(angle2/6-(angle/6)^2)/(angle/6);
		density109_alongz(uz,1) = density109_alongz(uz,1)+angle/6;
		density109_alongz(uz,2) = density109_alongz(uz,2)+sqrt(angle2/6-(angle/6)^2)/(angle/6);		
		densityv_alongz(uz)     = densityv_alongz(uz)+1;
	end
end
normalised_quantity_density = [];
for n = 1:numel(densityv_alongz)
    normalised_quantity_density(end+1) = density_anomalousv_alongz(n)/ (densityv_alongz(n)+density_anomalousv_alongz(n));
end
save('zmv.mat','zmv');
save('normalised_quantity_density.mat','normalised_quantity_density');

%% Plotting
fig_order_vrtx = figure('visible','off');
set(fig_order_vrtx,'position',[ 23        1616        1000         438]);
subplot(2,1,1)
hold on
plot(zmv,normalised_quantity_density,'x-','LineWidth',2); 
errorbar([197.7812,279.5963,369.6744,465.5747,576.5098], [0,0,0,0,0],[18.4692,18.4692,18.4692,18.4692,18.4692],'horizontal','x','LineWidth',2);
hold off
legend({'Vertex of order $\neq 4$/ Total','Fiber presence'},'Interpreter','Latex');
title('Counting of order $=4$ and order $\neq 4$ vertices','Interpreter','Latex');
xlabel('$z$ position [px]','Interpreter','Latex');
ylabel('Proportion of vertices order different of $4$','Interpreter','Latex');
subplot(2,1,2)
hold on
plot(zmv,(density109_alongz(:,1)./max(1,densityv_alongz))/les109degres  ,'x-','LineWidth',2)
plot(zmv,density109_anomalousv_alongz(:,1)./max(1,density_anomalousv_alongz)/les109degres,'x-','LineWidth',2)
%plot([min(zmv),max(zmv)], [90/les109degres, 90/les109degres])
%plot([min(zmv),max(zmv)], [mean((density109_alongz(:,1)./max(1,densityv_alongz))/les109degres ), mean((density109_alongz(:,1)./max(1,densityv_alongz))/les109degres )])
%plot([min(zmv),max(zmv)], [mean(density109_anomalousv_alongz(:,1)./max(1,density_anomalousv_alongz)/les109degres), mean(density109_anomalousv_alongz(:,1)./max(1,density_anomalousv_alongz)/les109degres)])
hold off
ylim([0,2])
title('Local mean dihedral angle distribution','Interpreter','Latex')
xlabel('$z$ position [px]','Interpreter','Latex');
ylabel('Dihedral angle normalised by $109.5^\circ$','Interpreter','Latex');

hold on
errorbar([197.7812,279.5963,369.6744,465.5747,576.5098], [0,0,0,0,0],[18.4692,18.4692,18.4692,18.4692,18.4692],'horizontal','x','LineWidth',2);
hold off
%legend({'Vertex of order $=4$','Vertex of order $\neq 4$','$\langle$ order $=4\rangle$','$\langle$ order $\neq 4\rangle$','Fiber presence'},'Interpreter','Latex');
legend({'Vertex of order $=4$','Vertex of order $\neq 4$','Fiber presence'},'Interpreter','Latex');
exportgraphics(fig_order_vrtx,strcat(true_name,'/order_vertex.png'));
disp('The total number of vertex on this plot is :');
disp(sum(densityv_alongz)+sum(density_anomalousv_alongz));
disp('The mean angle for vertex of order different from 4 is :');
disp(mean(density109_anomalousv_alongz(:,1)./max(1,density_anomalousv_alongz)/les109degres));

%%

histogram_dihedral = figure('Visible','off');
hold on
edges = linspace(0.4, 2, 11); 
histogram((density109_alongz(:,1)./max(1,densityv_alongz))/les109degres ,'BinEdges',edges);
histogram(density109_anomalousv_alongz(:,1)./max(1,density_anomalousv_alongz)/les109degres,'BinEdges',edges);
legend({'Vertex of order $=4$','Vertex of order $\neq 4$'},'Interpreter','Latex');
xlabel('Dihedral angle normalised by $109.5^\circ$','Interpreter','latex')
ylabel('Count','Interpreter','latex')
hold off
exportgraphics(histogram_dihedral,strcat(true_name,'/dihedral_distribution.png'));


%%
%{
test_x : array of x values of vertices_COM.
test_y : array of y values of vertices_COM.
test_z : array of z values of vertices_COM.
%}
test_list = zeros(numel(vertices_COM),3);
test_x = zeros(numel(vertices_COM),1);
test_y = zeros(numel(vertices_COM),1);
test_z = zeros(numel(vertices_COM),1);
for i=1:numel(vertices_COM)
    test_list(i,1) = vertices_COM{i}(1);
    test_list(i,2) = vertices_COM{i}(2);
    test_list(i,3) = vertices_COM{i}(3);
    test_x(i) = vertices_COM{i}(1);
    test_y(i) = vertices_COM{i}(2);
    test_z(i) = vertices_COM{i}(3);
end

%% Length of the struts as a function of z
%{
length_with_z           : strut length.
z_position_struts       : algebraic mean of the z position of a strut (its z center)
num_slices              : number of binning slices.
bin_slice_z             : the list of slices to do.
plot_slice_z            : array to plot slices latter on.
local_mean_length_strut : strut length in the binned z selected bin.
local_std_length_strut  : standard deviation of the local_mean_length_strut.
%}
length_with_z     = zeros(numel(struts),1);
z_position_struts = zeros(numel(struts),1);
for i = 1:numel(struts)
    length_with_z(i)     = sqrt((test_x(struts{i}(2))-test_x(struts{i}(1)))^2+(test_y(struts{i}(2))-test_y(struts{i}(1)))^2+(test_z(struts{i}(2))-test_z(struts{i}(1)))^2)  ;
    z_position_struts(i) = 0.5*(test_z(struts{i}(1))+test_z(struts{i}(2))) ;
end

num_slices    = 21;
bin_slice_z   = linspace(min(z_position_struts),max(z_position_struts),num_slices);
plot_slice_z  = [];
for k = 1:(numel(bin_slice_z)-1)
    plot_slice_z(end+1) = (bin_slice_z(k)+bin_slice_z(k+1))/2;
end
local_mean_length_strut = [];
local_std_length_strut  = [];
num_bin                 = zeros(num_slices-1,1);
for j = 1:(num_slices-1)
    cache_array  = [];
    for i = 1:numel(z_position_struts)
        if z_position_struts(i) > bin_slice_z(j) && z_position_struts(i) < bin_slice_z(j+1)
            cache_array(end+1) = length_with_z(i);
            num_bin(j) = num_bin(j)+1;
        end
    end
    local_mean_length_strut(end+1) = mean(cache_array);
    local_std_length_strut(end+1)  = std(cache_array);
end

num_bin                 = num_bin./numel(z_position_struts);
local_mean_length_strut = local_mean_length_strut.'; % .' is doing the transposition of the matrix same as transpose(A)
local_std_length_strut  = local_std_length_strut.';

%% Plot
%{
error_bar_mean : for plotting purposes.
%}
fig11 = figure('visible','off');
error_bar_mean = [] ;
for l = 1:numel(plot_slice_z)
    error_bar_mean(end+1) = plot_slice_z(1)-bin_slice_z(1);
end
hold on
errorbar(plot_slice_z,local_mean_length_strut,local_std_length_strut,local_std_length_strut,error_bar_mean,error_bar_mean,'.','LineWidth',2)
plot([min(plot_slice_z), max(plot_slice_z)],[av_st_length, av_st_length],'LineWidth',2)
legend({'$\langle l\rangle (z)$','$\langle l \rangle$'},'Location','southwest','Interpreter','Latex');
xlabel('$z$ position [pixel]','Interpreter','Latex');
hold off
exportgraphics(fig11,strcat(true_name,'/strut-length-fct-z.png'));

%% 
%{
list_pertinent_index_vertices   : the list of the vertices of order different of 0, i.e. existing vertices.
connectivity_order_vertex       : vertex order of each vertices before cleaning up the data set of the 0 order vertices.
connectivity_order_vertex_fixed : vertex order of each vertices after cleaning up the data set of the 0 order vertices.
%}
list_pertinent_index_vertices = [];
connectivity_order_vertex     = zeros(numel(vertices_COM),1);
for i = 1:numel(vertices_COM)
    cache_connectivity = 0;
    for j = 1:numel(struts)
        if struts{j,1}(1) == i
            cache_connectivity = cache_connectivity+1;
        end
        if struts{j,1}(2) == i
            cache_connectivity = cache_connectivity+1;
        end
    end
    if cache_connectivity ~= 0
        list_pertinent_index_vertices(end+1) = i;
    end
    connectivity_order_vertex(i) = cache_connectivity;
end
disp('The mean order of vertex before cleaning is :');
disp(mean(connectivity_order_vertex));

connectivity_order_vertex_fixed = zeros(numel(list_pertinent_index_vertices),1);
for i = 1:numel(list_pertinent_index_vertices)
    cache_connectivity = 0;
    for j = 1:numel(struts)
        if struts{j,1}(1) == list_pertinent_index_vertices(i)
            cache_connectivity = cache_connectivity+1;
        end
        if struts{j,1}(2) == list_pertinent_index_vertices(i)
            cache_connectivity = cache_connectivity+1;
        end
    end
    connectivity_order_vertex_fixed(i) = cache_connectivity;
end
disp('The mean order of vertex after cleaning is :');
disp(mean(connectivity_order_vertex_fixed));


%% Plot
histogram_order_vertex = figure('Visible','off');
histogram(connectivity_order_vertex);
xlabel('vertex order','Interpreter','latex');
exportgraphics(histogram_order_vertex,strcat(true_name,'/histogram_order_vertex.png'));

histogram_order_vertex2 = figure('Visible','off');
histogram(connectivity_order_vertex_fixed);
xlabel('vertex order','Interpreter','latex');
exportgraphics(histogram_order_vertex2,strcat(true_name,'/histogram_order_vertex_fixed.png'));

%% topological g(n)
%{
max_n_order         : maximum n-th neighbourhood to calculate the connectivity.
All_n_connectionn_n :
All_n_connectionn   :
topological_g_r     : array of the value of the topological g(n) for each value of n.
up_to               : making the calculation of g(n) as a mean over vertices from 1 up to "up_to". This calculation is lengthy, thus this value, but if one has great computing performance, one should use up_to = numel(list_pertinent_index_vertices).
%}
%{
max_n_order         = 7;
All_n_connectionn_n = cell(max_n_order,1);
All_n_connectionn   = cell(max_n_order,1);
topological_g_r     = zeros(max_n_order-1,1);
up_to               = 400;%numel(list_pertinent_index_vertices);
for k = 1:up_to
    index_start                  = list_pertinent_index_vertices(k);
    index                        = index_start;
    first_vertex_start           = vertices_COM{index_start,1};
    already_tagged_vertex        = {};
    already_tagged_vertex{end+1} = index_start;
    All_n_connectionn{1}         = [index];
    All_n_connectionn{2}         = find([struts{:}] == index);    
    already_done_connection      = [];
    already_done_connection      = [already_done_connection,index];
    setdif_contact               = [];
    All_n_connectionn_n{2}       = [];
    for i = 1:numel(All_n_connectionn{2})
        already_done_connection(end+1) = struts{fix(All_n_connectionn{2}(i)/2) +rem(All_n_connectionn{2}(i),2)}( rem(All_n_connectionn{2}(i),2)+1 );
        All_n_connectionn_n{2}(end+1)  = struts{fix(All_n_connectionn{2}(i)/2) +rem(All_n_connectionn{2}(i),2)}( rem(All_n_connectionn{2}(i),2)+1 );
    end
    All_n_connectionn_n{2}       = unique(All_n_connectionn_n{2});
    topological_g_r(1)           = topological_g_r(1)+numel(All_n_connectionn_n{2}); 
    setdif_contact               = [setdif_contact,All_n_connectionn_n{2}];
    for w = 2:max_n_order-1
        All_n_connectionn_n{w+1} = [];
        for i = 1:numel(All_n_connectionn_n{w})
            All_n_connectionn{w+1} = find([struts{:}] == All_n_connectionn_n{w}(i));
            for i = 1:numel(All_n_connectionn{w+1})
                already_done_connection(end+1)  = struts{fix(All_n_connectionn{w+1}(i)/2) +rem(All_n_connectionn{w+1}(i),2)}( rem(All_n_connectionn{w+1}(i),2)+1 );
                All_n_connectionn_n{w+1}(end+1) = struts{fix(All_n_connectionn{w+1}(i)/2) +rem(All_n_connectionn{w+1}(i),2)}( rem(All_n_connectionn{w+1}(i),2)+1 );
            end
        end
        All_n_connectionn_n{w+1} = unique(All_n_connectionn_n{w+1});
        All_n_connectionn_n{w+1} = setdiff(All_n_connectionn_n{w+1},[All_n_connectionn{1},setdif_contact]);
        already_done_connection  = unique(already_done_connection);
        topological_g_r(w)       = topological_g_r(w)+numel(All_n_connectionn_n{w+1});
        setdif_contact           = [setdif_contact, All_n_connectionn_n{w+1}];
    end   
end
for i = 1:numel(topological_g_r)
        (i) = (1/up_to)* topological_g_r(i)/(i^2);
end
%% Ploting
topological_g_r_fig=figure('visible','off');
hold on
plot(1:numel(topological_g_r),topological_g_r,'x-','LineWidth',2);
hold off
xlabel('number of step $n$','Interpreter','Latex');
ylabel('number of unique vertex','Interpreter','Latex');
ylim([0,4]);
title('Topological $g(n)$','Interpreter','Latex');
exportgraphics(topological_g_r_fig,strcat(true_name,'/topological_g_r.png'));
%}

%% The 3D plot of the struts


%{

[fi, rep] = uigetfile('*.tif');
chemin = fullfile(rep, '*.tif');
list = dir(chemin); 

mkdir(strcat(true_name,'/slices'));
for j = 1:1013
    fff=figure(Visible="off");
    image = fullfile(rep, list(j).name);
    I = imread(image);
    imshow(I);    
    hold on
    cache = 0 ;
    for i=1:numel(list_pertinent_index_vertices)
        if vertices_COM{list_pertinent_index_vertices(i),1}(3) == j
            %plot3([test_x(struts{i}(1)),test_x(struts{i}(2))],[test_y(struts{i}(1)),test_y(struts{i}(2))],[test_z(struts{i}(1)),test_z(struts{i}(2))])
            %plot([test_x(struts{i}(1)),test_x(struts{i}(2))],[test_y(struts{i}(1)),test_y(struts{i}(2))])
            scatter(vertices_COM{list_pertinent_index_vertices(i),1}(1),vertices_COM{list_pertinent_index_vertices(i),1}(2),'x','LineWidth',2);
            cache = cache+1;
        end
    end
    xlim([0, 1812])
    ylim([0, 1198])
    xlabel('$x$ position [pixel]','Interpreter','latex')
    ylabel('$y$ position [pixel]','Interpreter','latex')
    
    hold off
    if cache ~= 0
        name = strcat('/','slice',num2str(j),'.png');
        exportgraphics(fff,strcat(true_name,'/slices',name));
    end

end
%plot3(test_x,test_y,test_z,'x');
%}

