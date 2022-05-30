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
more off
clear
directory='';
load([directory,'Nodes_and_Struts_cleaned_15.mat'])
%matfile([directory,'Nodes_and_Struts_cleaned_15.mat']) <- content of mat file.

%%% pretraitement des struts
strutsmat=cell2mat(struts);
n_struts=numel(struts);
strutsmat_unfolded=[strutsmat(:,1);strutsmat(:,2)];

%%%% pretraitement des vertex
vertices_alive=unique(strutsmat_unfolded);
n_vertices=numel(vertices_alive);
disp(['nombre de vertex=',num2str(n_vertices)])

verticesraw=cell2mat(vertices_COM);
xg_v=verticesraw(vertices_alive,1);
yg_v=verticesraw(vertices_alive,2);
zg_v=verticesraw(vertices_alive,3);

xmin_v=min(xg_v);xmax_v=max(xg_v);
ymin_v=min(yg_v);ymax_v=max(yg_v);
zmin_v=min(zg_v);zmax_v=max(zg_v);


%%%% Strut lengths histogram
strutlengths=zeros(1,n_struts);
for j=1:n_struts
    j1=struts{j}(1);j2=struts{j}(2);
    strutlengths(j)=sqrt((vertices_COM{j1}(1)-vertices_COM{j2}(1))^2+...
(vertices_COM{j1}(2)-vertices_COM{j2}(2))^2+...
(vertices_COM{j1}(3)-vertices_COM{j2}(3))^2);
end
av_st_length=mean(strutlengths);
std_st_length=std(strutlengths);

fig_slh=figure;
histogram(strutlengths/av_st_length,'Normalization','pdf');
xlabel('$\ell/\langle\ell\rangle$=normalized strut length','Interpreter','latex');
ylabel('proba($\ell/\langle\ell\rangle$)','Interpreter','latex');
an=annotation('textbox','fitboxtotext','on');
an.String={['Average=',num2str(av_st_length)],['Std dev=',num2str(std_st_length)]};
an.Position=[0.4929    0.7595    0.2957    0.1296];
an.FitBoxToText='on';
drawnow
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% fluctuations de l'angle dièdre (écart aux 109 degres..)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
les109degres=acos(-1/3)/pi*180;
binning=30;
xmv=linspace(0.9*xmin_v,1.1*xmax_v,binning);
ymv=linspace(0.9*ymin_v,1.1*ymax_v,binning);
zmv=linspace(0.9*zmin_v,1.1*zmax_v,binning);


density109=zeros(binning,binning+1,binning);
density109_alongz=zeros(binning,2);
densityv_alongz=zeros(binning,1);
density_anomalousv_alongz=zeros(binning,1);

for j=1:numel(xg_v)
   % if(mod(j,floor(numel(xg_v)/10))==0) disp(num2str(j/n_struts));end
    ux=find(xmv>=xg_v(j),1);
    uy=find(ymv>=yg_v(j),1); 
    uz=find(zmv>=zg_v(j),1);
    
    vertexabsoluteindex=vertices_alive(j);
	tmplines=find(strutsmat_unfolded==vertexabsoluteindex);
	vertex_order=nnz(tmplines);
	if(vertex_order~=4)
		density_anomalousv_alongz(uz)=density_anomalousv_alongz(uz)+1;
	else
		angle=0;angle2=0;
		tmplines(tmplines>n_struts)=tmplines(tmplines>n_struts)-n_struts;
		for k=1:3
			j1=struts{tmplines(k)}(1);j2=struts{tmplines(k)}(2);
			if (j2==vertexabsoluteindex) j2=j1;end
			%on est obligé d'utiliser vertices_COM et non pas xg_v...
			s1x=vertices_COM{j2}(1)-vertices_COM{vertexabsoluteindex}(1);
			s1y=vertices_COM{j2}(2)-vertices_COM{vertexabsoluteindex}(2);
			s1z=vertices_COM{j2}(3)-vertices_COM{vertexabsoluteindex}(3);
			for kk=k+1:4
				j3=struts{tmplines(kk)}(1);j4=struts{tmplines(kk)}(2);
				if (j4==vertexabsoluteindex) j4=j3;end
				s2x=vertices_COM{j4}(1)-vertices_COM{vertexabsoluteindex}(1);
				s2y=vertices_COM{j4}(2)-vertices_COM{vertexabsoluteindex}(2);
				s2z=vertices_COM{j4}(3)-vertices_COM{vertexabsoluteindex}(3);
				lecos=(s1x*s2x+s1y*s2y+s1z*s2z)/sqrt(s1x^2+s1y^2+s1z^2)/sqrt(s2x^2+s2y^2+s2z^2);
				if(abs(lecos)>1) lecos=sign(lecos)*min(abs(lecos),1);end
				angle=angle+acos(lecos)/pi*180;
				angle2=angle2+(acos(lecos)/pi*180)^2;
			end
		end
		density109(ux,uy,uz)=density109(ux,uy,uz)+sqrt(angle2/6-(angle/6)^2)/(angle/6);
		density109_alongz(uz,1)=density109_alongz(uz,1)+angle/6;
		density109_alongz(uz,2)=density109_alongz(uz,2)+sqrt(angle2/6-(angle/6)^2)/(angle/6);		
		densityv_alongz(uz)=densityv_alongz(uz)+1;
	end
end
%%% visu
f=figure;
set(f,'position',[ 23        1616        1000         438]);
subplot(2,1,1)
plotyy(zmv,densityv_alongz,zmv,density_anomalousv_alongz);
title('raw counting of order 4 and order <>4 vertices');
xlabel('z');
subplot(2,1,2)
plotyy(zmv,density109_alongz(:,1)./max(1,densityv_alongz)/les109degres,zmv,density109_alongz(:,2)./max(1,densityv_alongz))
title('local mean dihedral angle and standard deviation')
xlabel('z');

