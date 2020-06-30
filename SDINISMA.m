% Runs the semi-discretisation method in the normalised time doamin for a
% standard 1DoF machining system and one with a semi-active inerter,
% parameterised with RA and RM [Static version to accompany ISMA 2020 paper]

%For both cases
k=240; %Approximation parameter, divisions per shorter period
N=k; %For case with no inerter only
zeta=0.005; %damping
Os = linspace(0.05,0.5,100); %Nondimensional spindle speeds
DoCs = linspace(0,0.15,50);  %Nondimensional depths of cut

%These ones for the case with inerter only
RM=0.2; %Modulation ratio 2piOm/60w_m
RA=0.1; %Amplitude m1/m0
%% Preallocation
A=[0 1; 0 0]; %Preallocate A matrix(EQUATION 9)
B=zeros(2,2); %Preallocate B matrix

D=eye(N-1); %identity matrix part of D
D=[zeros(3,N-1);D];D=[zeros(N+2,2),D,zeros(N+2,1)];D(3,1)=1; %concatenate on two empty rows for P,R, empty colums at start and end, add in the extra 1 in first column

%% Without inerter
beta=-2*zeta; %precalculation of beta

ss0=zeros(length(Os),length(DoCs)); %preallocation
doc0=zeros(length(Os),length(DoCs)); %preallocation
ei0=zeros(length(Os),length(DoCs)); %preallocation

for OmegaCount=1:length(Os)
    Omega = Os(OmegaCount); %Define local spindle speed
    ['Checking Omega=' num2str(Omega) '/' num2str(Os(end)) ' (1)'] %#ok<NOPTS>
    tau = 2*pi./Omega; %TODO: check the 2pi
    dt=tau./N; %timestep
    for wCount = 1:length(DoCs)
        w = DoCs(wCount);
        gamma=w; %precalculation of gamma
        alpha=-(1+w); %precalculation of alpha
        A(2,1)=alpha;
        A(2,2)=beta;
        Fi=eye(N+2); %preallocate Floquet matrix %CHECK THIS PADDING
        for i=1:N
            B(2,1)=gamma; %fill in A,B matrices
            P=expm(A*dt);R=(expm(A*dt)-eye(2))*inv(A)*B; %create P,R matrices as per Insperger04
            D(1:2,1:2)=P;D(1:2,N+2)=R(1:2,1); %coefficient matrix
            Fi=D*Fi; %Update floquet matrix
            D(1:2,1:2)=0;D(1:2,N+2)=0; %reset coefficient matrix
        end
        ss0(OmegaCount,wCount)=Omega; %matrix of spindle speeds
        doc0(OmegaCount,wCount)=w; %matrix of depth of cuts
        ei0(OmegaCount,wCount)=max(abs(eig(Fi))); %matrix of eigen values
    end
    
end


%% With inerter
if RM>=1 %Modulation period is longer than rotation period
    N=k/RM; %divisions per longer period (For case with inerter)
    
    D=eye(N-1); %identity matrix part of D
    D=[zeros(3,N-1);D];D=[zeros(N+2,2),D,zeros(N+2,1)];D(3,1)=1; %concatenate on two empty rows for P,R, empty colums at start and end, add in the extra 1 in first column
    
    ss=zeros(length(Os),length(DoCs)); %preallocation
    doc=zeros(length(Os),length(DoCs)); %preallocation
    ei=zeros(length(Os),length(DoCs)); %preallocation
    
    for OmegaCount=1:length(Os)
        Omega = Os(OmegaCount); %Define local spindle speed
        ['Checking Omega=' num2str(Omega) '/' num2str(Os(end))] %#ok<NOPTS>
        tau = 2*pi./Omega; %TODO: check the 2pi
        T=RM*tau; %Nondimensonalised major period is modulation  period
        dt=T./k; %timestep
        t=0:dt:T-dt; %time array from 0 to N-1
        chi=1./(1+RA.*cos(Omega.*t./RM)); %denominator of mass
        beta=-2*zeta*chi; %beta values for this spindle speed
        for wCount = 1:length(DoCs)
            w = DoCs(wCount);
            alpha=-(1+w)*chi; %precalculation of alpha
            gamma=w*chi; %precalculation of gamma
            Fi=eye(N+2); %preallocate Floquet matrix
            for i=1:k
                %fill in A,B matrices for each timestep
                A(2,1)=alpha(i);
                A(2,2)=beta(i);
                B(2,1)=gamma(i);
                P=expm(A*dt);R=(expm(A*dt)-eye(2))*inv(A)*B; %create P,R matrices as per Insperger04
                D(1:2,1:2)=P;D(1:2,N+2)=R(1:2,1); %coefficient matrix
                Fi=D*Fi; %Update floquet matrix
                D(1:2,1:2)=0;D(1:2,N+2)=0; %reset coefficient matrix
            end
            ss(OmegaCount,wCount)=Omega; %matrix of spindle speeds
            doc(OmegaCount,wCount)=w; %matrix of depth of cuts
            ei(OmegaCount,wCount)=max(abs(eig(Fi))); %matrix of eigen values
        end
    end
    
else %%for RM<1, Modulation period is shorter than rotation period 
    
    N=k; %T=tau    
    m=RM*k; 
    
    D=eye(N-1); %identity matrix part of D
    D=[zeros(3,N-1);D]; D=[zeros(N+2,2),D,zeros(N+2,1)];D(3,1)=1; %concatenate on two empty rows for P,R, empty colums at start and end, add in the extra 1 in first column
    
    ss=zeros(length(Os),length(DoCs)); %preallocation
    doc=zeros(length(Os),length(DoCs)); %preallocation
    ei=zeros(length(Os),length(DoCs)); %preallocation
    
    for OmegaCount=1:length(Os)
        Omega = Os(OmegaCount); %Define local spindle speed
%       ['Checking Omega=' num2str(Omega) '/' num2str(Os(end)) ' (' num2str(a*b+1) ')' ] %#ok<NOPTS>
        tau = 2*pi./Omega; %TODO: check the 2pi
        T=tau; %Nondimensonalised major period is just tau if tau>T_m
        dt=T./k; %timestep
        t=0:dt:T-dt; %time array from 0 to N-1
        chi=1./(1+RA.*cos(Omega.*t./RM)); %denominator of mass
        beta=-2*zeta*chi; %beta values for this spindle speed
        for wCount = 1:length(DoCs)
            w = DoCs(wCount);
            alpha=-(1+w)*chi; %precalculation of alpha
            gamma=w*chi; %precalculation of gamma
            Fi=eye(N+2); %preallocate Floquet matrix
            
            
            for i=1:m
                %This loop finds the floquet matrix across a single
                %rotationmodulation period
                
                %fill in A,B matrices for each timestep
                A(2,1)=alpha(i);
                A(2,2)=beta(i);
                B(2,1)=gamma(i);
                P=expm(A*dt);R=(expm(A*dt)-eye(2))*inv(A)*B; %create P,R matrices as per Insperger04
                D(1:2,1:2)=P;D(1:2,N+2)=R(1:2,1); %coefficient matrix
                Fi=D*Fi; %Update floquet matrix
                D(1:2,1:2)=0;D(1:2,N+2)=0; %reset coefficient matrix
            end
            %True floquet matrix is single cycle multiplied by
            %itself k/N=1/RM times
            Fi=Fi^(1/RM);
            
            ss(OmegaCount,wCount)=Omega; %matrix of spindle speeds
            doc(OmegaCount,wCount)=w; %matrix of depth of cuts
            ei(OmegaCount,wCount)=max(abs(eig(Fi))); %matrix of eigen values
        end
    end
    
end
hold on
contour(ss,doc,ei,[1 1],'r-')
 
%% Uncomment for plotting

% figure
% hold on
% contour(ss0,doc0,ei0,[1 1],'k--','displayname','Constant mass')
% contour(ss,doc,ei,[1 1],'r-','DisplayName',['$R_A=$' num2str(RA) ', $R_M=$' num2str(RM)] )
% plot(ss(ei<=1),doc(ei<=1),'go', 'displayname', 'Stable') % to show actual
% values studied
% plot(ss(ei>1),doc(ei>1),'rx', 'displayname', 'Unstable') % to show actual
% values studied
% 

%xlabel('$\tilde{\Omega}$', 'interpreter','latex')
% ylabel('$\tilde{w}$', 'interpreter','latex')
% title(['$R_A=$' num2str(RA) ', $R_M=$' num2str(RM)], 'interpreter','latex')
% hl = legend('show');
% set(hl, 'Interpreter','latex')
