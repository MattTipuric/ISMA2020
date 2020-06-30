%{
Time domain model to verify semi-deiscretisation results. Detects chatter
through loss of contact (i.e. force value vbecomes zero) so may overpredict
stability at boundary if not enough revolutions studied. [Version for publication for ISMA]
%}

%% parameters
zeta=0.005; %Damping ratio

hm=1e-6; %Feed per rev (Arbitrary, shouldn't affect dynamics)

SR=60; %Steps per revolution

Revs=1500; %	Number of revolutions to simulate
N=SR*Revs; %total number of steps

% close all
Os = linspace(0,1,100); %Nondimensionalised spindle speed
DoCs = linspace(0,0.5,75); %Nondimensionalised Depth of cut
US=zeros(length(Os),length(DoCs));

%% System with inerter
RA=0.1; %Amplitude ratio
RM=5; %Modulation ratio

for a=1:length(Os)
    num2str(a) %To check progress
    Om=Os(a);
    dt=2*pi/(SR*Om); %Length of timestep
    for b=1:length(DoCs)
        num2str(b)
        w=DoCs(b);
        chk=0;
        x=zeros(1,N);xdot=zeros(1,N);xddot=zeros(1,N); %Prealocation and zero conditions
        z=zeros(1,N); %Preallocation of normed displacement
        chi=zeros(1,N); %Preallocation of varying mass term
       
        xmin=ones(1,SR)*hm; %inital x_min is just depth of cut as initial surface displacement=0
        
        for n=SR+1:N %for each rotation after the first
            i=rem(n,SR)+1; %index of x_min
            
            z(n)=w*(xmin(i)-x(n-1)); %Nondimensionalised cutting depth (equivalent to force term in dimensional verasion)
            if    z(n)<0
                % this lets us know the first time the system loses contact
                if chk==0

                    US(a,b)=1; %Marks this combination as unstable, stops loop
                    break
                end
            end
            
            chi(n)=1/(1+RA*cos(Om*n*dt./RM)); %varying mass term at t=n*dt 
            xddot(n)=chi(n-1)*(z(n)-2*zeta*xdot(n-1)-x(n-1)); %acceleration 
            xdot(n)=xdot(n-1)+xddot(n)*dt; %velocity 
            x(n)=x(n-1)+xdot(n)*dt; %displacement 
            xmin(i)=min([hm+x(n) hm+xmin(i)]); %update xmin
        end
    end
end

[x,y,~]=find(US==0)

figure
plot(Os(x),DoCs(y),'.')
xlabel('$\tilde{\Omega}$','interpreter','latex')
ylabel('$\tilde{w}$','interpreter','latex')






