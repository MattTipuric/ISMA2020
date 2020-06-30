% Uses standard equations to run analytically find stability lobes for a
% stand 1DoF Machining system with discrete inertance variation

% clear
% close all
cmap=load('cmap.mat').cmap; %colour map
zeta=0.05; %damping ratio
RA=linspace(0,0.2,5); %Amplitude ratios
n=floor(length(cmap)/length(RA)); cmap=cmap(1:n:n*length(RA),:); %shortens colourmap

figure
Nmax=10; %How many lobes do you want to study?
NN=Nmax:-1:0; 

for a=1:length(RA)
    num2str(a)
    ra=RA(a); %local R_A value
    
    p=linspace(1-ra,2.5,5000); %Frequency ratio above w-wn (only interested in negtive part)
    p=p(2:end); %avoids zeo in the real part of transfer fn
    Re=((1-(1+ra).*p.^2)./((1-(1+ra).*p.^2).^2+4.*zeta.^2.*p.^2)); %real part of transfer function
    Im=((-2.*zeta.*p)./((1-(1+ra).*p.^2).^2+4.*zeta.^2.*p.^2));     %imaginary part of transfer function
    blim=-1./(2.*Re); %Limiting depth of cut
    blim(1)=blim(2);
    e=(2*pi-atan(Re./Im)); %Wave fraction
    
    
    OM=zeros(length(NN),length(p)); %spindle speeds
    for i=1:length(NN) %This loop should find intersections between the Nth curve and the N-1th curve
        N=NN(i); %Local number of complete waves
        Om=p./(N+e/(2*pi)); %Spindle speed, Hz
        OM(i,:)=Om; %Spindle speed for outputfrom loop
        
        % find intersection of curve with previous curve
        if i>1 %Don't plot first iteration as it theoretically extends to infinity
            s1=[blim',OM(i-1,:)']; %rows are observations, colomns are data type (x,y)
            s2=[blim',OM(i,:)'];%rows are observations, colomns are data type (x,y)
            [~,I(i-1)] = min(pdist2(s2,s1,'seuclidean','Smallest',1)); %This finds the index of the closest pair for the Nth curve. Use standardised euclidian because of orders of magnitude
            [~,I2(i)] = min(pdist2(s1,s2,'seuclidean','Smallest',1));%This finds closest x,y pair for N-1th curve
            %       plot(blim(I),Om(i-1,I),'rx',blim(I2),Om(i,I2),'bx','markersize', 20)%rows are observations, colomns are data type (x,y)
        end
    end
    
    % Values for output
    hold on
    for i=2:length(I)
        out(i).Om=OM(i,I2(i):I(i));
        out(i).blim=blim(I2(i):I(i));
        plot(OM(i,I2(i):I(i)),blim(I2(i):I(i)),'color',cmap(a,:)) %Plots this curve
        hold on
    end

end
%%
    ylabel('$\tilde{w}$','interpreter','latex')
    xlabel('$\tilde{\Omega}$','interpreter','latex')
    ylim([0 0.8])

colormap(cmap);
c=colorbar;
c.Label.String = '$R_A$';
c.Label.Interpreter = 'latex';
caxis([0 max(RA)])
