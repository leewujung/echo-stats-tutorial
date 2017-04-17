function theta = rand_piston_angle(N,bpa)

if isempty(bpa)
    theta_top = pi/2;
else
    theta_top = bpa/180*pi;
end

% position in the beam
count = 1;
theta = zeros(1,N);
while count <= N
    xx = rand(1);
    yy = rand(1);
    zz = rand(1);
    if sqrt(xx.^2+yy.^2+zz.^2)<1
        tmp = atan(sqrt(xx.^2+yy.^2)./zz);
        if tmp<=theta_top
            theta(count) = tmp;
            count = count+1;
        end
        %xy = sqrt(xx.^2+yy.^2);
        %theta(count) = atan(xy./(zz));
        %count = count +1;
    end
end
