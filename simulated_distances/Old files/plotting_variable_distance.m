function [dis3d,dis2d] = simulatevariabledistance(maxdist,values)

% maxdist=300;
% 
% values = 1000;
x=randsphere(values,3,maxdist);

dis3d = (x(:,1).^2) + (x(:,2).^2) + (x(:,3).^2);

dis2d = (x(:,1).^2)+(x(:,2).^2);

[n1,xout1] = hist((abs(x(:,1))),100);
figure(1)
bar(xout1,n1*100/sum(n1))
xlabel('X value');
ylabel('Percentage of samples');

[n2,xout2] = hist(abs((x(:,2))),100);
figure(2)
bar(xout2,n2*100/sum(n2))
xlabel('Y value');
ylabel('Percentage of samples');

[n3,xout3] = hist(abs((x(:,3))),100);
figure(3)
bar(xout3,n3*100/sum(n3))
xlabel('Y value');
ylabel('Percentage of samples');

[n,xout] = hist(sqrt(dis2d),100);
figure(4)
bar(xout,n*100/sum(n))
xlabel('2D Distance');
ylabel('Percentage of samples');

[n4,xout4] = hist(sqrt(dis3d),100);
figure(5)
bar(xout4,n4*100/sum(n4))
xlabel('3D Distance');
ylabel('Percentage of samples');

end

function X = randsphere(m,n,r)
 
X = randn(m,n);
s2 = sum(X.^2,2);
X = X.*repmat(r*(gammainc(s2/2,n/2).^(1/n))./sqrt(s2),1,n);
 
end

% figure(10)
% scatter3(x(:,1),x(:,2),x(:,3),'.')