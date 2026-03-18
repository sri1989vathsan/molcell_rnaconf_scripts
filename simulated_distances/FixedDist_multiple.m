function [] = FixedDist_multiple(dis,values,numconf)
% dis = 200; 
% values = 1000;
t = linspace(0,dis,numconf+1);

% figure(11)
% scatter3(sqrt(x(1,:)),sqrt(x(2,:)),sqrt(x(3,:)),'.')
dis2d=[];
for i = 2:length(t)
x = [];    
x = randfixedsum(3,round(values+t(i)*100),t(i).^2,0,t(i).^2);
dis2d = [dis2d x(1,:)+x(2,:)];
end
[n1,xout1] = hist(sqrt(x(1,:)),100);
figure(1)
bar(xout1,n1*100/sum(n1))
xlabel('X value');
ylabel('Percentage of samples');

[n2,xout2] = hist(sqrt(x(2,:)),100);
figure(2)
bar(xout2,n2*100/sum(n2))
xlabel('Y value');
ylabel('Percentage of samples');

[n3,xout3] = hist(sqrt(x(3,:)),100);
figure(3)
bar(xout3,n3*100/sum(n3))
xlabel('Y value');
ylabel('Percentage of samples');

[n,xout] = hist(sqrt(dis2d),20);
figure(4)
bar(xout,n*100/sum(n))
xlabel('2D Distance');
ylabel('Percentage of samples');

end

function [x,v] = randfixedsum(n,m,s,a,b)

% Check the arguments.
if (m~=round(m))|(n~=round(n))|(m<0)|(n<1)
 error('n must be a whole number and m a non-negative integer.')
elseif (s<n*a)|(s>n*b)|(a>=b)
 error('Inequalities n*a <= s <= n*b and a < b must hold.')
end

% Rescale to a unit cube: 0 <= x(i) <= 1
s = (s-n*a)/(b-a);

% Construct the transition probability table, t.
% t(i,j) will be utilized only in the region where j <= i + 1.
k = max(min(floor(s),n-1),0); % Must have 0 <= k <= n-1
s = max(min(s,k+1),k); % Must have k <= s <= k+1
s1 = s - [k:-1:k-n+1]; % s1 & s2 will never be negative
s2 = [k+n:-1:k+1] - s;
w = zeros(n,n+1); w(1,2) = realmax; % Scale for full 'double' range
t = zeros(n-1,n);
tiny = 2^(-1074); % The smallest positive matlab 'double' no.
for i = 2:n
 tmp1 = w(i-1,2:i+1).*s1(1:i)/i;
 tmp2 = w(i-1,1:i).*s2(n-i+1:n)/i;
 w(i,2:i+1) = tmp1 + tmp2;
 tmp3 = w(i,2:i+1) + tiny; % In case tmp1 & tmp2 are both 0,
 tmp4 = (s2(n-i+1:n) > s1(1:i)); % then t is 0 on left & 1 on right
 t(i-1,1:i) = (tmp2./tmp3).*tmp4 + (1-tmp1./tmp3).*(~tmp4);
end

% Derive the polytope volume v from the appropriate
% element in the bottom row of w.
v = n^(3/2)*(w(n,k+2)/realmax)*(b-a)^(n-1);

% Now compute the matrix x.
x = zeros(n,m);
if m == 0, return, end % If m is zero, quit with x = []
rt = rand(n-1,m); % For random selection of simplex type
rs = rand(n-1,m); % For random location within a simplex
s = repmat(s,1,m);
j = repmat(k+1,1,m); % For indexing in the t table
sm = zeros(1,m); pr = ones(1,m); % Start with sum zero & product 1
for i = n-1:-1:1  % Work backwards in the t table
 e = (rt(n-i,:)<=t(i,j)); % Use rt to choose a transition
 sx = rs(n-i,:).^(1/i); % Use rs to compute next simplex coord.
 sm = sm + (1-sx).*pr.*s/(i+1); % Update sum
 pr = sx.*pr; % Update product
 x(n-i,:) = sm + pr.*e; % Calculate x using simplex coords.
 s = s - e; j = j - e; % Transition adjustment
end
x(n,:) = sm + pr.*s; % Compute the last x

% Randomly permute the order in the columns of x and rescale.
rp = rand(n,m); % Use rp to carry out a matrix 'randperm'
[ig,p] = sort(rp); % The values placed in ig are ignored
x = (b-a)*x(p+repmat([0:n:n*(m-1)],n,1))+a; % Permute & rescale x

end

