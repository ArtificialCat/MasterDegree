% Fisher's linear discriminant. 
% : xi is column vector of which element is test metric. 
% Therefore size of row is the number of test metrics. 
% Number of column is the number of data sets. 
x1 = rand(2, 30) + 0.75.*ones(2,30); %[d1(:,c1) d1(:,c2)]';
x2 = rand(2, 30) + 0.3 .*ones(2,30); %[d2(:,c1) d2(:,c2)]';

m1 = mean(x1')'; 
m2 = mean(x2')'; 
m = m1 + m2; 
Sw1 = zeros(size(x1, 1), size(x1,1)); 
Sw2 = zeros(size(x1, 1), size(x1,1)); 
for i = 1:size(x1,2)
    Sw1 =  Sw1 + (x1(:,i)-m1)*(x1(:,i)-m1)';
end
for i = 1:size(x2,2)
    Sw2 =  Sw2 + (x2(:,i)-m2)*(x2(:,i)-m2)';
end

Sw = Sw1 + Sw2; 
w = Sw^(-1)*(m2-m1);
scatter(x1(1,:), x1(2,:), 10, 'ro');
hold on;
scatter(x2(1,:), x2(2,:), 10,'bo');
c = 0.5.*m; 
quiver(c(1,1), c(2,1), 1, -w(1,1)/w(2,1));
quiver(c(1,1), c(2,1), -1, w(1,1)/w(2,1));
%quiver(w(1,1),w(2,1), 0.5)
hold off;
figure; 
y1 = x1'*w; 
y2 = x2'*w; 
hist([y1 y2])