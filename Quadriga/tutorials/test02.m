numb = 1000;
test = zeros(3,8*numb+8);
for i = 1:numb+1
    test(:,8*(i-1)+1:8*(i-1)+8) = channel_sim();
end
scatter(test(1,:),test(2,:)); hold on;
top = [0 50; 0 50*sqrt(3)];
bottom = [1 50; 0 -50*sqrt(3)];
plot(top(1,:),top(2,:)); hold on;
plot(bottom(1,:),bottom(2,:));
hold off;